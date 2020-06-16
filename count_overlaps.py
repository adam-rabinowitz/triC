import argparse
import collections
from scipy.stats import hypergeom


# Function to print dictionary with title
def print_counter(
    title, d
):
    print(title)
    for key, value in d.items():
        if isinstance(key, tuple):
            key = '_'.join(key)
        print('  {}: {}'.format(key, str(value)))


# Function to determine if two fragments overlap
def is_overlapping(
    frag1, frag2
):
    # Check chromosome
    if frag1.chrom != frag2.chrom:
        overlap = False
    elif frag1.start >= frag2.end:
        overlap = False
    elif frag1.end <= frag2.start:
        overlap = False
    else:
        overlap = True
    return(overlap)


# Function to process fragment and target BED files
def process_fragments(
    fragment_bed, region_bed
):
    # Gnerate output variables
    frag = collections.namedtuple('frag', ['chrom', 'start', 'end'])
    fragments = collections.OrderedDict()
    targets = collections.OrderedDict()
    # Extract fragments
    with open(fragment_bed, 'rt') as infile:
        for line in infile:
            chrom, start, end, index = line.strip().split('\t')
            fragments[int(index)] = frag(chrom, int(start), int(end))
    # Find targets
    with open(region_bed, 'rt') as infile:
        for line in infile:
            chrom, start, end, name = line.strip().split('\t')
            targets[name] = set()
            target = frag(chrom, int(start), int(end))
            for index, fragment in fragments.items():
                if is_overlapping(fragment, target):
                    targets[name].add(index)
    # Return data
    return(fragments, targets)


# Function to extract loops from file. Function expects a path to a tab
# delimited file with three columns: bait name, region1 name region2 name
def parse_loops(
    loop_file
):
    loops = collections.OrderedDict()
    with open(loop_file) as infile:
        for line in infile:
            bait, region1, region2 = line.strip().split('\t')
            try:
                loops[bait].append((region1, region2))
            except KeyError:
                loops[bait] = [(region1, region2)]
    return(loops)


# Function to calculate genomc distance vetween 2 fragments
def calculate_distance(
    frag1, frag2
):
    # Calculate distances
    if frag1.chrom != frag2.chrom:
        distance = float('inf')
    else:
        distance = frag2.start - frag1.end
    # Assert distance is positive and return
    assert(distance > 0)
    return(distance)


# Function to demultiplex ligations
def demultiplex_ligations(
    ligations, fragments
):
    # Generate output variables
    counter = collections.OrderedDict()
    probes = collections.OrderedDict()
    # Loop through file and extract data
    with open(ligations, 'rt') as infile:
        for line in infile:
            # Extract index and fragment data
            bait, name, locations, indices = line.strip().split('\t')
            index_list = list(map(int, indices.split('_')))
            assert(index_list == sorted(index_list))
            fragment_list = [fragments[i] for i in index_list]
            # Add data to output
            try:
                counter[bait] += len(fragment_list)
            except KeyError:
                counter[bait] = len(fragment_list)
            try:
                probes[bait][name] = (index_list, fragment_list)
            except KeyError:
                probes[bait] = collections.OrderedDict()
                probes[bait][name] = (index_list, fragment_list)
    return(probes, counter)


# Function to extract trimers from ligation file
def extract_trimers(
    ligations, min_distance
):
    # Generate output variables
    counter = collections.OrderedDict([
        ('single', 0), ('proximal', 0), ('trimer', 0)
    ])
    trimers = collections.OrderedDict()
    # Open file and loop through lines
    for name, (indices, fragments) in ligations.items():
        # Calculate distance
        distances = [
            calculate_distance(*f) for f in zip(fragments[:-1], fragments[1:])
        ]
        distal = [d >= min_distance for d in distances]
        # Generate output and return
        if len(fragments) < 2:
            counter['single'] += 1
        elif sum(distal) < 1:
            counter['proximal'] += 1
        else:
            counter['trimer'] += 1
            trimers[name] = set(indices)
    return(trimers, counter)


# Find overlaps
def count_overlaps(
    trimers, indices1, indices2
):
    # Create output variables
    counter = collections.OrderedDict([
        ('total', 0), ('n1', 0), ('n2', 0), ('n12', 0)
    ])
    # Loop through trimers
    for ligations in trimers.values():
        counter['total'] += 1
        # Find overlaps
        if len(indices1.intersection(ligations)) > 0:
            counter['n1'] += 1
            if len(indices2.intersection(ligations)) > 0:
                counter['n2'] += 1
                counter['n12'] += 1
        elif len(indices2.intersection(ligations)) > 0:
            counter['n2'] += 1
    return(counter)


# Function to calculate hypergeometric pvalue for upper tail
def hypergeom_test_upper(
    total, n1, n2, n12
):
    n12 -= 1
    pvalue = hypergeom.sf(k=n12, M=total, n=n1, N=n2)
    return(pvalue)


# Function to calculate hypergeometric pvalue for lower tail
def hypergeom_test_lower(
    total, n1, n2, n12
):
    # Adjust successes
    pvalue = hypergeom.cdf(k=n12, M=total, n=n1, N=n2)
    return(pvalue)


# Function to calculate hypergeometric pvalue for lower tail
def hypergeom_median(
    total, n1, n2
):
    # Adjust successes
    median = hypergeom.median(M=total, n=n1, N=n2)
    return(median)


# Function to calculate pvalues
def hypergeom_test_twosided(
    total, n1, n2, n12
):
    # Calculate lower and upper tails and median
    lower = hypergeom_test_lower(
        total=total, n1=n1, n2=n2, n12=n12
    )
    upper = hypergeom_test_upper(
        total=total, n1=n1, n2=n2, n12=n12
    )
    median = hypergeom_median(
        total=total, n1=n1, n2=n2
    )
    # Create output and return
    hyper = collections.namedtuple('hyper', ['median', 'lower', 'upper'])
    return(hyper(median, min(1.0, 2*lower), min(1.0, 2*upper)))


# Calculate enrichemnt
def perform_independence_test(
    loops, regions, trimers, path
):
    # Generate header
    header = '\t'.join([
        'bait', 'region1', 'region2', 'total_ligations', 'region1_ligations',
        'region2_ligations', 'region12_ligations', 'median12_ligations',
        'lower_pvalue', 'upper_pvalue'
    ])
    # Open file and loop through counts
    with open(path, 'wt') as outfile:
        outfile.write(header + '\n')
        for probe, region_list in loops.items():
            for region1, region2 in region_list:
                counts = count_overlaps(
                    trimers[probe], regions[region1], regions[region2]
                )
                hyper_results = hypergeom_test_twosided(**counts)
                line_list = (
                    [probe, region1, region2] +
                    list(map(str, counts.values())) +
                    list(map(str, hyper_results))
                )
                outfile.write('\t'.join(line_list) + '\n')


if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='statistically examine three way loops'
    )
    parser.add_argument(
        '--ligations', required=True, help='ligation file'
    )
    parser.add_argument(
        '--digest', required=True, help='genome digestion BED file'
    )
    parser.add_argument(
        '--proximal', required=True, type=int, help='minimum ligation distance'
    )
    parser.add_argument(
        '--regions', required=True, help='region of interest BED file'
    )
    parser.add_argument(
        '--loops', required=True, help='file containing loop anchors'
    )
    parser.add_argument(
        '--output', required=True, help='output prefix'
    )
    args = parser.parse_args()
    # Read in fragments and find indices of regions of interest
    fragments, regions = process_fragments(
        args.digest, args.regions
    )
    loops = parse_loops(args.loops)
    # Read in ligations
    distal_ligations, ligation_counter = demultiplex_ligations(
        args.ligations, fragments
    )
    print_counter('Demultiplex', ligation_counter)
    # Get trimers
    trimers = collections.OrderedDict()
    for probe, ligations in distal_ligations.items():
        trimers[probe], trimer_counts = extract_trimers(ligations, args.proximal)
        print_counter('\n{}'.format(probe), trimer_counts)
    # Count overlaps
    perform_independence_test(
        loops, regions, trimers, args.output
    )
