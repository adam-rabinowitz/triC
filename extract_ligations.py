import argparse
import collections
import intervaltree
import gzip
import os
import pyBigWig
import pysam
import re


# Function to print dictionary with title
def print_counter(
    title, d
):
    print(title)
    for key, value in d.items():
        print('  {}: {}'.format(key, str(value)))


# Creates interval trees containing restriction fragments in bed file
def create_fragment_trees(
    bed
):
    # Create output variables
    fragment_trees = collections.defaultdict(intervaltree.IntervalTree)
    fragment = collections.namedtuple(
        'fragment', ['chr', 'start', 'end', 'index']
    )
    # Loop through lines of input bed file
    with open(bed, 'rt') as infile:
        for line in infile:
            # Extract fragment data
            chrom, start, end, index = line.strip().split('\t')
            start, end, index = int(start), int(end), int(index)
            # Add fragment to intervaltree if positive length
            if start < end:
                interval = intervaltree.Interval(
                    start, end, fragment(chrom, start, end, index)
                )
                fragment_trees[chrom].add(interval)
            # Skip zero length fragments
            elif start == end:
                pass
            # Raise error for negative length fragments
            else:
                raise ValueError('negative length interval')
            # Adjust chromsome lengths
    return(fragment_trees)


# Function to extract chromosome lengths from bed file
def get_chromosome_lengths(
    bed
):
    # Create output variables
    lengths = collections.OrderedDict()
    # Loop through lines of bed file
    with open(bed, 'rt') as infile:
        for line in infile:
            # Extract fragment data
            chrom, start, end = line.strip().split('\t')[:3]
            end = int(end)
            # Store lengths
            try:
                lengths[chrom] = max(lengths[chrom], end)
            except KeyError:
                lengths[chrom] = end
    # Convert lengths to list of tuples and return
    lengths = list(lengths.items())
    return(lengths)


# Function find intervals overlapping query sequence
def find_overlaps(
    fragment_trees, chrom, start, end
):
    intervals = fragment_trees[chrom].overlap(start, end)
    fragments = [x[2] for x in intervals]
    return(fragments)


# Function to calculate genomc distance vetween 2 fragments
def calculate_distance(
    frag1, frag2
):
    # Calculate distances
    if frag1.chr != frag2.chr:
        distance = float('inf')
    else:
        distance = frag2.start - frag1.end
    # Assert distance is positive and return
    assert(distance > 0)
    return(distance)


# Function to demultiplex intervals
def parse_probes(
    bed, fragment_trees, proximity
):
    # Create output variables
    bait = collections.namedtuple('bait', ['fragment', 'proximal'])
    baits = collections.OrderedDict()
    with open(bed, 'rt') as infile:
        for line in infile:
            # Extract probe location
            chrom, start, end, name = line.strip().split('\t')
            start, end = int(start), int(end)
            # Find probe fragment index
            probe_fragments = find_overlaps(
                fragment_trees, chrom, start, end
            )
            if len(probe_fragments) != 1:
                raise ValueError("no unique bait for {}".format(name))
            probe_fragment = probe_fragments[0]
            # Find proximity overlaps
            proximity_start = probe_fragment.start - proximity
            proximity_end = probe_fragment.end + proximity
            proximity_fragments = find_overlaps(
                fragment_trees, chrom, proximity_start, proximity_end
            )
            proximity_indices = set([x.index for x in proximity_fragments])
            # Add data to output
            baits[name] = bait(
                probe_fragment, proximity_indices
            )
    return(baits)


# Function to extract all reads from a sam file
def parse_sam(
    sam, pattern='^(.*?):RF(\\d+):(\\d+)$'
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('unmapped', 0), ('singleton', 0), ('duplicates', 0), ('unique', 0)
    ])
    previous_locations = set()
    reads = collections.defaultdict(list)
    # Create regx and pattern to save data
    regx = re.compile(pattern)
    # Loop through reads in SAM file
    with pysam.Samfile(sam) as infile:
        for read in infile:
            # Skip unmapped reads
            if read.is_unmapped:
                counter['unmapped'] += 1
            # Process mapped reads
            else:
                # Get read name and chunk
                name = regx.match(read.query_name).groups()[0]
                chunk = regx.match(read.query_name).groups()[1:3]
                chunk = tuple(map(int, chunk))
                # Get read location
                location = (
                    read.reference_name, read.reference_start,
                    read.reference_end, '-' if read.is_reverse else '+'
                )
                # Add to dictionary
                reads[name].append((chunk, location))
    # Loop through reads and determine if location are unique
    for name in list(reads.keys()):
        # Extract read locations
        read_data = reads[name]
        read_data = sorted(read_data, key=lambda x: x[0])
        locations = [x[1] for x in read_data]
        locations_str = '_'.join(
            ['{}_{}_{}_{}'.format(*x) for x in locations]
        )
        # Count and store unique reads
        if len(locations) == 1:
            counter['singleton'] += 1
            del reads[name]
        elif locations_str in previous_locations:
            counter['duplicates'] += len(locations)
            del reads[name]
        else:
            counter['unique'] += len(locations)
            reads[name] = locations
            previous_locations.add(locations_str)
    return(reads, counter)


# Function to map reads to genomic fragments
def map_ligations(
    reads, fragment_trees
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('zero', 0), ('single', 0), ('multiple', 0)
    ])
    ligations = collections.defaultdict(list)
    # Generate regx for extracting counts
    for name, reads in reads.items():
        for read in reads:
            chrom, start, end = read[0:3]
            fragments = find_overlaps(fragment_trees, chrom, start, end)
            if len(fragments) == 0:
                counter['zero'] += 1
            elif len(fragments) > 1:
                counter['multiple'] += 1
            else:
                ligations[name].append(fragments[0])
                counter['single'] += 1
    return(ligations, counter)


# Function to split fragments by bait
def demultiplex_ligations(
    ligations, probes
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('no baits', 0), ('multiple baits', 0)
    ])
    captured = {}
    # Loop through probes and add to log and output
    probe_indices = set()
    probe_names = {}
    for name in probes:
        counter[name] = 0
        index = probes[name].fragment.index
        probe_indices.add(index)
        probe_names[index] = name
    # Loop though intervals
    for read_name, fragments in ligations.items():
        # Get intervals overlapping baits
        fragment_indices = set([x.index for x in fragments])
        common_indices = fragment_indices.intersection(probe_indices)
        # Skip intervals without baits
        if len(common_indices) == 0:
            counter['no baits'] += len(fragments)
        # Skip intervals with mutiple baits
        elif len(common_indices) > 1:
            counter['multiple baits'] += len(fragments)
        # Process unique baits
        else:
            # Get probe name
            probe_index = common_indices.pop()
            probe_name = probe_names[probe_index]
            counter[probe_name] += len(fragments)
            captured[(read_name, probe_name)] = fragments
    return(captured, counter)


# Function to extract fragment ligations
def remove_bait_proximal(
    ligations, probes
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('duplicate', 0), ('bait', 0), ('proximal', 0), ('distal', 0)
    ])
    distal = {}
    # Loop through dictionaries
    for (read_name, probe_name), fragments in ligations.items():
        # Get probe information
        bait_index = probes[probe_name].fragment.index
        proximal_indices = probes[probe_name].proximal
        # Create fragment dictionary
        fragment_dict = {x.index: x for x in fragments}
        fragment_indices = set(fragment_dict.keys())
        duplicates = len(fragments) - len(fragment_indices)
        counter['duplicate'] += duplicates
        # Remove bait and count
        fragment_indices.remove(bait_index)
        counter['bait'] += 1
        # Find distal
        proximal_set = fragment_indices.intersection(proximal_indices)
        distal_set = fragment_indices.difference(proximal_indices)
        counter['proximal'] += len(proximal_set)
        counter['distal'] += len(distal_set)
        # Idenitfy distal fragments and store
        if len(distal_set) > 0:
            distal_fragments = [fragment_dict[x] for x in distal_set]
            distal[(read_name, probe_name)] = distal_fragments
    return(distal, counter)


# Function to extract trimers from ligation file
def extract_trimers(
    ligations, proximity
):
    # Generate output variables
    counter = collections.OrderedDict([
        ('single', 0), ('proximal', 0), ('trimer', 0)
    ])
    trimers = collections.OrderedDict()
    # Open file and loop through lines
    for (read_name, probe_name), fragments in ligations.items():
        # Calculate distance
        fragments = sorted(fragments, key=lambda x: x.index)
        distances = [
            calculate_distance(*f) for f in zip(fragments[:-1], fragments[1:])
        ]
        distal = [d >= proximity for d in distances]
        # Generate output and return
        frag_no = len(fragments)
        if len(fragments) < 2:
            counter['single'] += frag_no
        elif sum(distal) < 1:
            counter['proximal'] += frag_no
        else:
            counter['trimer'] += frag_no
            trimers[(read_name, probe_name)] = fragments
    return(trimers, counter)


# Function to save ligations to file
def save_ligations(
    ligations, probes, path
):
    # Process probes
    probe_str = {}
    for probe_name in probes:
        probe_fragment = probes[probe_name].fragment
        probe_location = "{}:{}-{}".format(*probe_fragment[0:3])
        probe_index = probe_fragment.index
        probe_str[probe_name] = "{}\t{}\t{}".format(
            probe_name, probe_location, probe_index
        )
    # Create open function
    if path.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    # Open file and loop through ligations
    with open_func(path, 'wt') as outfile:
        for (read_name, probe_name), fragments in ligations.items():
            # Generate fragment string
            fragment_locations = [
                "{}:{}-{}".format(*f[0:3]) for f in fragments
            ]
            location_str = '_'.join(fragment_locations)
            # Generate index string
            index_str = '_'.join(map(str, [x.index for x in fragments]))
            # Create output line and write to file
            line = '{}\t{}\t{}\t{}\n'.format(
                read_name, probe_str[probe_name], location_str, index_str
            )
            outfile.write(line)


# Generates bigwig file for ligations from a single capture
def generate_bigwig(
    ligations, lengths, prefix
):
    # Count each fragment for each bait
    probe_counts = {}
    for (read_name, probe_name), fragments in ligations.items():
        for fragment in fragments:
            try:
                probe_counts[probe_name][fragment] += 1
            except KeyError:
                probe_counts[probe_name] = collections.defaultdict(int)
                probe_counts[probe_name][fragment] += 1
    # Extract fragments and sort
    for probe_name, fragment_counts in probe_counts.items():
        fragments = list(fragment_counts.keys())
        fragments.sort(key=lambda x: x.index)
        # Extract data
        chroms = [x.chr for x in fragments]
        starts = [x.start for x in fragments]
        ends = [x.end for x in fragments]
        values = [float(fragment_counts[x]) for x in fragments]
        # Create output path
        path = '.'.join([prefix, probe_name, 'bw'])
        # Open bigwig, write and close
        bw = pyBigWig.open(path, 'w')
        bw.addHeader(lengths)
        bw.addEntries(
            chroms, starts, ends=ends, values=values
        )
        bw.close()


if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='idenitfy ligation events from BAM/SAM file'
    )
    parser.add_argument(
        '--sam', required=True, type=str, help='BAM/SAM alignment file'
    )
    parser.add_argument(
        '--digest', required=True, type=str, help='genome digestion BED file'
    )
    parser.add_argument(
        '--probes', required=True, type=str, help='probe BED file'
    )
    parser.add_argument(
        '--proximal', required=True, type=int, help='proximity distance'
    )
    parser.add_argument(
        '--prefix', required=True, type=str, help='output prefix'
    )
    args = parser.parse_args()
    # Generate output directory
    abs_path = os.path.abspath(args.prefix)
    abs_dir = os.path.dirname(abs_path)
    if not os.path.isdir(abs_dir):
        os.makedirs(abs_dir)
    # Create interval trees and get chromosome lengths
    fragment_trees = create_fragment_trees(args.digest)
    chrom_lengths = get_chromosome_lengths(args.digest)
    # Get probe data
    probes = parse_probes(args.probes, fragment_trees, args.proximal)
    # Create read dictionary
    reads, mapped_counter = parse_sam(args.sam)
    print_counter('Alignment', mapped_counter)
    # Find identity of ligated fragments
    ligations, overlap_counter = map_ligations(reads, fragment_trees)
    print_counter('\nAssignment', overlap_counter)
    # Demultiplex data
    captured, probe_counter = demultiplex_ligations(ligations, probes)
    print_counter('\nDemultiplex', probe_counter)
    # Extract distal ligations
    distal, distal_counter = remove_bait_proximal(captured, probes)
    print_counter('\nDistal', distal_counter)
    # Extract trimers
    trimers, trimer_counter = extract_trimers(distal, args.proximal)
    print_counter('\nTrimers', trimer_counter)
    # Save ligations to file
    save_ligations(
        distal, probes, args.prefix + '.all_ligations.txt.gz'
    )
    save_ligations(
        trimers, probes, args.prefix + '.trimeric_ligations.txt.gz'
    )
    # Save bigwigs to file
    generate_bigwig(
        distal, chrom_lengths, args.prefix + '.all_ligations'
    )
    generate_bigwig(
        trimers, chrom_lengths, args.prefix + '.trimeric_ligations'
    )
