import argparse
import collections
import gzip
import pyBigWig
import re


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


def read_trimers(
    path
):
    # Create output variable
    trimers = {}
    # Generate read function
    if path.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    # Open file and loop through entries
    with open_func(path, 'rt') as infile:
        for line_no, line in enumerate(infile):
            trimer_data = line.strip().split('\t')
            bait_location = trimer_data[2]
            prey_locations = trimer_data[4].split(',')
            all_locations = [bait_location] + prey_locations
            all_locations = [
                re.split(':|-', x) for x in all_locations
            ]
            all_locations = [
                (x[0], int(x[1]), int(x[2])) for x in all_locations
            ]
            trimers[line_no] = all_locations
    return(trimers)


# Function to extract regions for analsysis from text file
def parse_regions(
    path
):
    # Create output variable
    regions = []
    # Create regular expression
    regx = re.compile('(chr.*?):(\\d+)-(\\d+)')
    with open(path) as infile:
        for line in infile:
            region1, region2, label = line.strip().split('\t')
            # Process regions
            region1 = regx.match(region1).groups()
            region1 = (region1[0], int(region1[1]), int(region1[2]))
            region2 = regx.match(region2).groups()
            region2 = (region2[0], int(region2[1]), int(region2[2]))
            # Store data
            regions.append((region1, region2, label))
    return(regions)


# Function to find overlaps between chromosome intervals
def find_overlaps(
    query, subjects
):
    indices = set()
    for index, subject in enumerate(subjects):
        if query[0] != subject[0]:
            continue
        if query[1] >= subject[2]:
            continue
        if query[2] <= subject[1]:
            continue
        indices.add(index)
    return(indices)


# Find trimers overlapping two specified regions
def find_trimer_fragments(
    region1, region2, trimers
):
    # Create output variable
    all_fragments = []
    counts = collections.OrderedDict([
        ('region1', 0), ('region2', 0), ('both', 0)
    ])
    # Find overlaps
    for trimer in trimers.values():
        # Check for overlaps with region1
        region1_overlaps = find_overlaps(region1, trimer)
        region1_match = bool(region1_overlaps)
        region2_overlaps = find_overlaps(region2, trimer)
        region2_match = bool(region2_overlaps)
        # Count matches
        if region1_match:
            counts['region1'] += 1
        if region2_match:
            counts['region2'] += 1
        if region1_match & region2_match:
            counts['both'] += 1
        else:
            continue
        # Extract additional fragments to both overlapping regions
        overlap_indices = region1_overlaps.union(region2_overlaps)
        fragments = [
            fragment
            for index, fragment in enumerate(trimer)
            if index not in overlap_indices
        ]
        all_fragments.extend(fragments)
    return(all_fragments, counts)


# Function to generate bigwig files
def create_bed(
    region1, region2, label, path
):
    # Count each fragment for each bait
    line1 = '{}\t{}\t{}\t{}.region1\n'.format(
        *region1, label
    )
    line2 = '{}\t{}\t{}\t{}.region2\n'.format(
        *region2, label
    )
    with open(path, 'wt') as outfile:
        outfile.write(line1)
        outfile.write(line2)


# Function to generate bigwig files
def create_bigwig(
    fragments, chrom_lengths, path
):
    # Count each fragment for each bait
    fragment_counts = collections.defaultdict(int)
    for fragment in fragments:
        fragment_counts[fragment] += 1
    # Extract unique fragments and sort
    chroms = [x[0] for x in chrom_lengths]
    unique_fragments = list(fragment_counts.keys())
    sorted_fragments = sorted(
        unique_fragments, key=lambda x: (chroms.index(x[0]), x[1])
    )
    # Extract data
    chroms = [x[0] for x in sorted_fragments]
    starts = [x[1] for x in sorted_fragments]
    ends = [x[2] for x in sorted_fragments]
    values = [float(fragment_counts[x]) for x in sorted_fragments]
    # Open bigwig, write and close
    bw = pyBigWig.open(path, 'w')
    bw.addHeader(chrom_lengths)
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
        '--trimers', required=True, type=str, help='file containing trimers'
    )
    parser.add_argument(
        '--digest', required=True, type=str, help='genome digestion BED file'
    )
    parser.add_argument(
        '--regions', required=True, type=str, help='file containing trimer anchors'
    )
    parser.add_argument(
        '--prefix', required=True, type=str, help='output prefix'
    )
    args = parser.parse_args()
    # Extract chromosome lengths, trimers and region anchors
    chrom_lengths = get_chromosome_lengths(args.digest)
    trimers = read_trimers(args.trimers)
    regions = parse_regions(args.regions)
    # Loop through regions and extract bigwigs
    for region1, region2, label in regions:
        # Extract fragments
        fragments, counts = find_trimer_fragments(
            region1, region2, trimers
        )
        print(
            '{}:\n\tregion1: {}\n\tregion2: {}\n\tboth:    {}'.format(
                label, *counts.values()
            )
        )
        # Create bed file
        bed_path = '.'.join([args.prefix, label, 'bed'])
        create_bed(
            region1, region2, label, bed_path
        )
        # Create bigwig
        bw_path = '.'.join([args.prefix, label, 'bw'])
        create_bigwig(
            fragments, chrom_lengths, bw_path
        )
