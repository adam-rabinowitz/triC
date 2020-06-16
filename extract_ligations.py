import argparse
import collections
import intervaltree
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


# Creates interval trees containing restriction fragments in fasta file
def parse_fragments(
    bed
):
    # Create output variables
    trees = collections.defaultdict(intervaltree.IntervalTree)
    lengths = collections.OrderedDict()
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
                trees[chrom].add(interval)
            # Skip zero length fragments
            elif start == end:
                pass
            # Raise error for negative length fragments
            else:
                raise ValueError('negative length interval')
            # Adjust chromsome lengths
            lengths[chrom] = max(lengths.get(chrom, 0), end)
    return(trees, lengths)


# Function find intervals overlapping query sequence
def find_overlaps(
    fragments, chrom, start, end
):
    intervals = fragments[chrom].overlap(start, end)
    tuples = [x[2] for x in intervals]
    return(tuples)


# Function to demultiplex intervals
def parse_probes(
    bed, fragments, proximity
):
    # Create output variables
    bait = collections.namedtuple('bait', ['location', 'proximal'])
    baits = collections.OrderedDict()
    with open(bed, 'rt') as infile:
        for line in infile:
            # Extract probe location
            chrom, start, end, name = line.strip().split('\t')
            start, end = int(start), int(end)
            # Find probe fragment index
            probe_overlap = find_overlaps(
                fragments, chrom, start, end
            )
            if len(probe_overlap) != 1:
                raise ValueError("no unique bait for {}".format(name))
            probe_overlap = probe_overlap[0]
            # Find proximity overlaps
            proximity_start = probe_overlap.start - proximity
            proximity_end = probe_overlap.end + proximity
            proximity_overlaps = find_overlaps(
                fragments, chrom, proximity_start, proximity_end
            )
            proximity_indices = set([x.index for x in proximity_overlaps])
            # Add data to output
            baits[name] = bait(
                probe_overlap, proximity_indices
            )
    return(baits)


# Function to extract all reads from a sam file
def parse_sam(
    sam, pattern='^(.*?):RF(\\d+):(\\d+)$'
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('unmapped', 0), ('mapped', 0)
    ])
    reads = collections.defaultdict(list)
    # Create regx and pattern to save data
    regx = re.compile('^(.*?):RF(\\d+):(\\d+)$')
    with pysam.Samfile(sam) as infile:
        for read in infile:
            # Skip unmaped reads
            if read.is_unmapped:
                counter['unmapped'] += 1
            # Process unmapped reads
            else:
                counter['mapped'] += 1
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
    return(reads, counter)


# Function removes reads with a single or duplicate alignments
def remove_duplicates(
    reads
):
    # Gnerate log and output variables
    counter = collections.OrderedDict([
        ('duplicates', 0), ('unique', 0)
    ])
    filtered = collections.defaultdict(list)
    # Loop through reads and dtermine if location are unique
    previous_locations = set()
    for name in reads.keys():
        # Extracr read locations
        read_data = reads[name]
        read_data.sort()
        locations = '_'.join(['_'.join(map(str, x[1])) for x in read_data])
        # Remove singletons
        if locations in previous_locations:
            counter['duplicates'] += len(read_data)
        else:
            counter['unique'] += len(read_data)
            previous_locations.add(locations)
            filtered[name] = read_data
    return(filtered, counter)


# Function to map reads to genomic fragments
def map_ligations(
    reads, fragments
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('zero', 0), ('single', 0), ('multiple', 0)
    ])
    ligations = collections.defaultdict(list)
    # Generate regx for extracting counts
    for name, read_data in reads.items():
        for read in read_data:
            chrom, start, end = read[1][0:3]
            overlaps = find_overlaps(fragments, chrom, start, end)
            if len(overlaps) == 0:
                counter['zero'] += 1
            elif len(overlaps) > 1:
                counter['multiple'] += 1
            else:
                ligations[name].append(overlaps[0])
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
    demultiplex = collections.OrderedDict()
    # Loop through probes and add to log and output
    probe_indices = set()
    probe_names = {}
    for name in probes.keys():
        counter[name] = 0
        demultiplex[name] = {}
        index = probes[name].location.index
        probe_indices.add(index)
        probe_names[index] = name
    # Loop though intervals
    for name, fragments in ligations.items():
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
            demultiplex[probe_name][name] = fragments
    return(demultiplex, counter)


# Function to extract fragment ligations
def remove_proximal(
    ligations, proximal
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('proximal', 0), ('duplicate', 0), ('ligated', 0)
    ])
    distal = {}
    # Loop through dictionaries
    for name, fragments in ligations.items():
        # Extract unique fragment indices
        fragment_indices = set([x.index for x in fragments])
        fragment_dict = {x.index: x for x in fragments}
        duplicates = len(fragments) - len(fragment_indices)
        counter['duplicate'] += duplicates
        # Find distal
        proximal_indices = fragment_indices.intersection(proximal)
        distal_indices = fragment_indices.difference(proximal)
        counter['proximal'] += len(proximal_indices)
        # Process distal
        if len(distal_indices) > 0:
            counter['ligated'] += len(distal_indices)
            index_list = list(distal_indices)
            index_list.sort()
            # Conjoinn ligations and store
            fragment_list = [fragment_dict[x] for x in index_list]
            distal[name] = fragment_list
    return(distal, counter)


# Function to save ligations to file
def save_ligations(
    probe_ligations, path
):
    with open(path, 'wt') as outfile:
        for probe, ligations in probe_ligations.items():
            for name, fragments in ligations.items():
                indices = [str(x.index) for x in fragments]
                locations = ['{}:{}-{}'.format(*x[0:3]) for x in fragments]
                line = '{}\t{}\t{}\t{}\n'.format(
                    probe, name, '_'.join(locations), '_'.join(indices)
                )
                outfile.write(line)


# Generates bigwig file for ligations from a single capture
def generate_bigwig(
    ligations, lengths, bigwig
):
    # Extract header from chromsome lengths
    header = list(lengths.items())
    # Count each fragment
    fragment_counts = collections.defaultdict(int)
    for reads in ligations.values():
        for read in reads:
            fragment_counts[read] += 1
    # Extract fragments and sort
    fragments = list(fragment_counts.keys())
    fragments.sort(key=lambda x: x[3])
    # Extract data
    chroms = [x.chr for x in fragments]
    starts = [x.start for x in fragments]
    ends = [x.end for x in fragments]
    values = [float(fragment_counts[x]) for x in fragments]
    # Open bigwig, write and close
    bw = pyBigWig.open(bigwig, 'w')
    bw.addHeader(header)
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
    # Create interval tree
    fragments, chrom_lengths = parse_fragments(args.digest)
    probes = parse_probes(args.probes, fragments, args.proximal)
    # Create read dictionary
    reads, mapped_counter = parse_sam(args.sam)
    print_counter('Mapping', mapped_counter)
    # Remove duplicates
    unique_reads, duplicate_counter = remove_duplicates(reads)
    print_counter('\nDuplicates', duplicate_counter)
    # Find identity of ligated fragments
    ligations, overlap_counter = map_ligations(unique_reads, fragments)
    print_counter('\nFragments', overlap_counter)
    # Demultiplex data
    probe_ligations, probe_counter = demultiplex_ligations(ligations, probes)
    print_counter('\nDemultiplex', probe_counter)
    # Extract distal ligations
    distal_ligations = collections.OrderedDict()
    for probe in probe_ligations.keys():
        distal_ligations[probe], distal_counter = remove_proximal(
            probe_ligations[probe], probes[probe].proximal
        )
        print_counter('\n{}'.format(probe), distal_counter)
    # Save ligations to text file
    ligation_path = args.prefix + '.captured_ligations.txt'
    save_ligations(distal_ligations, ligation_path)
    # Save ligations to bigwig
    for probe in distal_ligations:
        path = '.'.join([args.prefix, probe, 'bw'])
        generate_bigwig(distal_ligations[probe], chrom_lengths, path)
