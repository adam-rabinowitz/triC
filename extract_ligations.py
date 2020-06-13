import argparse
import collections
import intervaltree
import pysam
import re


# Function to print dictionary with title
def printCounter(
    title, d
):
    print(title)
    for key, value in d.items():
        print('  {}: {}'.format(key, str(value)))


# Creates interval trees containing restriction fragments in fasta file
def createFragmentTrees(
    bed
):
    # Create output variables
    trees = collections.defaultdict(intervaltree.IntervalTree)
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
    return(trees)


# Function find intervals overlapping query sequence
def findOverlappingFragments(
    fragments, chrom, start, end
):
    intervals = fragments[chrom].overlap(start, end)
    tuples = [x[2] for x in intervals]
    return(tuples)


# Function to demultiplex intervals
def readProbes(
    bed, fragments, proximity
):
    # Create output variables
    bait = collections.namedtuple('bait', ['location', 'proximal'])
    baitDict = collections.OrderedDict()
    with open(bed, 'rt') as infile:
        for line in infile:
            # Extract probe location
            chrom, start, end, name = line.strip().split('\t')
            start, end = int(start), int(end)
            # Find probe fragment index
            probeOverlap = findOverlappingFragments(
                fragments, chrom, start, end
            )
            if len(probeOverlap) != 1:
                raise ValueError("no unique bait for {}".format(name))
            # Find proximity overlaps
            proximityStart = probeOverlap[0].start - proximity
            proximityEnd = probeOverlap[0].end + proximity
            proximityOverlaps = findOverlappingFragments(
                fragments, chrom, proximityStart, proximityEnd
            )
            proximityIndices = set([x.index for x in proximityOverlaps])
            # Add data to output
            baitDict[name] = bait(
                probeOverlap[0], proximityIndices
            )
    return(baitDict)


# Function to extract all reads from a sam file
def parseSamFile(
    sam, pattern='^(.*?):RF(\\d+):(\\d+)$'
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('unmapped', 0), ('mapped', 0)
    ])
    readDict = collections.defaultdict(list)
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
                readDict[name].append((chunk, location))
    return(readDict, counter)


# Function removes reads with a single or duplicate alignments
def removeDuplicates(
    readDict
):
    # Gnerate log and output variables
    counter = collections.OrderedDict([
        ('duplicates', 0), ('unique', 0)
    ])
    filterDict = collections.defaultdict(list)
    # Loop through reads and dtermine if location are unique
    previousLocations = set()
    for name in readDict.keys():
        # Extracr read locations
        reads = readDict[name]
        reads.sort()
        locations = '_'.join(['_'.join(map(str, x[1])) for x in reads])
        # Remove singletons
        if locations in previousLocations:
            counter['duplicates'] += len(reads)
        else:
            counter['unique'] += len(reads)
            previousLocations.add(locations)
            filterDict[name] = reads
    return(filterDict, counter)


# Function to map reads to fragments
def mapReadsToFragments(
    readDict, fragments
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('zero', 0), ('single', 0), ('multiple', 0)
    ])
    fragmentDict = collections.defaultdict(list)
    # Generate regx for extracting counts
    for name, reads in readDict.items():
        for read in reads:
            chrom, start, end = read[1][0:3]
            overlaps = findOverlappingFragments(fragments, chrom, start, end)
            if len(overlaps) == 0:
                counter['zero'] += 1
            elif len(overlaps) > 1:
                counter['multiple'] += 1
            else:
                fragmentDict[name].append(overlaps[0])
                counter['single'] += 1
    return(fragmentDict, counter)


# Function to split fragments by bait
def demultiplexAlignments(
    fragmentDict, probes
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('no baits', 0), ('multiple baits', 0)
    ])
    demultiDict = collections.OrderedDict()
    # Loop through probes and add to log and output
    probeIndices = set()
    probeNames = {}
    for name in probes.keys():
        counter[name] = 0
        demultiDict[name] = {}
        index = probes[name].location.index
        probeIndices.add(index)
        probeNames[index] = name
    # Loop though intervals
    for name, fragments in fragmentDict.items():
        # Get intervals overlapping baits
        fragmentIndices = set([x.index for x in fragments])
        commonIndices = fragmentIndices.intersection(probeIndices)
        # Skip intervals without baits
        if len(commonIndices) == 0:
            counter['no baits'] += len(fragments)
        # Skip intervals with mutiple baits
        elif len(commonIndices) > 1:
            counter['multiple baits'] += len(fragments)
        # Process unique baits
        else:
            # Get probe name
            probeIndex = commonIndices.pop()
            probeName = probeNames[probeIndex]
            counter[probeName] += len(fragments)
            demultiDict[probeName][name] = fragments
    return(demultiDict, counter)


# Function to extract fragment ligations
def extractDistalLigations(
    ligations, proximal
):
    # Generate log and output variables
    counter = collections.OrderedDict([
        ('proximal', 0), ('duplicate', 0), ('ligated', 0)
    ])
    distal = collections.OrderedDict()
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
def saveLigations(
    ligationDict, path
):
    with open(path, 'wt') as outfile:
        for name, fragments in ligationDict.items():
            indices = [str(x.index) for x in fragments]
            locations = ['{}:{}-{}'.format(*x[0:3]) for x in fragments]
            line = '{}\t{}\t{}\n'.format(
                name, '_'.join(locations), '_'.join(indices)
            )
            outfile.write(line)


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
    # Create interval tree
    fragments = createFragmentTrees(args.digest)
    probes = readProbes(args.probes, fragments, args.proximal)
    # Create read dictionary
    readDict, mappedCounter = parseSamFile(args.sam)
    printCounter('Mapping', mappedCounter)
    # Remove duplicates
    uniqueReadDict, duplicateCounter = removeDuplicates(readDict)
    printCounter('\nDuplicates', duplicateCounter)
    # Find location of reads
    fragmentDict, overlapCounter = mapReadsToFragments(uniqueReadDict, fragments)
    printCounter('\nFragments', overlapCounter)
    # Demultiplex data
    demultiDict, captureCounter = demultiplexAlignments(fragmentDict, probes)
    printCounter('\nDemultiplex', captureCounter)
    # Extract ligations
    ligationDict = collections.OrderedDict()
    for probe in demultiDict.keys():
        ligationDict[probe], ligationCounter = extractDistalLigations(
            demultiDict[probe], probes[probe].proximal
        )
        printCounter('\n{}'.format(probe), ligationCounter)
    # Print ligationś
    for probe in ligationDict:
        path = '.'.join([args.prefix, probe, 'txt'])
        saveLigations(ligationDict[probe], path)
