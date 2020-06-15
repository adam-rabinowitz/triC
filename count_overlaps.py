import collections
import numpy as np


# Function to read genomic fragments from a bed file
def get_fragments(
    bed
):  
    # Create output variables
    fragment = collections.namedtuple('fragment', ['chrom', 'start', 'end'])
    fragments = {}
    # Loop through bed file
    with open(bed, 'rt') as infile:
        for line in infile:
            chrom, start, end, index = line.strip().split()
            fragments[int(index)] = fragment(
                chrom, int(start), int(end)
            )
    return(fragments)


# Function to calculate genomc distance vetween 2 fragments
def calculate_distance(
    frag1, frag2
):
    # Calculate distances
    if frag1.chrom != frag2.chrom:
        distance = np.Inf
    else:
        distance = frag2.start - frag1.end
    # Assert distance is positive and return
    assert(distance > 0)
    return(distance)


# Function to extract trimers from ligation file
def generate_trimers(
    ligations, fragments, min_distance
):  
    # Open file and loop through lines
    with open('test.scyl.txt', 'rt') as infile:
        for line in infile:
            # Extract indices
            indices = line.strip().split('\t')[2]
            index_list = list(map(int, indices.split('_')))
            fragment_list = [fragments[x] for x in index_list]
            # Calculate distances and distal ligations
            distances = []
            for frags in zip(fragment_list[:-1], fragment_list[1:]):
                distances.append(calculate_distance(*frags))
            distal = [x >= min_distance for x in distances]
            # Generate output and return
            if len(fragment_list) < 2:
                output = ('single', None)
            elif sum(distal) < 1:
                output = ('proximal', None)
            else:
                output = ('trimer', (index_list, fragment_list))
            yield(output)


fragments = get_fragments('dm6.dpnII.bed')
min_distance = 1000
counter = collections.Counter()
for call, data in generate_trimers(
    'test.scyl.txt', fragments, 1000
):
    counter[call] += 1
print(counter)