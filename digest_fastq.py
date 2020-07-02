import argparse
import collections
import gzip


# Function to generate fastq files
def fastq_generator(
        path
):
    # Generate fastq tuple
    read = collections.namedtuple('read', ['name', 'sequence', 'quality'])
    # Create open function
    if path.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    # Create variable to parse file
    fastq, count = [], 0
    # Open file and loop through
    with open_func(path, 'rt') as infile:
        for line in infile:
            count += 1
            if count == 3:
                continue
            fastq.append(line.strip())
            if count == 4:
                yield(read(*fastq))
                fastq, count = [], 0
    if count != 0:
        raise ValueError('incomplete fastq records')


# Function split sequence into restriction fragments
def find_fragments(
    sequence, motif
):
    # Find sequence length and set search parameters
    sequence_length = len(sequence)
    search_start = 1
    search_end = sequence_length - 1
    # Find start locations of fra
    fragment_starts = [0]
    while True:
        next = sequence.find(motif, search_start, search_end)
        if next == -1:
            break
        else:
            fragment_starts.append(next)
            search_start = next + 1
    # Find fragment ends
    motif_length = len(motif)
    fragment_ends = [start + motif_length for start in fragment_starts[1:]]
    fragment_ends.append(sequence_length)
    # Create list containing paired starts and end
    fragments = zip(fragment_starts, fragment_ends)
    return(fragments)


# Function to split fastq read based on occurence of motif
def split_read(
    read, motif, min_length
):
    # Process read name
    old_name = read.name.split(' ', 1)
    read_number = old_name[1][0]
    new_name_prefix = old_name[0] + ':RF' + read_number
    assert(new_name_prefix.endswith(('1', '2')))
    # Find fragments and loop through them
    fragment_list = []
    for count, (start, end) in enumerate(
        find_fragments(read.sequence, motif), 1
    ):
        # Skip short fragments
        fragment_length = end - start
        if fragment_length < min_length:
            continue
        # Generate fastq data from fragment
        new_name = new_name_prefix + ':' + str(count)
        new_sequence = read.sequence[start:end]
        new_quality = read.quality[start:end]
        fragment = read._replace(
            name=new_name, sequence=new_sequence, quality=new_quality
        )
        fragment_list.append(fragment)
    return(fragment_list)


# Create open function
def split_fastq(
    infastq_list, outfastq, motif, min_length
):
    # Convert motif to uper to match fastq
    motif = motif.upper()
    # Select function to open fastq file
    if outfastq.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    # Open output file and loop through reads in input files
    with open_func(outfastq, 'wt') as outfile:
        for infastq in infastq_list:
            for read in fastq_generator(infastq):
                # Split reads into fragments and save
                fragments = split_read(read, motif, min_length)
                for fragment in fragments:
                    out_string = '{}\n{}\n+\n{}\n'.format(*fragment)
                    outfile.write(out_string)


# Perform digestion
if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='digest multiple input fastq files'
    )
    parser.add_argument(
        '--infastq', required=True, nargs='+', help='input fastq files'
    )
    parser.add_argument(
        '--motif', required=True, help='restriction enzyme motif'
    )
    parser.add_argument(
        '--minlength', required=True, type=int, help='minimum fragment length'
    )
    parser.add_argument(
        '--outfastq', required=True, help='output fastq file'
    )
    args = parser.parse_args()
    # Digest fastq files
    split_fastq(
        infastq_list=args.infastq, outfastq=args.outfastq, motif=args.motif,
        min_length=args.minlength
    )
