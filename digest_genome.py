import argparse
import os
import pysam


# Function indexes FASTA file if no index is present
def indexFasta(
    fasta
):
    faidx = fasta + '.fai'
    if not os.path.isfile(faidx):
        pysam.faidx(fasta)


# Function split sequence into restriction fragments
def findFragments(
    sequence, motif
):
    # Convert sequence and motif to upper
    sequence = sequence.upper()
    motif = motif.upper()
    # Find all fragment starts
    searchStart = 0
    motifStarts = []
    while True:
        next = sequence.find(motif, searchStart)
        if next == -1:
            break
        else:
            motifStarts.append(next)
            searchStart = next + 1
    # Generate fragments
    motifLength = len(motif)
    fragStarts = [0] + [x + motifLength for x in motifStarts]
    fragEnds = motifStarts + [len(sequence)]
    fragments = zip(fragStarts, fragEnds)
    return(fragments)


# Creates interval trees containing restriction fragments in fasta file
def createFragmentBed(
    fasta, motif, bed
):
    # Open fasta file and get chromsome sequences
    index = 1
    with pysam.FastaFile(fasta) as infile, open(bed, 'wt') as outfile:
        for chrom in infile.references:
            sequence = infile.fetch(chrom)
            # Find fragments and save to file
            for start, end in findFragments(sequence, motif):
                line = '{}\t{}\t{}\t{}\n'.format(chrom, start, end, index)
                outfile.write(line)
                index += 1


if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='digest genome fasta file using motif'
    )
    parser.add_argument(
        '--fasta', required=True, type=str, help='genome FASTA file'
    )
    parser.add_argument(
        '--motif', required=True, type=str, help='RE motif used to digest genome'
    )
    parser.add_argument(
        '--bed', required=True, type=str, help='output bed file'
    )
    args = parser.parse_args()
    # Index fasta and create bed file
    indexFasta(args.fasta)
    fragments = createFragmentBed(args.fasta, args.motif, args.bed)
