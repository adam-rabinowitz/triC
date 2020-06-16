# triC - scripts for the analysis of triC datasets

digest_genome.py:
    
    Script digests a genome, using a supplied motif, and generates a BED file
    of the digested fragments. The script takes three arguments:

    --fasta - Path to genome FASTA file
    --motif - Sequence recognised by restriction enzyme
    --bed - Path to output bed file

    The digested fragments do not include the motif that would norally be
    present on both ends of the fragment. The fourth column of the BED file is
    the name/index of the fragment.

exract_ligations.py:

    Script identifies ligation events assocaited with a set of supplied
    capture-c probes. The script takes five arguments:

    --sam - Path to iput SAM file. SAM file was generated as described in the
        following repository (https://github.com/oudelaar/TriC)
    --digest - Digested genome generated using the digest_genome.py script.
    --probes - Path to a 4 column BED file containing the location of the
        capture-c probes. The fourth column should be the name of the captured
        region.
    --proximal - The size of the upstream and downstream region to exclude from
        the ligation analysis.
    --prefix - The prefix of the output ligation files.

    An output text file of four columns is generated containing the ligations
    for all captured fragments. The four columns are, probe name, FASTQ read name,
    location of the ligated fragments and the index of the ligated fragments.
    
    Additionaly a bigwig is generated for each probes showing the count of captured
    fragments.
