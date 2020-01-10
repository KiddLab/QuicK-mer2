# QuicK-mer2
k-mer based analysis for paralog specific copy number estimation

To compile the software, clone the repository, move the QuicK-mer2 directory and type make.

You may receive some warnings. This will create an executable named quicKmer2.  


## Usage
The basic functionality of quickMer2 is described by executing the program with no options.

'''
./quicKmer2 
QuicK-mer2
Operation modes: 
	index	Index a bed format kmer list
	count	CNV estimate from library
	search	Search K-kmer in genome
	est	GC normalization into copy number
	sparse	Fractionate indexed kmer for memory reduction or regenerate GC control/Window

Simple operation:
1. Construct a dictionary from fasta using "search"
2. Count depth from sample fasta/fastq "count"
3. Estimate copy number with "est"
'''