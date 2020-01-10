# QuicK-mer2
k-mer based analysis for paralog specific copy number estimation

To compile the software, clone the repository, move the QuicK-mer2 directory and type make.

You may receive some warnings. This will create an executable named quicKmer2.  


## Usage
The basic functionality of quickMer2 is described by executing the program with no options.

```
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
```

The typical flow is to first create a set of unique kmers from a genome reference fasta sequence using
the search command. The search command has several options

```
./quicKmer2 search

quicKmer2 search [Options] ref.fa

Options:
-h		Show this help information
-k [num]	Size of K-mer. Must be between 3-32. Default 30
-t [num]	Number of threads for edit distance search
-s [num]	Size of hash dictionary. Can use suffix G,M,K
-e [num]	Edit distance search. 0, 1 or 2
-d [num]	Edit distance depth threshold to keep. 1..255, Default 100
-w [num]	Output window definition size. Default 1000
-c [filename]	Input bedfile for GC control regions
```

The -c option specifies a bedfile of regions that are not expected to be copy-number variable
across the analyzed samples. Sex chromosomes, unplaced chromosome sequences, known segmental duplications
and known copy number varints should be excluded.  Exclusion files and be converted to inclusion files
using appropriate bedtools commands.


