# Tutorial for running QuicK-mer2
This is a tutorial for running QuicK-mer2.  It is based on human data from the 1000 
Genomes Project.  For more information, see [our paper](https://www.mdpi.com/2073-4425/11/2/141):

```
Rapid, Paralog-Sensitive CNV Analysis of 2457 Human Genomes Using QuicK-mer2. Shen F, Kidd JM.
Genes. 2020 Jan 29;11(2). pii: E141. doi:10.3390/genes11020141.
```


## Step 1 Download and compile quicKmer2

Setup a directory to download and work from.  Clone the repository and compile the required core program components.
There may be some compile warnings.

```
git clone https://github.com/KiddLab/QuicK-mer2.git
cd QuicK-mer2/
make
```

The QuicK-mer2 directory needs to be in your path so that the correction script can find required
utilities. You can add the directory to your path temporarily using 

```
export PATH=$PWD:$PATH
```

For shared installation, refer to your cluster user guide.

## Step 2 Download prebuilt reference and k-mer index files

Pre-computed masking information is available for some commonly used reference genomes in https://kiddlabshare.med.umich.edu/QuicK-mer/QuicK-mer2-refs/

For this tutorial we will use a version of the human GRCh38 assembly.

Make a new directory and download and unzip the required files. First, change you directory to be outside of the 
QuicK-mer2 git repository. Then:

```
mkdir GRCh38
cd GRCh38

wget https://kiddlabshare.med.umich.edu/QuicK-mer/QuicK-mer2-refs/GRCh38/GRCh38_BSM.fa.bed.gz
wget https://kiddlabshare.med.umich.edu/QuicK-mer/QuicK-mer2-refs/GRCh38/GRCh38_BSM.fa.gz
wget https://kiddlabshare.med.umich.edu/QuicK-mer/QuicK-mer2-refs/GRCh38/GRCh38_BSM.fa.qgc.gz
wget https://kiddlabshare.med.umich.edu/QuicK-mer/QuicK-mer2-refs/GRCh38/GRCh38_BSM.fa.qm.gz
wget https://kiddlabshare.med.umich.edu/QuicK-mer/QuicK-mer2-refs/GRCh38/include-control.bed.gz
gunzip *gz
```

This directory now contains several files to use with the  pipeline.  The genome reference
is a version of GRCh38 that does not include alternative or HLA haplotypes.  Files were created
as described in the publication with k=30, e=2, d=100, and w=1000.  

`GRCh38_BSM.fa` is the genome reference sequence
`GRCh38_BSM.fa.bed` give the coordinates for 1000 kmer windows
`GRCh38_BSM.fa.qgc` is a file that gives information on local GC content
`GRCh38_BSM.fa.qm` is the kmer index file
`include-control.bed` defines the utilized control regions.

## Step 3 Download a 1000 Genome file for testing

We will test the pipeline using a 30X coverage sample from the 1000 genomes project.
We will make a new directory and download the sample.  First, we will move out of the 
GRCh38 directory then make a new directory to hold the sample data.

The sequences are saved in the CRAM format.  To decompress them you must also have
the exact same reference fasta used when making the CRAM. This differs slightly 
from the file used for k-mer analysis as it included decoy and alternative sequences.  We
will additionally download the required fasta file and index it.

```
mkdir sample-data
cd sample-data
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
```
Depending on your network speed, it may take >~30 minutes to download the required files.


## Step 4 Run quickmer count

In this step the kmer count for each kmer in the computed index will be determined in the cram file.
Reads (formatted as fasta) are extracted from the CRAM using samtools.  Options for the CRAM reference fasta
and omitting unneeded decompression are included.  The results are then piped to the QuicKmer count program
and saved into a output file.

For this example, the quickmer count program uses six threads and required sufficient RAM (~) to load
the index.  You may need to run this on a compute node on your cluster following your own local 
procedures. First we will make a directory for storing the output.

```
mkdir sample-out
```

At this point, your layout should look like the following:

```
[jmkidd@gl-login1 qm2-tutorial]$ ls -lh
total 224K
drwxrwxr-x 2 jmkidd kiddlab 172 May  5 12:41 GRCh38
drwxrwxr-x 3 jmkidd kiddlab 223 May  5 12:24 QuicK-mer2
drwxrwxr-x 2 jmkidd kiddlab 160 May  5 12:52 sample-data
drwxrwxr-x 2 jmkidd kiddlab   0 May  5 13:24 sample-out
[jmkidd@gl-login1 qm2-tutorial]$ 
[jmkidd@gl-login1 qm2-tutorial]$ ls -lh *
GRCh38:
total 67G
-rw-rw-r-- 1 jmkidd kiddlab 3.0G Jan 21 14:39 GRCh38_BSM.fa
-rw-rw-r-- 1 jmkidd kiddlab  99M Jan 21 14:39 GRCh38_BSM.fa.bed
-rw-rw-r-- 1 jmkidd kiddlab 4.3G Jan 21 14:47 GRCh38_BSM.fa.qgc
-rw-rw-r-- 1 jmkidd kiddlab  49G Jan 21 19:40 GRCh38_BSM.fa.qm
-rw-rw-r-- 1 jmkidd kiddlab 882K Jan 21 16:25 include-control.bed

QuicK-mer2:
total 992K
-rw-rw-r-- 1 jmkidd kiddlab 1.5K May  5 12:24 lowess.py
-rw-rw-r-- 1 jmkidd kiddlab   96 May  5 12:24 makefile
-rwxrwxr-x 1 jmkidd kiddlab 117K May  5 12:24 quicKmer2
-rw-rw-r-- 1 jmkidd kiddlab  42K May  5 12:24 QuicKmer.c
-rw-rw-r-- 1 jmkidd kiddlab 3.9K May  5 12:24 README.md
-rwxrwxr-x 1 jmkidd kiddlab 2.1K May  5 12:24 smooth_GC_mrsfast.py

sample-data:
total 22G
-rw-rw-r-- 1 jmkidd kiddlab 3.1G May  5 12:49 GRCh38_full_analysis_set_plus_decoy_hla.fa
-rw-rw-r-- 1 jmkidd kiddlab 151K May  5 12:52 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
-rw-rw-r-- 1 jmkidd kiddlab  15G May  5 13:11 NA12878.final.cram

sample-out:
total 0

```

Then run the following command:

```
samtools view -F 3840 -T sample-data/GRCh38_full_analysis_set_plus_decoy_hla.fa --input-fmt-option required_fields=0x202 sample-data/NA12878.final.cram | \
awk '{print ">\n"$10}' | quicKmer2 count -t 6 GRCh38/GRCh38_BSM.fa /dev/fd/0 sample-out/NA12878.qm2 ;
```

This command used 6 threads for quicKmer2 count.  I run it on a compute node that has 7 available CPU (6 for quicKmer2 and 1 for samtools). Sufficient
RAM should be available to load the index files and other associated data, ~50 GB should be sufficient. The command should complete in ~25 minutes
and  will print the following information to standard out:

```
[Option] Set 6 threads
Hash Size: 0x100000000
First location: 0xC9ACBD81
Read 0x100000000 hash
Read 1G kmers
Read 2G kmers
...
Read 76G kmers
Counting elapse 690 sec, total 81843317281 kmers
Pileup finish
Read chain file 4294967296 entries
Mean sequencing depth: 25.23
Exit quicK-mer2 count
```

When completed, it will create two new files in sample-out/:
`NA12878.qm2.bin` this is the kmer count file
`NA12878.qm2.txt` is a summary of the depth x GC content profile


## Step 5 Run quickmer estimate

Next, perform GC correction and convert raw kmer counts to copy-number estimates in 1k windows.

```
quicKmer2 est GRCh38/GRCh38_BSM.fa sample-out/NA12878.qm2 sample-out/NA12878.qm2.CN.1k.bed
```

This makes the file `NA12878.qm2.png`, which is a graph of the depth by GC content and calculated correction curve
as well as `NA12878.qm2.CN.1k.bed`, which is the estimated depth in the indicated windows.

## Step 6 Prepare output for display in the browser
To display in the browser, a version without decoy and EBV sequences should be generated.  This can be done
from the unix command line:

```
grep -v decoy sample-out/NA12878.qm2.CN.1k.bed | grep -v chrEBV > sample-out/NA12878.qm2.CN.1k.bed.browser
```

The file `sample-out/NA12878.qm2.CN.1k.bed.browser` is now actually in bedGraph format.  It can be converted
to BigWig format (using bedGraphToBigWig) for easy display as a custom UCSC track.


To create a heatmap view file, the `make-colortrack-fordisplay.py` script can be used.  

```
make-colortrack-fordisplay.py --cn sample-out/NA12878.qm2.CN.1k.bed.browser --name NA12878
```

This makes the file NA12878.qm2.CN.1k.bed.browser.bedColor which is a heatmap view of the QuicKmer-2 results as a bed9 files which can be converted 
to bigBed format for easy display.




