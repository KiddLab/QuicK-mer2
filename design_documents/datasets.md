# Brief description of datasets for the QuicK-mer2 extension

1. Testing data: A sample QuicK-mer2 output bed file to test on (Sample DM09, belonging to GMKF Charge Data, private patient data not linked here)
The output bedfile contains four columns: chromosome, start coordinate, end coordinate, copy number
Validation: already ran through QuicK-mer2, and test file was originally sampled in a collaborative project (https://pubmed.ncbi.nlm.nih.gov/29300383/), but data is not public

2. Required comparison data: All samples involved in the 1000 Genomes Project, processed through QuicK-mer2
Each bed file will be processed into a final numpy table, where row = sample and column = window coordinates (pre-decided based on reference genome, hg38), and each entry is a sample's copy number at that window
This table will be created during this process and provided for future use

Samples for the 1000 Genomes Project can be found here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/, whose parent host is the International Genome Sample Resource (IGSR), the host of the 1000 Genomes Project (https://www.internationalgenome.org/data/). 