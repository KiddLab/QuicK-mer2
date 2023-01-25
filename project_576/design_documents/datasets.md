# Brief description of datasets for the QuicK-mer2 extension

1. A sample QuicK-mer2 output bed file to test on (Sample DM09, belonging to GMKF Charge Data)

2. All samples involved in the 1000 Genomes Project, processed through QuicK-mer2
Each bed file will be processed into a final numpy table, where row = sample and column = window coordinates (pre-decided based on reference genome, hg38), and each entry is a sample's copy number at that window
This table will be created during this process and provided for future use