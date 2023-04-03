A tutorial for running the QuicK-mer2 Human Rarity Extension

-A sample file has been provided: DM09_subset.qm2.bed in the qm2_human_rarity folder

	-This is a human patient sample that has already been processed by QuicK-mer2, and copy number has been identified in relation to GrCH38

	- This sample has been trimmed down from the normal size for ease of processing - it includes only the first two chromosomes plus the beginning of chr3

Follow the tutorial_notebook.ipynb for additional details for processing DM09_subset!

1. Add the following import: from qm2_human_rarity import compare_against_1000

2. Read in the file using compare_against_1000.read_in_qm2
	- This needs two arguments: the filename, and a boolean for whether or not this is a full GrCH38 file (for this subset, this is False)
	- This produces two outputs: a numpy table of one row and X columns, where X = the number of windows in your file, and each entry is the copy number at that state; the second output is a dictionary where key = column index and value = the actual human coordinate

3. Calculate duplications: compare_against_1000.find_dups(output_table, output_dict)
	- The arguments are the numpy table and dictionary produced by step 2
	- The result is a list of lists of pairs, where each pair is [Copy number, window index], each inner list is a single duplication containing all the windows fulfilling the condition (cn > 2.5), and the final list is a compilation of all duplications found in the sample.

4. Calculate deletions: compare_against_1000.find_deletions(output_table, output_dict)
	- The arguments are the numpy table and dictionary produced by step 2
	- The result is a list of lists of pairs, where each pair is [Copy number, window index], each inner list is a single deletion containing all the windows fulfilling the condition (cn < 1.5), and the final list is a compilation of all deletions found in the sample.

5. Comparing duplications and deletions using compare_against_1000.compare_1000_genomes
	- The arguments are the numpy dictionary produced by step 2, a path to the provided 99_normalcy_range_tenk_genomes.npy, and the dups and dels provided by steps 3 and 4 (both optional arguments)
	- The 99_normalcy numpy table is an array where each column is a QuicK-mer2 window sized to GrCH38, and values are the copy number range that encompasses 99% of the copy number range seen in the 1000 Genomes dataset
	- The result is a dictionary of duplications and a dictionary of deletions, where the key is the variant and the value is True or False, with True indicating that this variation is rare, and outside of the boundary of 99% of the 1000 Genomes. 


6. Output duplications and deletions as a file using compare_against_1000.write_dups_and_dels(output_dict, dups, deletions, sample_name)
	- The output_dict is the same one from Step 2
	- Dups are produced by Step 3; this is optional, can do just deletions
	- Deletions are produced by Step 4; this is optional, can do just duplications
	- sample_name is an optional name, defaults to "sample"; can pass "DM09_subset" as the sample name
	- This function will write two files, sample_name_duplications.bed and sample_name_deletions.bed 