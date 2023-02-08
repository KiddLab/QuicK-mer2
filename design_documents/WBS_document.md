QuicK-mer2 Extension for Human Disease Application

- Activity 1 : Set up Github
  - [X] Task 1.1: Set up repo
  - [X] Task 1.2: Create high-level description
  - [X] Task 1.3: Create code structure 

- Activity 2: Build prototype 1

- Activity 2.1: Refine design to describe tasks
  - [ ] Task 2.1.1: Break project into activities
  - [ ] Task 2.2.2: Break project into tasks

- Activity 2.2: Test with QuicK-mer2 sample
  - [X] Task 2.2.1: Identify sample (DM09)
  - [X] Task 2.2.2: Download sample DM09 from the CHARGE patient dataset as a CRAM
  - [X] Task 2.2.3: Run the QuicK-mer2 pipeline on DM09, and retrieve resulting BED file (pre-processing step)

- Activity 2.3: Build database of 1000 genomes qm2 files (only done once, provided for future usage)
  - [ ] Task 2.3.1: Download all samples for the 1000 Genomes Project (see datasets.md)
  - [ ] Task 2.3.2: Process all samples through QuicK-mer2 and collect resulting BED files
  - [ ] Task 2.3.3: Build the comparison table, where rows = samples and columns = QuicK-mer2 windows based on GrCH38 (pre-defined in the original QuicK-mer2 GitHub) using Python 

- Activity 2.4: Define output for submission using Python code
  - [X] Task 2.4.1: Identify regions of deletions or duplications 

		- [X] The first implementation of this function will utilize a test matrix to identify at least three consecutive windows in the patient sample with either copy number < 1.5 (deletions) or > 2.5 (duplications) by iterating through each window and merging consecutive regions together until the condition is broke

		- [ ] The second implementation will scale up to use DM09, our test sample identified in Activity 2.2

		- [ ] The third implementation will allow reading in from a list of samples and identifying all deletions/duplications

  - [ ] Task 2.4.2: For each identified region of deletion or duplication, check for similar deletions/duplications in the 1000 Genomes Dataset table (Task 2.2.2.3)
  - [ ] Task 2.4.3: Identify all samples that also contain the deletion or duplication, and output a statitics of rarity (% of samples that also contain the region) 
  - [ ] Task 2.4.3.1: Possibility: identify any connected populations that may or may not contain the deletion or duplication 
  - [ ] Task 2.4.4: Produce first output: two BED files, one containing deletion regions and one containing duplication regions, with chromosome, start, end, median CN, and rarity %

- Activity 2.5: Save output
  - Task 2.5.1: [ ] Provide 1000 Genomes Table as a numpy array (.npy)
  - Task 2.5.2: [ ] Provide output patient processing as a BED file (.bed)

- Activity 3: Build onto CAVATICA 
