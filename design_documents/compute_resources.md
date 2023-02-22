# Computing Resources Needed for QuicK-mer2 extension

Using a text editor (pyCharm) to write the code; sample data downloaded from Kidd Lab access to Great Lakes


At most, Great Lakes or computer cluster access
--needed to run QuicK-mer2 on samples in the first place
----A useful minimum for running the parent QuicK-mer2 is approximately 7 CPUs with 6 GB per CPU (6 threads is ideal, per the Github's recommendation)
--extension of QuicK-mer2 should be low computing power - primarily Python script that does file parsing and list comparison
--Most computing power needed in the extension will be required to open and parse the 1000 Genomes numpy table, which will be a couple million columns


CAVATICA:
-in the second part of this project, we plan to deploy QuicK-mer2 onto Gabriella Miller Kids' First's own computing program for widespread usage
-this will not require any computational resources, due to deployment on their specific cloud resource

Computing time:
-comparison to the 1000 genomes project should only take a few seconds * number of samples needed for comparison
-ex: testing all 200 patients of a human cohort = approximately 2000 seconds

