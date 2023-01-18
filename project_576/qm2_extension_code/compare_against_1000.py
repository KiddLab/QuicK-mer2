import numpy
import numpy as np

print("compare a sample's deletions and duplications against 1000 genomes for rarity")
# read in a QM2 bed file and convert it into a numpy table

# the number of QM2 windows produced by a sample that has been processed to GrCH38
human_window_count = 1000000


def read_in_qm2(filename):
    f = open(filename, 'rt')
    inFile = f.readlines()
    sample_table = np.zeros((1, human_window_count))
    table_key = 0
    sample_dict = {}
    for line in inFile:
        line = line.rstrip().split()
        coord = (line[0], line[1], line[2])
        sample_dict[coord] = float(
            line[3])  # key = coordinates in bed format (chro, start, end); value = copy number in float
        sample_table[1, table_key] = float(line[3])  # assigning the correct column with the copy number as a float
        table_key += 1
    f.close()
