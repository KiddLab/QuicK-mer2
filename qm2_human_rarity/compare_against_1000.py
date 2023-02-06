import numpy as np

print("compare a sample's deletions and duplications against 1000 genomes for rarity")
# read in a QM2 bed file and convert it into a numpy table

# the number of QM2 windows produced by a sample that has been processed to GrCH38
human_window_count = 1000000

def read_in_qm2(filename):
    f = open(filename, 'rt')
    inFile = f.readlines()
    sample_table = np.zeros((1, human_window_count), dtype=np.single)
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


def build_1000_genomes(file_list):
    pass


def find_deletions(sample_table):
    dele = False
    final_dels = []
    current = []
    window_index = 0
    for cn in sample_table:
        if dele is True:
            if cn < 1.5:
                current.append([cn, window_index])
            else:
                if len(current) >= 3:
                    final_dels.append(current)
                current = []
                dele = False
        else:
            if cn < 1.5:
                dele = True
                current.append([cn, window_index])

        window_index += 1

    return final_dels


def find_dups(sample_table):
    dup = False
    final_dups = []
    current = []
    window_index = 0
    for cn in sample_table:
        if dup is True:
            if cn > 2.5:
                current.append([cn, window_index])
            else:
                if len(current) >= 3:
                    final_dups.append(current)
                current = []
                dup = False
        else:
            if cn > 2.5:
                dup = True
                current.append([cn, window_index])

        window_index += 1

    return final_dups

