import numpy as np


def read_in_qm2(filename):
    """
    Read in a BED file produced by QuicK-mer2, and create a numpy array with one row and 2300498 windows, where each
    value is the copy number for that particular window
    Return: sample_table, the numpy array
    """
    # set variables
    human_window_count = 2300498  # the number of QM2 windows produced by a sample that has been processed to GrCH38
    sample_table = np.zeros((1, human_window_count), dtype=np.single)
    table_key = 0
    sample_dict = {}

    # open file
    f = open(filename, 'rt')
    inFile = f.readlines()

    # parse file
    for line in inFile:
        line = line.rstrip().split()
        coord = (line[0], line[1], line[2])
        sample_dict[coord] = float(line[3])  # key = coordinates in bed format (chro, start, end); value = copy number in float
        sample_table[1, table_key] = float(line[3])  # assigning the correct column with the copy number as a float
        table_key += 1
        if table_key == human_window_count:
            print("The windows of this sample exceeds the number of GrCH38 windows.")
            break
    f.close()

    return sample_table


def build_1000_genomes(file_list):
    pass


def find_deletions(sample_table):
    """
    :param sample_table: the read-in QM2 file converted in a numpy array
    :return: final_dels, a list which contains lists which contains pairs, where each pair is a copy number that is below
    1.5 and its corresponding window index, each list is a series of at least 3 consecutive windows that all have their
    copy number below 1.5, and the overall list that contains the series of deletions found in sample_table
    """
    # tested and implemented
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
    """
    :param sample_table: the read-in QM2 file converted in a numpy array
    :return: final_dups, a list which contains lists which contains pairs, where each pair is a copy number that is above
    2.5 and its corresponding window index, each list is a series of at least 3 consecutive windows that all have their
    copy number above 2.5, and the overall list that contains the series of duplications found in sample_table
    """
    # tested and implemented
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
