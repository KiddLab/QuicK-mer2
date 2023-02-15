from qm2_human_rarity import compare_against_1000
import numpy as np
import os.path


def test_read_in_qm2_fake_sample():
    # creates a test array with 2.4 as all entries
    if os.path.isfile("test_file.bed"):
        pass
    else:
        f = open("test_file.bed", 'w')
        for i in range(2300498):
            f.write("chro" + '\t' + "1" + '\t' + '2' + '\t' + '2.4' + '\n')
        f.close()
    test_file = "test_file.bed"
    output_table, output_dict = compare_against_1000.read_in_qm2(test_file)
    test_answer = np.zeros((1, 2300498), dtype=np.single)
    test_answer.fill(2.4)
    assert (test_answer == output_table).all(), "Sample should be a np.table of length 2300498 and have 2.4 as all entries."


def test_read_in_qm2_real_sample():
    # utilizes sample DM09 to test array reading
    test_file = "qm2_human_rarity/DM09.qm2.CN.1k.bed"
    output_table, output_dict = compare_against_1000.read_in_qm2(test_file)
    first_five = np.array([0.909472, 1.037988, 0.615303, 0.604325, 1.049525], dtype=np.single)
    last_five = np.array([2.553247, 1.825867, 6.551373, 2.797043, 8.206348], dtype=np.single)
    assert (output_table[0, :5] == first_five).all(), "First five lines should match."
    assert (output_table[0, -5:] == last_five).all(), "Last five lines should match."
