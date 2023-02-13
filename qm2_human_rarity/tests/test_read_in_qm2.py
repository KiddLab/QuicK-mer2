from qm2_human_rarity import compare_against_1000
import numpy as np


def test_read_in_qm2_fake_sample():
    f = open("test_file.bed", 'w')
    for i in range(2300498):
        f.write("chro" + '\t' + "1" + '\t' + '2' + '\t' + '2.4' + '\n')
    f.close()
    test_file = "test_file.bed"
    output_table = compare_against_1000.read_in_qm2(test_file)
    test_answer = np.zeros((1, 2300498), dtype=np.single)
    test_answer.fill(2.4)
    assert (test_answer == output_table).all(), "Sample should be a np.table of length 2300498 and have 2.4 as all entries."

def test_read_in_qm2_real_sample():
    test_file = "../../../../DM09.qm2.CN.1k.bed"
    output_table = compare_against_1000.read_in_qm2(test_file)
