from qm2_human_rarity import compare_against_1000
import numpy as np
import os.path
import pytest


def test_read_in_qm2_real_subset():
    # utilizes a subset of sample DM09 to test array reading
    test_file = "qm2_human_rarity/DM09_subset.qm2.bed"
    output_table, output_dict = compare_against_1000.read_in_qm2(test_file, False)
    first_one = 0.909472
    last_one = 1.929796
    first_table = round(float(output_table[0]), 6)
    last_table = round(float(output_table[-1]), 6)
    assert first_table == first_one, "First cn should match."
    assert last_table == last_one, "Last cn should match."
