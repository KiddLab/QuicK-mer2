from qm2_human_rarity import compare_against_1000
import numpy as np


def test_find_dups():
    test_array = np.array([1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 2, 1])
    test_answer = [[[3, 2], [4, 3], [5, 4], [6, 5]]]
    dup_windows = compare_against_1000.find_dups(test_array)
    assert test_answer == dup_windows, "The correct test array is [[[3, 2], [4, 3], [5, 4], [6, 5]]]."
