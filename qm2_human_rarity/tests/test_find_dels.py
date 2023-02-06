from qm2_human_rarity import compare_against_1000
import numpy as np


def test_find_dels():
    test_array = np.array([1.1, 1.2, 2, 1.3, 1.4, 1.45, 5, 1.3, 1.3, 1.4, 1.2, 3.1, 4.2])
    test_answer = [[[1.3, 3], [1.4, 4], [1.45, 5]], [[1.3, 7], [1.3, 8], [1.4, 9], [1.2, 10]]]
    dup_windows = compare_against_1000.find_deletions(test_array)
    assert test_answer == dup_windows, "The correct test array is [[[1.3, 3], [1.4, 4], [1.45, 5]], [1.3, 7], [1.3, 8], [1.4, 9], [1.2, 10]]."
