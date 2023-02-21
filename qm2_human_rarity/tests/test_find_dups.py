from qm2_human_rarity import compare_against_1000
import numpy as np


def test_find_dups():
    test_array = np.array([1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 2, 1])
    test_chro_dict = {}
    for i in range(13):
        test_chro_dict[i] = ["chr1", 0, 1]
    test_answer = [[[3, 2], [4, 3], [5, 4], [6, 5]]]
    dup_windows = compare_against_1000.find_dups(test_array, test_chro_dict)
    assert test_answer == dup_windows, "The correct test array is [[[3, 2], [4, 3], [5, 4], [6, 5]]]."


def test_find_dups_chro_break_and_end():
    test_array = np.array([1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 1, 4, 5, 6])
    test_chro_dict = {}
    for i in range(17):
        if i < 10:
            test_chro_dict[i] = ["chr1", 0, 1]
        else:
            test_chro_dict[i] = ["chr2", 0, 1]
    test_answer = [[[3, 2], [4, 3], [5, 4], [6, 5]], [[4, 14], [5, 15], [6, 16]]]
    dup_windows = compare_against_1000.find_dups(test_array, test_chro_dict)
    assert test_answer == dup_windows, "The correct test array is [[[3, 2], [4, 3], [5, 4], [6, 5]], [[4, 14], [5, 15], [6, 16]]]."


def test_find_dups_none_found():
    test_array = np.array([1, 1, 1, 1, 1, 1])
    test_chro_dict = {}
    for i in range(13):
        test_chro_dict[i] = ["chr1", 0, 1]
    test_answer = "No duplications found in this sample."
    dup_windows = compare_against_1000.find_dups(test_array, test_chro_dict)
    assert test_answer == dup_windows, "There should be no duplications in this test."


def test_find_dups_empty():
    test_answer = "No table provided."
    test_chro_dict = {}
    for i in range(13):
        test_chro_dict[i] = ["chr1", 0, 1]
    dup_windows = compare_against_1000.find_dups(None, test_chro_dict)
    assert test_answer == dup_windows, "There should be no table provided in this test."
