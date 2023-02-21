from qm2_human_rarity import compare_against_1000
import numpy as np


def test_find_dels():
    test_array = np.array([1.1, 1.2, 2, 1.3, 1.4, 1.45, 5, 1.3, 1.3, 1.4, 1.2, 3.1, 4.2])
    test_chro_dict = {}
    for i in range(13):
        test_chro_dict[i] = ["chr1", 0, 1]
    test_answer = [[[1.3, 3], [1.4, 4], [1.45, 5]], [[1.3, 7], [1.3, 8], [1.4, 9], [1.2, 10]]]
    del_windows = compare_against_1000.find_deletions(test_array, test_chro_dict)
    assert test_answer == del_windows, "The correct test array is [[[1.3, 3], [1.4, 4], [1.45, 5]], [1.3, 7], [1.3, 8], [1.4, 9], [1.2, 10]]."

def test_find_dels_chro_break():
    test_array = np.array([1.1, 1.2, 2, 1.3, 1.4, 1.45, 5, 1.3, 1.3, 1.4, 1.2, 3.1, 4.2])
    test_chro_dict = {}
    for i in range(13):
        if i < 9:
            test_chro_dict[i] = ["chr1", 0, 1]
        else:
            test_chro_dict[i] = ["chr2", 0, 1]
    test_answer = [[[1.3, 3], [1.4, 4], [1.45, 5]]]
    del_windows = compare_against_1000.find_deletions(test_array, test_chro_dict)
    assert test_answer == del_windows, "The correct test array is [[[1.3, 3], [1.4, 4], [1.45, 5]] due to chromosome break."

def test_find_dels_end_of_file(): #add implement for end of file deletion
    test_array = np.array([1.1, 1.2, 2, 1.3, 1.4, 1.45, 5, 1.3, 1.3, 1.4, 1.2, 1.3])
    test_chro_dict = {}
    for i in range(13):
        if i < 9:
            test_chro_dict[i] = ["chr1", 0, 1]
        else:
            test_chro_dict[i] = ["chr2", 0, 1]
    test_answer = [[[1.3, 3], [1.4, 4], [1.45, 5]], [[1.4, 9], [1.2, 10], [1.3, 11]]]
    del_windows = compare_against_1000.find_deletions(test_array, test_chro_dict)
    assert test_answer == del_windows, "The correct test array is [[[1.3, 3], [1.4, 4], [1.45, 5]], [[1.4, 9], [1.2, 10], [1.3, 11]]] due to end of file deletion."

def test_find_dels_none_found():
    test_array = np.array([2, 2, 2, 2, 2, 2])
    test_chro_dict = {}
    for i in range(13):
        test_chro_dict[i] = "chr1"
    test_answer = "No deletions found in this sample."
    del_windows = compare_against_1000.find_deletions(test_array, test_chro_dict)
    assert test_answer == del_windows, "There should be no deletions in this test."


def test_find_dels_empty():
    test_answer = "No table provided."
    test_chro_dict = {}
    for i in range(13):
        test_chro_dict[i] = "chr1"
    del_windows = compare_against_1000.find_deletions(None, test_chro_dict)
    assert test_answer == del_windows, "There should be no table provided in this test."
