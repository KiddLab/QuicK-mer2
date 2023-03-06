from qm2_human_rarity import compare_against_1000
import pytest

test_file = "qm2_human_rarity/DM09_subset.qm2.bed"
output_table, output_dict = compare_against_1000.read_in_qm2(test_file, False)


def test_get_coords():
    return_coords = compare_against_1000.get_coords(0, 3, output_dict)
    return_answer = ["chr1", '0', '82642']
    assert return_coords == return_answer, "The correct return is [chr1, 0, 82642]."


def test_get_coords_flipped():
    with pytest.raises(Exception):
        compare_against_1000.get_coords(4, 3, output_dict)


def test_get_coords_no_chro_match():
    with pytest.raises(Exception):
        compare_against_1000.get_coords(0, 800000, output_dict)
