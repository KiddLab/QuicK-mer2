from qm2_human_rarity import compare_against_1000
import numpy as np

test_file = "qm2_human_rarity/DM09_subset.qm2.bed"
output_table, output_dict = compare_against_1000.read_in_qm2(test_file, False)
final_dups = compare_against_1000.find_dups(output_table, output_dict)
final_dels = compare_against_1000.find_deletions(output_table, output_dict)


def test_compare_1000_genomes():
    normal_file_path = "../99_normalcy_range_tenk_genomes.npy"
    normal_dups, normal_dels = compare_against_1000.compare_1000_genomes(output_dict, normal_file_path, final_dups, final_dels)
    assert normal_dups['chr1:123494757-123585956'] is True, "This duplication should be rare."
    assert normal_dels['chr1:6378151-6385794'] is True, "This deletion should be rare."
