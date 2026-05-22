from bit.modules.general import get_package_path
from bit.modules.filter_kofamscan_results import filter_kofamscan_results

test_ko_annots = get_package_path("tests/data/ko-annots.txt")


def test_filter_kofamscan_results(tmp_path):

    output_file = tmp_path / "output.tsv"

    filter_kofamscan_results(str(test_ko_annots), str(output_file))

    lines = output_file.read_text().strip().split("\n")

    assert lines[0] == "gene_ID\tKO_ID\tKO_function"
    assert lines[1] == "k119_6520_1\tK11959\turea transport system substrate-binding protein"
    assert lines[2] == "k119_19560_1\tNA\tNA"
