from pathlib import Path
from bit.modules.lineage_to_tsv import get_rank, convert_lineage_to_tsv


def test_get_rank():
    input_lineage = "p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales"
    expected_remaining_lineage = "c__Gammaproteobacteria;o__Enterobacterales"
    expected_rank = "Proteobacteria"

    result_lineage, result_rank = get_rank(input_lineage, "p__")

    assert result_lineage == expected_remaining_lineage
    assert result_rank == expected_rank


def test_convert_lineage_to_tsv(tmp_path):
    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "output.tsv"

    input_file.write_text("seq1\troot;d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Sulfurovaceae;g__Sulfurovum;s__Sulfurovum lithotrophicum\n")

    convert_lineage_to_tsv(str(input_file), str(output_file), make_taxid=False)

    lines = output_file.read_text().strip().split("\n")
    assert lines[0] == "seq_ID\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies"
    assert lines[1] == "seq1\tBacteria\tCampylobacterota\tCampylobacteria\tCampylobacterales\tSulfurovaceae\tSulfurovum\tSulfurovum lithotrophicum"
