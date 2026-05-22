from unittest.mock import patch, MagicMock
from bit.modules.ncbi.get_lineage_from_taxids import get_lineage_from_taxids


def test_get_lineage_from_taxids_with_strain(tmp_path):

    input_file = tmp_path / "taxids.txt"
    output_file = tmp_path / "output.tsv"
    input_file.write_text("562\n338191\n")

    # simulate what taxonkit reformat2 outputs:
    # col1=taxid, col2=raw ;-separated lineage (dropped by cut -f 1,3-), col3+=extracted fields
    mock_reformat2_output = (
        "562\tBacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli\t"
        "Bacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli\tNA\n"
        "338191\tArchaea;Nitrososphaerota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Nitrosopumilus\t"
        "Archaea\tNitrososphaerota\tNitrososphaeria\tNitrosopumilales\tNitrosopumilaceae\tNitrosopumilus\tNA\tNA\n"
    ).encode()

    mock_p1 = MagicMock()
    mock_p1.stdout = MagicMock()
    mock_p1.returncode = 0

    mock_p2 = MagicMock()
    mock_p2.communicate.return_value = (mock_reformat2_output, b"")
    mock_p2.returncode = 0

    with patch("bit.modules.ncbi.get_lineage_from_taxids.get_ncbi_tax_data"), \
         patch("bit.modules.ncbi.get_lineage_from_taxids.subprocess.Popen", side_effect=[mock_p1, mock_p2]):

        get_lineage_from_taxids(str(input_file), str(output_file), include_strain=True)

    lines = output_file.read_text().strip().split("\n")

    assert lines[0] == "taxid\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain"
    assert lines[1] == "562\tBacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli\tNA"
    assert lines[2] == "338191\tArchaea\tNitrososphaerota\tNitrososphaeria\tNitrosopumilales\tNitrosopumilaceae\tNitrosopumilus\tNA\tNA"


def test_get_lineage_from_taxids_no_strain(tmp_path):

    input_file = tmp_path / "taxids.txt"
    output_file = tmp_path / "output.tsv"
    input_file.write_text("562\n338191\n")

    mock_reformat2_output = (
        "562\tBacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli\t"
        "Bacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli\n"
        "338191\tArchaea;Nitrososphaerota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Nitrosopumilus\t"
        "Archaea\tNitrososphaerota\tNitrososphaeria\tNitrosopumilales\tNitrosopumilaceae\tNitrosopumilus\tNA\n"
    ).encode()

    mock_p1 = MagicMock()
    mock_p1.stdout = MagicMock()
    mock_p1.returncode = 0

    mock_p2 = MagicMock()
    mock_p2.communicate.return_value = (mock_reformat2_output, b"")
    mock_p2.returncode = 0

    with patch("bit.modules.ncbi.get_lineage_from_taxids.get_ncbi_tax_data"), \
         patch("bit.modules.ncbi.get_lineage_from_taxids.subprocess.Popen", side_effect=[mock_p1, mock_p2]):

        get_lineage_from_taxids(str(input_file), str(output_file), include_strain=False)

    lines = output_file.read_text().strip().split("\n")

    assert lines[0] == "taxid\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies"
    assert lines[1] == "562\tBacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli"
    assert lines[2] == "338191\tArchaea\tNitrososphaerota\tNitrososphaeria\tNitrosopumilales\tNitrosopumilaceae\tNitrosopumilus\tNA"
