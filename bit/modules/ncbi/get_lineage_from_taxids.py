import subprocess
from bit.modules.ncbi.get_ncbi_tax_data import get_ncbi_tax_data
from bit.modules.general import report_failure


def get_lineage_from_taxids(taxids_file, output_file, include_strain=False):

    # ensure NCBI tax data is present (downloads if needed)
    get_ncbi_tax_data(quiet=True)

    if include_strain:
        pattern = r'{domain|superkingdom}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}\t{strain|subspecies|no rank}'
        header = "taxid\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\n"
    else:
        pattern = r'{domain|superkingdom}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}'
        header = "taxid\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"

    with open(taxids_file) as taxids_fh:
        p1 = subprocess.Popen(
            ["taxonkit", "lineage"],
            stdin=taxids_fh,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p2 = subprocess.Popen(
            ["taxonkit", "reformat2", "-r", "NA", "-f", pattern],
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p1.stdout.close()
        stdout, _ = p2.communicate()
        p1.wait()

    if p1.returncode != 0 or p2.returncode != 0:
        report_failure("taxonkit encountered an error. Is it installed and is the NCBI taxonomy data available?")

    with open(output_file, "w") as out:
        out.write(header)
        for line in stdout.decode().splitlines():
            if not line.strip():
                continue
            fields = line.split("\t")
            cut_line = "\t".join([fields[0]] + fields[2:])
            out.write(cut_line.replace(";", "\t") + "\n")
