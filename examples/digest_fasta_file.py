"""
FASTA Digestion Example
========================
Parse protein sequences from a FASTA file and digest each one with trypsin.
"""

import tempfile

import peptacular as pt

FASTA_TEXT = """\
>sp|P1|EXAMPLE1 Example protein one
MPEPTIDEKAGVSEQR
>sp|P2|EXAMPLE2 Example protein two
MSEQKGARVTDEPTIDER
"""


def run():
    # Write the FASTA text to a real file so pt.parse_fasta can read it from disk
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(FASTA_TEXT)
        fasta_path = f.name

    records = pt.parse_fasta(fasta_path)

    for record in records:
        protein = pt.parse(record.sequence)
        print(f"\n{record.header} ({record.sequence})")
        for span in protein.digest(pt.Proteases.TRYPSIN, missed_cleavages=1, min_len=4):
            peptide = protein[span]
            print(f"  {peptide.serialize()}  mass={peptide.mass():.4f}")


if __name__ == "__main__":
    run()
