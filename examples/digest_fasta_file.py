"""
FASTA Digestion Example
========================
Parse protein sequences from a FASTA file and digest each one with trypsin.
"""

import os
import tempfile

import peptacular as pt

FASTA_TEXT = """\
>sp|P1|EXAMPLE1 Example protein one
MPEPTIDEKAGVSEQR
>sp|P2|EXAMPLE2 Example protein two
MSEQKGARVTDEPTIDER
"""


def run():
    # Write the FASTA text to a real file so pt.parse_fasta can read it from disk. Use a
    # temporary directory so the file is cleaned up automatically when the block exits.
    with tempfile.TemporaryDirectory() as tmp_dir:
        fasta_path = os.path.join(tmp_dir, "example.fasta")
        with open(fasta_path, "w") as f:
            f.write(FASTA_TEXT)

        records = pt.parse_fasta(fasta_path)

        for record in records:
            protein = pt.parse(record.sequence)
            print(f"\n{record.header} ({record.sequence})")
            for span in protein.digest(pt.Proteases.TRYPSIN, missed_cleavages=1, min_len=4):
                peptide = protein[span]
                print(f"  {peptide.serialize()}  mass={peptide.mass():.4f}")


if __name__ == "__main__":
    run()
