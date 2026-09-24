"""
FASTA Digestion Example
========================
Read proteins from a FASTA file with fastatacular and digest each one with trypsin.

peptacular does not read files. Its sequence functions accept any object with a
``sequence`` string attribute, so a fastatacular ``SequenceEntry`` goes straight in.
Install the reader with ``pip install fastatacular``.
"""

import os
import tempfile

from fastatacular import read_fasta

import peptacular as pt

FASTA_TEXT = """\
>sp|P1|EXAMPLE1 Example protein one
MPEPTIDEKAGVSEQR
>sp|P2|EXAMPLE2 Example protein two
MSEQKGARVTDEPTIDER
"""


def run():
    # Write the FASTA text to a real file, in a temporary directory that is removed on exit.
    with tempfile.TemporaryDirectory() as tmp_dir:
        fasta_path = os.path.join(tmp_dir, "example.fasta")
        with open(fasta_path, "w") as f:
            f.write(FASTA_TEXT)

        entries = read_fasta(fasta_path)

    for entry in entries:
        print(f"\n{entry.accession} {entry.pname} ({entry.sequence})")
        for peptide, span in pt.digest(entry, pt.Protease.TRYPSIN, missed_cleavages=1, min_len=4):
            print(f"  {peptide:<20} {span.start:>3}-{span.end:<3} mass={pt.mass(peptide):.4f}")


if __name__ == "__main__":
    run()
