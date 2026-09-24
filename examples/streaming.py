"""Stream proteins, collect mass results, and inspect an input diagnostic."""

import io
from contextlib import closing

from fastatacular import FastaReader

import peptacular as pt


def run():
    fasta = io.StringIO(">first\nPEPTIDE\n>second\nMKR\n")
    # FastaReader yields entries lazily; iter_batch takes them as-is (anything with a .sequence str).
    with FastaReader(fasta) as proteins:
        with closing(pt.iter_batch("mass", proteins, batch_size=2, errors="collect")) as results:
            for result in results:
                assert result.ok
                print(result.index, result.input.identifier, round(result.value, 4))

    results = pt.batch("mass", ["PEPTIDE", "PEP[UnknownModification]TIDE"], errors="collect")
    assert results[0].ok
    assert results[1].error.code == "unresolved_modification"
    print(results[1].error.message)

    diagnostic = pt.diagnose("PEP[+42]TIDE", "comp")
    assert diagnostic.code == "unavailable_composition"
    print(diagnostic.message)


if __name__ == "__main__":
    run()
