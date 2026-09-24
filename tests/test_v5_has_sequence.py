"""Sequence functions accept any object with a str ``sequence`` attribute (a FASTA or PEFF entry)."""

from dataclasses import dataclass

import pytest

import peptacular as pt


@dataclass(frozen=True)
class Entry:
    """Stands in for fastatacular.SequenceEntry / a PEFF entry: peptacular must not import them."""

    identifier: str
    sequence: str


ENTRY = Entry("sp|P1|TEST", "MKPEPTIDERSTK")
SEQ = ENTRY.sequence

CALLS = {
    "mass": lambda s: pt.mass(s),
    "mz": lambda s: pt.mz(s, 2),
    "comp": lambda s: pt.comp(s),
    "digest": lambda s: pt.digest(s, "trypsin", missed_cleavages=1),
    "cleavage_sites": lambda s: pt.cleavage_sites(s, "trypsin"),
    "fragment": lambda s: [f.mz for f in pt.fragment(s, ["b", "y"], [1])],
    "fast_fragment": lambda s: pt.fast_fragment(s, ["b", "y"], [1]),
    "isotopic_distribution": lambda s: pt.isotopic_distribution(s),
    "sequence_length": lambda s: pt.sequence_length(s),
    "strip_mods": lambda s: pt.strip_mods(s),
    "reverse": lambda s: pt.reverse(s),
    "calc_property": lambda s: pt.calc_property(s, "hphob_kyte_doolittle"),
    "is_subsequence": lambda s: pt.is_subsequence("PEP", s),
    "coverage": lambda s: pt.coverage(s, ["PEP"]),
}


@pytest.mark.parametrize("name", sorted(CALLS))
def test_entry_matches_its_sequence_string(name):
    call = CALLS[name]
    assert call(ENTRY) == call(SEQ)


def test_protocol_is_runtime_checkable_and_exported():
    assert "HasSequence" in pt.__all__
    assert isinstance(ENTRY, pt.HasSequence)
    assert not isinstance(42, pt.HasSequence)


def test_list_of_entries_uses_the_batch_path():
    entries = [ENTRY, Entry("b", "PEPTIDE")]
    assert pt.mass(entries) == pt.mass([SEQ, "PEPTIDE"])


def test_batch_and_diagnose_accept_entries_and_keep_the_original_input():
    [result] = pt.batch("mass", [ENTRY])
    assert result.ok and result.input is ENTRY
    assert result.value == pt.mass(SEQ)
    assert pt.diagnose(Entry("x", "PEP[Oxidation]TIDE"), "mass") is None


def test_sequence_attribute_is_read_as_proforma():
    assert pt.mass(Entry("m", "PEM[Oxidation]TIDE")) == pt.mass("PEM[Oxidation]TIDE")


@pytest.mark.parametrize("bad", [object(), Entry("x", 5), 3.2])  # ty: ignore[invalid-argument-type]
def test_objects_without_a_str_sequence_are_rejected(bad):
    with pytest.raises(TypeError, match="sequence"):
        pt.mass(bad)
    [result] = pt.batch("mass", [bad], errors="collect")
    assert result.error.code == "invalid_input"


def test_real_fastatacular_entries_digest(tmp_path):
    fastatacular = pytest.importorskip("fastatacular")
    path = tmp_path / "p.fasta"
    path.write_text(">sp|P1|A_HUMAN Test\nMKPEPTIDERSTK\n")
    [entry] = fastatacular.read_fasta(path)
    assert pt.digest(entry, "trypsin") == pt.digest(SEQ, "trypsin")
