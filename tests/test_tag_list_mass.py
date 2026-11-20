"""The mass of a pipe-separated tag list comes from the first tag that carries a mass.

ProForma 2.1 Final, section 7.5 (Support for the joint representation of experimental data and
its interpretation), on what a pipe-separated list may hold:

    "Other combinations between CV/ontology names, accession numbers, and delta masses using
    synonyms are allowed, though they MUST be synonymous terms."

and on who decides which tag a reader uses:

    "It is however acknowledged that readers can choose to implement the parsing in different
    ways. Some tools may always take CV terms, others could take delta masses, and so on."

Section 6.6 (The information tag) shows that a group label and an information tag may be written
in either order inside one bracket:

    ELV[+42#g1|INFO:newly discovered]I[#g1|INFO:awesome]S
    [+42#g1]?ELV[INFO:newly discovered|#g1]I[INFO:awesome|#g1]S

Note what 6.6 does not show: in both examples the bracket that begins with ``INFO:`` holds no
mass-bearing tag at all, the ``+42`` being a separate modification. So 6.6 is evidence that a
massless tag may be written first, not that the specification writes one in front of a
mass-bearing tag. The licence for reading the mass from a later tag is 7.5's sentence above.

peptacular skips the tags that carry no mass of their own (an information tag, a group or
position label such as ``#g1`` with its optional score, and the placement tags of a static
modification rule, ``Position:``, ``Limit:``, ``CoMKP`` and ``CoMUP``, each of which refuses to
be weighed) and takes the first one that has a mass. Which tag wins among the tags that do carry
a mass is a reader's choice and is unchanged: the first of those, in written order. Tags are
never summed.

Expected values from pyteomics 5.0.1, run outside peptacular on 2026-11-20:

    >>> pyteomics.mass.calculate_mass(sequence="PEPTMIDE")                       # 930.4004489397
    >>> pyteomics.mass.calculate_mass(sequence="PEPTMIDE", average=True)         # 931.0199184788603
    >>> pyteomics.mass.calculate_mass(sequence="PEPTIDE")                        # 799.35996402671
    >>> pyteomics.mass.nist_mass["O"][16][0]                                     # 15.99491461956
    >>> pyteomics.mass.calculate_mass(sequence="PEPTM", ion_type="b", charge=1)  # 556.24356063384
    >>> pyteomics.mass.calculate_mass(sequence="IDE", ion_type="y", charge=1)    # 376.1714412394

Oxidation is the addition of one oxygen atom, so the oxygen mass above is the expected delta for
``[Oxidation]`` and twice it for ``[Formula:O2]``.
"""

import pytest

import peptacular as pt
from peptacular.diagnostics import CompositionError, PeptacularError
from peptacular.proforma_components import ModificationTags

PEPTMIDE = 930.400449
PEPTMIDE_AVERAGE = 931.019918
PEPTIDE = 799.359964
OXYGEN = 15.994915
B5_PEPTM = 556.243561
Y3_IDE = 376.171441

# Each row is (the tag list, the delta of the tag that should win, that tag written on its own).

# Lists whose first tag carries no mass. Every one of these reads 0.0 on main (the mass of the
# information tag or the group label), so each row fails there and passes here.
FIRST_TAG_CARRIES_NO_MASS = [
    ("[INFO:a|+10]", 10.0, "[+10]"),
    ("[INFO:a|Oxidation]", OXYGEN, "[Oxidation]"),
    ("[#g1|+10]", 10.0, "[+10]"),
    ("[INFO:a|#g1|+10]", 10.0, "[+10]"),
    ("[INFO:a|Formula:O2]", 2 * OXYGEN, "[Formula:O2]"),
    ("[#g1|INFO:a|Obs:+15.995]", 15.995, "[Obs:+15.995]"),
]

# Lists whose first tag does carry a mass. Which mass-bearing tag wins is not changed by this
# fix, so every row here reads the same on main and here.
FIRST_TAG_CARRIES_A_MASS = [
    ("[+10|INFO:a]", 10.0, "[+10]"),
    ("[Oxidation|+20]", OXYGEN, "[Oxidation]"),
    ("[+20|Oxidation]", 20.0, "[+20]"),
    ("[+10|+20]", 10.0, "[+10]"),
    ("[Obs:+15.995|Oxidation]", 15.995, "[Obs:+15.995]"),
    ("[Oxidation|Obs:+15.995]", OXYGEN, "[Oxidation]"),
    ("[Formula:O2|Oxidation]", 2 * OXYGEN, "[Formula:O2]"),
]

# No tag carries a mass, so the list weighs nothing, as it always has, and the first tag stands.
NO_TAG_CARRIES_A_MASS = [("[INFO:a|INFO:b]", 0.0, None), ("[INFO:a|#g1]", 0.0, None)]

ALL_ROWS = FIRST_TAG_CARRIES_NO_MASS + FIRST_TAG_CARRIES_A_MASS + NO_TAG_CARRIES_A_MASS


@pytest.mark.parametrize(("tags", "delta", "winner"), ALL_ROWS)
def test_mass_comes_from_the_first_mass_bearing_tag(tags: str, delta: float, winner: str | None) -> None:
    """pyteomics' PEPTMIDE plus the delta the winning tag states."""
    assert pt.mass(f"PEPTM{tags}IDE") == pytest.approx(PEPTMIDE + delta, abs=1e-5)


@pytest.mark.parametrize(("tags", "delta", "winner"), ALL_ROWS)
def test_fragment_ions_come_from_the_first_mass_bearing_tag(tags: str, delta: float, winner: str | None) -> None:
    """b5 holds the modified M, y3 does not, so only b5 moves. Both values are pyteomics'."""
    seq = f"PEPTM{tags}IDE"
    b_ions = pt.fragment(seq, ion_types="b", charges=1)
    y_ions = pt.fragment(seq, ion_types="y", charges=1)
    assert b_ions[4].mz == pytest.approx(B5_PEPTM + delta, abs=1e-5)
    assert y_ions[2].mz == pytest.approx(Y3_IDE, abs=1e-5)


@pytest.mark.parametrize(("tags", "delta", "winner"), [row for row in ALL_ROWS if row[2] in (None, "[+10]", "[+20]", "[Obs:+15.995]")])
def test_average_mass_of_a_delta_mass_tag(tags: str, delta: float, winner: str | None) -> None:
    """A delta mass is the same number in both mass modes, so pyteomics' average PEPTMIDE serves."""
    assert pt.mass(f"PEPTM{tags}IDE", monoisotopic=False) == pytest.approx(PEPTMIDE_AVERAGE + delta, abs=1e-3)


@pytest.mark.parametrize(("tags", "delta", "winner"), ALL_ROWS)
def test_every_mass_entry_point_reads_the_winning_tag(tags: str, delta: float, winner: str | None) -> None:
    """A list gives the same numbers as its winning tag written on its own."""
    pytest.importorskip("numpy")
    import numpy as np

    seq = f"PEPTM{tags}IDE"
    plain = f"PEPTM{winner}IDE" if winner else "PEPTMIDE"
    assert pt.mass(seq) == pytest.approx(pt.mass(plain), abs=1e-9)
    assert pt.mass(seq, monoisotopic=False) == pytest.approx(pt.mass(plain, monoisotopic=False), abs=1e-9)
    assert pt.mz(seq + "/2") == pytest.approx(pt.mz(plain + "/2"), abs=1e-9)
    assert [f.mz for f in pt.fragment(seq)] == pytest.approx([f.mz for f in pt.fragment(plain)], abs=1e-9)
    fast, fast_plain = pt.fast_fragment(seq), pt.fast_fragment(plain)
    assert sorted(fast) == sorted(fast_plain)
    for key, mzs in fast.items():
        assert mzs == pytest.approx(fast_plain[key], abs=1e-9)
    assert np.allclose(pt.fragment_arrays([seq])["mz"], pt.fragment_arrays([plain])["mz"], atol=1e-9)


def _outcome(call) -> object:
    """The call's value, or the name of the error it raises, so two spellings can be compared."""
    try:
        return call()
    except Exception as exc:  # noqa: BLE001 - the error type is the observation
        return type(exc).__name__


@pytest.mark.parametrize(("tags", "delta", "winner"), ALL_ROWS)
def test_every_composition_entry_point_reads_the_winning_tag(tags: str, delta: float, winner: str | None) -> None:
    """The composition paths too, including the refusal a delta mass earns on both spellings.

    A delta-mass tag has no composition, so ``comp`` and ``isotopic_distribution`` raise
    ``CompositionError`` for it. The point here is that the list and its winning tag get the same
    answer, refusal included: before this fix ``[INFO:a|+10]`` returned the unmodified
    composition and an envelope 10 Da light.
    """
    seq = f"PEPTM{tags}IDE"
    plain = f"PEPTM{winner}IDE" if winner else "PEPTMIDE"
    assert _outcome(lambda: pt.comp(seq)) == _outcome(lambda: pt.comp(plain))
    assert _outcome(lambda: round(pt.mass(seq, calculate_with_composition=True), 9)) == _outcome(
        lambda: round(pt.mass(plain, calculate_with_composition=True), 9)
    )
    assert _outcome(lambda: [round(d.mass, 9) for d in pt.isotopic_distribution(seq, max_isotopes=3)]) == _outcome(
        lambda: [round(d.mass, 9) for d in pt.isotopic_distribution(plain, max_isotopes=3)]
    )


@pytest.mark.parametrize(
    "seq",
    [
        "[INFO:a|+10]-PEPTIDE",
        "PEPTIDE-[INFO:a|+10]",
        "{INFO:a|+10}PEPTIDE",
        "[INFO:a|+10]?PEPTIDE",
    ],
)
def test_every_placement_of_a_modification_reads_the_list_the_same_way(seq: str) -> None:
    """Terminal, labile and unknown-position modifications, all 10 Da over pyteomics' PEPTIDE."""
    assert pt.mass(seq) == pytest.approx(PEPTIDE + 10.0, abs=1e-5)


def test_a_static_modification_reads_the_list_the_same_way() -> None:
    assert pt.mass("<[INFO:a|+10]@M>PEPTMIDE") == pytest.approx(PEPTMIDE + 10.0, abs=1e-5)


def test_an_interval_reads_the_list_the_same_way() -> None:
    """Interval modifications have a mass and an m/z. fragment and fast_fragment refuse them."""
    assert pt.mass("PEP(TM)[INFO:a|+10]IDE") == pytest.approx(PEPTMIDE + 10.0, abs=1e-5)
    assert pt.mz("PEP(TM)[INFO:a|+10]IDE/2") == pytest.approx(pt.mz("PEP(TM)[+10]IDE/2"), abs=1e-9)


def test_a_bare_x_is_covered_by_a_mass_bearing_tag_anywhere_in_the_list() -> None:
    """check_unknown_residues reads the whole list, so the order of the tags cannot matter."""
    assert pt.mass("PEPX[INFO:a|+10]TIDE") == pytest.approx(PEPTIDE + 10.0, abs=1e-5)
    assert pt.mass("PEPX[#g1|+10]TIDE") == pytest.approx(PEPTIDE + 10.0, abs=1e-5)
    assert pt.mass("PEPX[INFO:a|+10]TIDE") == pytest.approx(pt.mass("PEPX[+10|INFO:a]TIDE"), abs=1e-9)
    assert len(pt.fragment("PEPX[INFO:a|+10]TIDE")) == len(pt.fragment("PEPX[+10|INFO:a]TIDE"))
    # A list of massless tags still leaves the X bare.
    with pytest.raises(PeptacularError, match="X has no mass of its own"):
        pt.mass("PEPX[INFO:a|#g1]TIDE")


def test_composition_comes_from_the_same_tag() -> None:
    assert pt.comp("PEPTM[INFO:a|Oxidation]IDE") == pt.comp("PEPTM[Oxidation]IDE")
    assert pt.comp("PEPTM[INFO:a|Formula:O2]IDE") == pt.comp("PEPTM[Formula:O2]IDE")
    # A delta mass has no composition, whichever position it is written in. Before this fix
    # ``[INFO:a|+10]`` returned the unmodified composition instead, 10 Da light.
    with pytest.raises(CompositionError):
        pt.comp("PEPTM[INFO:a|+10]IDE")
    with pytest.raises(CompositionError):
        pt.comp("PEPTM[+10|INFO:a]IDE")


def test_charge_comes_from_the_same_tag() -> None:
    """ProForma 2.1 section 11.5's example, with an information tag written in front of it."""
    with_info = pt.parse("PEPT[INFO:a|Formula:Zn:z+2]IDE/2")
    plain = pt.parse("PEPT[Formula:Zn:z+2]IDE/2")
    assert with_info.mass() == pytest.approx(plain.mass(), abs=1e-9)
    assert with_info.mz() == pytest.approx(plain.mz(), abs=1e-9)
    assert round(with_info.mass() / with_info.mz()) == 4


def test_a_glycan_tag_behind_an_information_tag() -> None:
    assert pt.mass("PEPTN[INFO:a|Glycan:HexNAc1]IDE") == pytest.approx(pt.mass("PEPTN[Glycan:HexNAc1]IDE"), abs=1e-9)


def test_serialisation_keeps_the_written_order() -> None:
    """Reading the mass from a later tag must not reorder or drop anything."""
    for seq in ["PEPTM[INFO:a|+10]IDE", "PEPTM[+10|INFO:a]IDE", "[INFO:a|#g1|+10]?PEPTIDE", "PEPTM[#g1|INFO:a|Obs:+15.995]IDE"]:
        assert pt.parse(seq).serialize() == seq


@pytest.mark.parametrize(("tags", "delta", "winner"), ALL_ROWS)
def test_mass_tag_and_has_mass_tag(tags: str, delta: float, winner: str | None) -> None:
    """The property the rule lives in, read directly."""
    mod = ModificationTags.from_string(tags.strip("[]"))
    assert mod.has_mass_tag is (winner is not None)
    if mod.has_mass_tag:
        assert mod.mass_tag.get_mass() == pytest.approx(delta, abs=1e-5)
        assert mod.get_mass() == pytest.approx(delta, abs=1e-5)
    else:
        # No tag carries a mass, so the first tag is returned and the list weighs nothing.
        assert mod.mass_tag is mod.first_tag
        assert mod.get_mass() == 0.0


def test_mass_tag_is_the_first_mass_bearing_tag_not_the_first_tag() -> None:
    mod = ModificationTags.from_string("INFO:a|#g1|+10")
    assert mod.first_tag is not mod.mass_tag
    assert mod.first_tag.serialize() == "INFO:a"
    assert mod.mass_tag.serialize() == "+10"
    assert len(mod) == 3


PLACEMENT_TAGS = ["Position:M", "Limit:2", "CoMKP", "CoMUP"]


@pytest.mark.parametrize("placement", PLACEMENT_TAGS)
@pytest.mark.parametrize("where", ["first", "middle", "last"])
def test_a_placement_tag_carries_no_mass_wherever_it_is_written(placement: str, where: str) -> None:
    """``Position:``, ``Limit:``, ``CoMKP`` and ``CoMUP`` all refuse to be weighed.

    Only ``Position:`` was skipped at first, so ``[INFO:a|Limit:2|+10]`` raised "PlacementTag has
    no mass" on the mass path and returned the unmodified value on the composition path.
    """
    tags = {
        "first": f"[{placement}|INFO:a|+10]",
        "middle": f"[INFO:a|{placement}|+10]",
        "last": f"[INFO:a|+10|{placement}]",
    }[where]
    seq = f"PEPTM{tags}IDE"
    assert pt.mass(seq) == pytest.approx(PEPTMIDE + 10.0, abs=1e-5)
    assert pt.mass(seq, monoisotopic=False) == pytest.approx(PEPTMIDE_AVERAGE + 10.0, abs=1e-3)
    assert pt.fragment(seq, ion_types="b", charges=1)[4].mz == pytest.approx(B5_PEPTM + 10.0, abs=1e-5)
    # The winning tag is a delta mass, which has no composition in any position.
    with pytest.raises(CompositionError):
        pt.comp(seq)


@pytest.mark.parametrize("placement", PLACEMENT_TAGS)
def test_a_placement_tag_in_front_of_a_named_modification(placement: str) -> None:
    """A named modification has a composition, so both paths give the same number."""
    seq = f"PEPTM[{placement}|Oxidation]IDE"
    assert pt.mass(seq) == pytest.approx(PEPTMIDE + OXYGEN, abs=1e-5)
    assert pt.mass(seq, calculate_with_composition=True) == pytest.approx(PEPTMIDE + OXYGEN, abs=1e-5)
    assert pt.comp(seq) == pt.comp("PEPTM[Oxidation]IDE")


@pytest.mark.parametrize("tags", ["[Position:M]", "[Limit:2]", "[CoMKP]", "[CoMUP]", "[CoMKP|CoMUP]", "[Limit:2|INFO:a]"])
def test_a_list_of_placement_tags_alone_is_unchanged(tags: str) -> None:
    """No tag carries a mass, so the first tag stands and the refusal it gives is kept.

    An information tag written first still reads as zero (``[INFO:a|Limit:2]`` weighs nothing),
    which is what peptacular has always done for a list that carries no mass at all.
    """
    with pytest.raises(PeptacularError, match="PlacementTag has no mass"):
        pt.mass(f"PEPTM{tags}IDE")
    assert pt.mass("PEPTM[INFO:a|Limit:2]IDE") == pytest.approx(PEPTMIDE, abs=1e-5)


def test_an_unknown_modification_behind_an_information_tag_is_reported() -> None:
    """The mass used to come from the information tag, so the unknown name went unnoticed."""
    with pytest.raises(pt.UnknownModificationError, match="NotAMod"):
        pt.mass("PEPTM[INFO:a|NotAMod]IDE")
    with pytest.raises(pt.UnknownModificationError, match="NotAMod"):
        pt.mass("PEPTM[NotAMod|INFO:a]IDE")


def test_a_secondary_ambiguity_site_weighs_nothing_extra() -> None:
    """A ``#g1`` reference marks a possible site; the mass belongs to the group's one modification.

    ProForma 2.1 section 7.6.2 writes the modification once with its label and the other sites as
    bare labels, so neither the bare label nor a label joined to an information tag may add mass.
    """
    plain = pt.mass("EM[Oxidation]EVTS")
    for seq in [
        "EM[Oxidation#g1]EVT[#g1]S",
        "EM[Oxidation#g1]EVT[INFO:x|#g1]S",
        "EM[Oxidation#g1]EVT[#g1|INFO:x]S",
        "EM[Oxidation#g1(0.8)]EVT[#g1(0.2)]S",
    ]:
        assert pt.mass(seq) == pytest.approx(plain, abs=1e-9), seq
        assert pt.parse(seq).serialize() == seq


# ---------------------------------------------------------------- validation of every tag

# ProForma 2.1 Final, section 7.5 and section 6.6, joined tag lists the specification prints as
# proper usage. Every one of these must stay valid.
SPEC_VALID_LISTS = [
    "ELVIS[U:Phospho|+79.966331]K",
    "ELVIS[U:Phospho|Obs:+79.978]K",
    "ELVIS[Phospho|O-phospho-L-serine]K",
    "ELVIS[UNIMOD:21|MOD:00046]K",
    "ELVIS[UNIMOD:21|Phospho]K",
    "ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966]K",
    "ELVIS[Obs:+79.966|Phospho|Sulfo]K",
    "ELVIS[Phospho|INFO:newly discovered]K",
    "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome [new-to-science] #cool]K",
    "ELN[Glycan:Hex1|INFO:WURCS=2.0/1,1,0/[a1211h-1b_1-4]/1/]IS",
    "ELV[+42#g1|INFO:newly discovered]I[#g1|INFO:awesome]S",
    "[+42#g1]?ELV[INFO:newly discovered|#g1]I[INFO:awesome|#g1]S",
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
    "[INFO:FLANK:ELV]-PRECIS-[INFO:FLANK:ELY]",
]

# One list per position kind, each with an unknown modification name written second.
UNKNOWN_SECOND = [
    "PEPTM[INFO:a|NoSuchMod]IDE",
    "[INFO:a|NoSuchMod]-PEPTIDE",
    "PEPTIDE-[INFO:a|NoSuchMod]",
    "{INFO:a|NoSuchMod}PEPTIDE",
    "[INFO:a|NoSuchMod]?PEPTIDE",
    "PEP(TM)[INFO:a|NoSuchMod]IDE",
    "<[INFO:a|NoSuchMod]@M>PEPTMIDE",
]

# The same kinds with the unknown name written third, and with a malformed accession instead.
UNKNOWN_THIRD = [
    "PEPTM[INFO:a|INFO:b|NoSuchMod]IDE",
    "[INFO:a|INFO:b|NoSuchMod]-PEPTIDE",
    "PEPTIDE-[INFO:a|INFO:b|NoSuchMod]",
    "{INFO:a|INFO:b|NoSuchMod}PEPTIDE",
    "[INFO:a|INFO:b|NoSuchMod]?PEPTIDE",
    "PEP(TM)[INFO:a|INFO:b|NoSuchMod]IDE",
    "<[INFO:a|INFO:b|NoSuchMod]@M>PEPTMIDE",
    "PEPTM[INFO:a|UNIMOD:99999999]IDE",
    "PEPTM[INFO:a|INFO:b|MOD:99999999]IDE",
]


@pytest.mark.parametrize("seq", UNKNOWN_SECOND + UNKNOWN_THIRD)
def test_validation_reports_a_bad_tag_in_any_position(seq: str) -> None:
    """Validation checks every tag of a list, not only the first. Each row is True on main."""
    assert pt.validate(seq) is False
    with pytest.raises(PeptacularError):
        pt.parse(seq, validate=True)
    # The mass path agrees, because the mass comes from that same tag.
    with pytest.raises(PeptacularError):
        pt.mass(seq)


@pytest.mark.parametrize("seq", SPEC_VALID_LISTS)
def test_the_specifications_own_joined_lists_stay_valid(seq: str) -> None:
    assert pt.validate(seq) is True
    assert pt.parse(seq, validate=True).serialize() == seq


def test_the_specifications_own_invalid_info_tag_is_now_refused() -> None:
    """Section 6.6, on three comments it prints as wrong:

        "The first two comments would be invalid because of an unpaired bracket and the last is
        invalid because of the use of the pipe character:
        ELVIS[Phospho|INFO:newly]discovered]K
        ELVIS[Phospho|INFO:newly[discovered]K
        ELVIS[Phospho|INFO:newly|discovered]K"

    The third splits into ``Phospho``, ``INFO:newly`` and ``discovered``, and that last tag is
    not a modification anyone knows. It read as valid while only the first tag was checked.
    """
    assert pt.validate("ELVIS[Phospho|INFO:newly|discovered]K") is False
    # The other two never parsed, because of the unpaired bracket.
    assert pt.validate("ELVIS[Phospho|INFO:newly]discovered]K") is False
    assert pt.validate("ELVIS[Phospho|INFO:newly[discovered]K") is False


def test_all_tags_false_still_checks_the_first_tag_only() -> None:
    """The parameter stays accepted, for a caller who wants the older, narrower check."""
    mod = ModificationTags.from_string("INFO:a|NoSuchMod")
    assert mod.validate() is not None
    assert mod.validate(all_tags=True) is not None
    assert mod.validate(all_tags=False) is None
    assert ModificationTags.from_string("NoSuchMod|INFO:a").validate(all_tags=False) is not None
