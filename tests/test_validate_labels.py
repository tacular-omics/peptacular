"""The ``#label`` rules of ``validate()``, against ProForma 2.1 Final (9 June 2026).

Quotes below were re-read from the specification's text layer on 2026-11-20:

- 7.6.2 (ii): "A single preferred location for the modification MUST be specified, so that
  the sequence can be easily rendered in visualization tools."
- 7.6.2: "Multiple modifications of unknown position can be defined in a single
  definition. These MUST use different labels. Labels that only differ in case MUST NOT be
  used."
- 7.9: "Such a label MUST not overlap with labels used for ambiguous modifications
  (section 7.6.2) or cross-linkers (section 9.2)."
- 9.2: "Crosslinked sites MUST be represented by immediately following the modification
  notation using the prefix #XL, followed by an arbitrary label consisting of alphanumeric
  characters."
- 12.11: "A labile modification (see section 6.4) and modifications of unknown position
  (see section 7.6) have to be defined as part of a peptidoform."

Appendix I, the formal grammar, on where a label may be written:

- ``modOrLabel = mod | ("[", modLabel, "]");``
- ``mod = "[", modInternal, "]";``
- ``modInternal = modSingle, [modLabel], {"|", modSingle, [modLabel]};``
- ``modLabel = "#", (modLabelXL | modLabelBranch | modLabelAmbiguous);`` with
  ``modLabelAmbiguous = ALPHANUMERIC, {ALPHANUMERIC}, ["(", SIGNEDNUMBER, ")"];`` and
  ``ALPHANUMERIC = LETTER|DIGIT;``

Every string in ``SPEC_VALID`` is copied from the specification, with the section that
prints it, and the specification presents each one as valid.

``SPEC_REFUSED_FOR_AN_UNKNOWN_MODIFICATION`` holds the label-bearing spec strings that
peptacular refuses for the modification rather than for its label, on main as well as
here, so they are not expected to validate. They are not the only spec strings refused
that way, and the cause has nothing to do with labels: the same XL-MOD name without any
label, ``EVTK[X:formaldehyde]LEK[XLMOD:02041]SEF``, is refused too.

Strings joined by ``//`` are left out: peptacular cannot parse a cross-linked pair into
one annotation, so every such string is already invalid, before this change.
"""

import pytest

import peptacular as pt
from peptacular.diagnostics import PeptacularError, UnknownModificationError

SPEC_VALID = [
    ("5.3.2", "PEP[Oxidation#1]M[#1]AT"),
    ("5.3.2", "PEP[Oxidation#1(0.95)]M[#1(0.05)]AT"),
    ("5.3.2", "P(EP)[Oxidation#1(0.95)]M[#1(0.05)]AT"),
    ("5.3.2", "{Phospho#g1}EM[U:Oxidation]EVS[#g1]"),
    ("5.3.2", "{Sulfo#s}HT[#s]T[#s]S[#s]FRGC[#s]"),
    ("7.6.2", "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK"),
    ("7.6.2", "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1][Formyl]PEK"),
    ("7.6.2", "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1][Formyl#g2]PEK[#g2]"),
    ("7.6.4", "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK"),
    ("7.6.4", "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.09)]ES[#s1(0.90)]PEK"),
    ("7.6.5", "PRT(ESFRMS)[+19.0523#g1(0.01)]ISK[#g1(0.99)]"),
    ("7.6.5", "PR[#g1(0.91)]T(EC[Carbamidomethyl]FRMS)[+19.05233#g1(0.09)]ISK"),
    ("7.7", "EN(?G[+100.00#g1]G[#g1]E[#g1])KL"),
    ("7.7", "[+100.00#g1]?EN(?G[#g1]G[#g1]E[#g1])KL"),
    ("7.9", "{HexNAc#g1}PEPN[#g1]ITE"),
    ("7.9", "{Hex#g1}{Sulfo#s}[#g1]-PN[#g1]KS[#s]C[#s]M"),
    ("7.9", "{Glycan:Hex5HexNAc5#g1}PEPSTAT[#g1]IS[#g1]T[#g1]ICS[#g1]S[#g1]T[#g1]RIKES[#g1]IT[#g1]ES[#g1]"),
    ("9.2.1", "EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]"),
    ("9.2.1", "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR"),
    ("9.2.1", "EMEVTK[XLMOD:02001#XL1]SESPEK"),
    ("9.2.3", "EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD"),
    ("9.2.3", "EVTSEKC[L-cystine (cross-link)#XL1]LEMSC[#XL1]EFD"),
    ("9.2.3", "EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD"),
    ("9.2.3", "EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD"),
    ("9.2.3", "EVTSEKC[Dehydro#XL1]LEMSC[#XL1][Dehydro]EFD"),
    ("11.2", "PETIEM[Dioxidation#1][Oxidation#2]REM[#1][#2]REM[#2]RM[#1]PEPTIDE"),
    ("11.2", "[Oxidation#1][Oxidation#2]?PTIM[#1][#2]ERM[#1][#2]ERM[#1][#2]ERM[#1][#2]DE"),
]

SPEC_REFUSED_FOR_AN_UNKNOWN_MODIFICATION = [
    # 9.2.3 (iii): the term "is not present in the Unimod OBO file, see section 6.2.1.2,
    # and so SHOULD not be used".
    ("9.2.3", "EVTSEKC[UNIMOD:2020#XL1]LEMSC[#XL1]EFD"),
    ("9.2.3", "EVTSEKC[Xlink:Disulfide#XL1]LEMSC[#XL1]EFD"),
    # The XL-MOD name of the 5.3.4 feature table, refused with or without its label.
    ("5.3.4", "EVTK[X:formaldehyde#XL1]LEK[#XL1]SEF"),
    # Section 12.1, among the use cases "NOT currently supported in the current version".
    ("12.1", "LEIK[N6-(L-asparagyl)-L-lysine#XL1]KIPHDN[#XL1]"),
]

# One label never defined, in each position kind. Every one of these returns True on main.
NEVER_DEFINED = [
    ("internal, two references", "PEPT[#g1]IDES[#g1]"),
    ("internal, one reference", "PEPTIDES[#g1]"),
    ("cross-linker label", "PEPK[#XL1]IDE"),
    ("branch label", "D[#BRANCH]ATER"),
    ("N-terminal", "[#g1]-PEPTIDE"),
    ("C-terminal", "PEPTIDE-[#g1]"),
    ("unknown position", "[#g1]?PEPTIDE"),
    ("interval", "PEPT(IDE)[#g1]S"),
    ("second label of a pair", "PEPT[Phospho#g1]S[#g1]IDEK[#g2]"),
    ("score on an undefined label", "PEPT[#g1(0.90)]IDES[#g1(0.10)]"),
    ("reference in the other case", "PEPT[Phospho#g1]S[#G1]IDE"),
]

# Both labels of each pair are defined, so only the case rule can refuse these.
CASE_ONLY = [
    ("two ambiguous groups", "PEPT[Phospho#g1]S[#g1]IDES[Formyl#G1]K[#G1]"),
    ("two unknown-position groups", "[Phospho#g1][Formyl#G1]?PEPTIDE"),
    ("labile against ambiguous", "{Phospho#g1}PEPT[Formyl#G1]IDEK[#G1]S[#g1]"),
]

LABILE_OVERLAP = [
    ("labile and an ambiguous group", "{Glycan:Hex#g1}PEPT[Phospho#g1]IDES[#g1]"),
    ("labile and a cross-linker", "{Glycan:Hex#XL1}PEPTK[XLMOD:02001#XL1]IDEK[#XL1]"),
    ("labile and a named modification", "{Phospho#s}PEPT[Formyl#s]IDE"),
]

LABILE_REFERENCED = [
    ("7.9", "{HexNAc#g1}PEPN[#g1]ITE"),
    ("7.9", "{Glycan:Hex#g1}PEPTIDET[#g1]"),
    ("7.9", "{Hex#g1}{Sulfo#s}[#g1]-PN[#g1]KS[#s]C[#s]M"),
]


@pytest.mark.parametrize(("section", "sequence"), SPEC_VALID)
def test_every_label_example_the_specification_presents_as_valid_still_validates(section: str, sequence: str) -> None:
    assert pt.validate(sequence) is True, f"section {section}: {sequence}"


@pytest.mark.parametrize(("section", "sequence"), SPEC_REFUSED_FOR_AN_UNKNOWN_MODIFICATION)
def test_the_spec_strings_refused_for_the_modification_and_not_for_the_label(section: str, sequence: str) -> None:
    assert pt.validate(sequence) is False, f"section {section}: {sequence}"
    with pytest.raises(UnknownModificationError):
        pt.mass(sequence)


def test_that_cause_is_the_modification_name_and_not_the_label() -> None:
    # The same XL-MOD name with no label at all, from the 5.3.4 feature table.
    assert pt.validate("EVTK[X:formaldehyde]LEK[XLMOD:02041]SEF") is False
    with pytest.raises(UnknownModificationError):
        pt.mass("EVTK[X:formaldehyde]LEK[XLMOD:02041]SEF")


@pytest.mark.parametrize(("kind", "sequence"), NEVER_DEFINED)
def test_a_label_that_is_never_defined_is_refused(kind: str, sequence: str) -> None:
    assert pt.validate(sequence) is False, kind
    with pytest.raises(PeptacularError, match="never defined"):
        pt.parse(sequence).validate_annotation()


def test_a_never_defined_label_loses_the_modification_in_silence() -> None:
    # 886.39199243098 is pyteomics 5.0.1 mass.calculate_mass(sequence="PEPTIDES"),
    # measured 2026-11-20: the labelled group in the sequence contributes nothing.
    assert pt.mass("PEPT[#g1]IDES[#g1]") == pytest.approx(886.39199243098, abs=1e-6)
    assert pt.mass("PEPTIDES") == pytest.approx(886.39199243098, abs=1e-6)


def test_a_defined_label_without_any_reference_stays_valid() -> None:
    # Section 9.2.1 calls this a "dead end" crosslink, and 7.6.2 says a group with exactly
    # one position "should be treated as a normal modification".
    assert pt.validate("EMEVTK[XLMOD:02001#XL1]SESPEK") is True
    assert pt.validate("PEPT[Phospho#g1]IDE") is True
    assert pt.validate("[Phospho#g1]?PEPTIDE") is True
    assert pt.validate("{Glycan:Hex#g1}PEPTIDE") is True


# A label written as its own element of a bracket that also holds a modification tag. The
# grammar of Appendix I does not admit this spelling while section 6.6 prints it, so
# peptacular keeps the reading it has always given it and these stay valid, as on main.
LABEL_AS_A_PIPE_ELEMENT = [
    "PEPT[Phospho|#g1]IDES[#g1]",
    "PEPT[#g1|Phospho]IDES[#g1]",
    "PEPT[+79.9663|#g1]IDES[#g1]",
    "ELV[+42#g1|INFO:newly discovered]I[#g1|INFO:awesome]S",
    "[+42#g1]?ELV[INFO:newly discovered|#g1]I[INFO:awesome|#g1]S",
]

# One label, two modifications. Every one of these returns True on main.
DEFINED_TWICE = [
    ("two labile modifications", "{Hex#g1}{HexNAc#g1}PEPT[#g1]IDE"),
    ("the same labile modification twice", "{Hex#g1}{Hex#g1}PEPT[#g1]IDE"),
    ("unknown position and a residue", "[Phospho#g1]?PEPT[Formyl#g1]IDE"),
    ("labile and unknown position", "{Hex#g1}[Phospho#g1]?PEPTIDE"),
    ("a multiplied unknown-position modification", "[Oxidation#1]^2?PEPTIDE"),
    ("the same residue modification twice", "PEPT[Phospho#g1][Phospho#g1]IDE"),
    ("two labels as pipe elements", "PEPT[Phospho|#g1]IDES[Formyl|#g1]K"),
]


@pytest.mark.parametrize("sequence", LABEL_AS_A_PIPE_ELEMENT)
def test_a_label_written_as_its_own_pipe_element_defines_it(sequence: str) -> None:
    assert pt.validate(sequence) is True, sequence
    assert pt.parse(sequence).serialize() == sequence


def test_a_label_in_a_bracket_of_its_own_is_a_reference() -> None:
    # Appendix I reads "[#g1]" as a modOrLabel with no modSingle of its own, so the
    # modification in the bracket before it does not carry the label.
    assert pt.validate("PEPT[Oxidation][#g1]IDE") is False
    assert pt.validate("PEPT[Oxidation][#g1]IDES[Phospho#g1]") is True
    # A bracket holding labels, scores and information tags only stays a reference.
    assert pt.validate("PEPT[Phospho#g1]IDES[INFO:a|#g1]") is True
    assert pt.validate("PEPT[Phospho#g1]IDES[#g1(0.10)]") is True


@pytest.mark.parametrize(("kind", "sequence"), DEFINED_TWICE)
def test_one_label_on_two_modifications_is_refused(kind: str, sequence: str) -> None:
    # Section 7.6.2: "These MUST use different labels", and (ii) a single preferred
    # location MUST be specified.
    assert pt.validate(sequence) is False, kind


def test_two_labels_differing_only_in_case_are_compared_as_ascii() -> None:
    # Labels are ALPHANUMERIC in the Appendix I grammar, so the comparison is lower(),
    # not casefold(), and the eszett that casefold() turns into "ss" is not a clash.
    assert pt.validate("PEPT[Phospho#ß]IDES[#ß]S[Formyl#ss]K[#ss]") is True
    assert pt.validate("PEPT[Phospho#g1]S[#g1]IDES[Formyl#G1]K[#G1]") is False


@pytest.mark.parametrize(("kind", "sequence"), CASE_ONLY)
def test_labels_that_differ_only_in_case_are_refused(kind: str, sequence: str) -> None:
    assert pt.validate(sequence) is False, kind
    with pytest.raises(PeptacularError, match="differing only in case"):
        pt.parse(sequence).validate_annotation()


def test_the_same_label_in_one_case_is_not_a_case_clash() -> None:
    assert pt.validate("PEPT[Phospho#g1]S[#g1]IDES[Formyl#g2]K[#g2]") is True
    assert pt.validate("PEPT[Phospho#G1]S[#G1]IDE") is True


@pytest.mark.parametrize(("kind", "sequence"), LABILE_OVERLAP)
def test_a_labile_label_reused_for_another_modification_is_refused(kind: str, sequence: str) -> None:
    assert pt.validate(sequence) is False, kind
    with pytest.raises(PeptacularError, match="MUST not overlap"):
        pt.parse(sequence).validate_annotation()


@pytest.mark.parametrize(("section", "sequence"), LABILE_REFERENCED)
def test_a_labile_label_referenced_in_the_sequence_stays_valid(section: str, sequence: str) -> None:
    assert pt.validate(sequence) is True, f"section {section}: {sequence}"


def test_the_specifications_own_invalid_group_is_still_refused() -> None:
    # Section 7.6.2: "The following example is not valid because a single preferred
    # location MUST be chosen for a modification".
    assert pt.validate("EM[Oxidation]EVT[#g1]S[Phospho#g1]ES[Phospho#g1]PEK") is False


def test_a_valid_string_is_not_a_promise_of_a_mass() -> None:
    # The docstrings of pt.validate and validate_annotation say this; BZJX are valid
    # residues (section 7.1) and an empty sequence is valid ProForma.
    for sequence in ["PEPXTIDE", "PEPBTIDE", "PEPZTIDE", ""]:
        assert pt.validate(sequence) is True, sequence
    for sequence in ["PEPXTIDE", "PEPBTIDE", "PEPZTIDE", ""]:
        with pytest.raises(PeptacularError):
            pt.mass(sequence)
    assert pt.mass("PEPX[+0]TIDE") == pytest.approx(799.35996402671, abs=1e-6)


def test_parse_with_validation_on_runs_the_field_checks_only() -> None:
    # parse(validate=True) validates each field as it is set, so it cannot run a check that
    # reads the whole annotation: a label may be defined in an interval that is set after the
    # positions referencing it, and slicing a labelled group legitimately leaves references
    # behind. The three new rules are therefore run by validate() and validate_annotation()
    # only, while validate_ambiguous_labels keeps firing there as it did before, now over
    # labile and unknown-position modifications as well.
    pt.parse("PEPT[#g1]IDES[#g1]", validate=True)
    pt.parse("PEPT[Phospho#g1]S[#g1]IDES[Formyl#G1]K[#G1]", validate=True)
    assert pt.validate("PEPT[#g1]IDES[#g1]") is False
    assert pt.validate("PEPT[Phospho#g1]S[#g1]IDES[Formyl#G1]K[#G1]") is False
    with pytest.raises(PeptacularError, match="more than one concrete modification"):
        pt.parse("{Glycan:Hex#g1}PEPT[Phospho#g1]IDES[#g1]", validate=True)
