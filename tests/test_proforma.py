import pytest

import peptacular as pt


@pytest.fixture
def basic_mod():
    return pt.Mod("HelloWorld", 1)


def test_parse():
    # Test the case where there are no modifications in the peptide _sequence
    proforma = "PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE")
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_parse_empty():
    # Test the case where there are no modifications in the peptide _sequence
    proforma = ""
    expected_output = pt.ProFormaAnnotation(sequence="")
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_single_mod(basic_mod: pt.Mod):
    # Test the case where there is a single modification in the peptide _sequence
    proforma = "PEPTIDE[HelloWorld]"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", internal_mods={6: {"HelloWorld": 1}})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_single_double_mod():
    # Test the case where there is a single modification in the peptide _sequence
    proforma = "PEPTIDE[Hello][World]"
    expected_output = pt.ProFormaAnnotation(
        sequence="PEPTIDE",
        internal_mods={6: {"Hello": 1, "World": 1}},
    )
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_single_double_multi_mod():
    # Test the case where there is a single modification in the peptide _sequence
    proforma = "PEPTIDE[Hello]^2[World]^3"
    # assert raises ValueError
    try:
        pt.ProFormaAnnotation.parse(proforma)
        raise AssertionError("Expected ValueError not raised")
    except ValueError:
        pass


def test_nterm_mod(basic_mod: pt.Mod):
    # Test the case where there is a single N-Term mod
    proforma = "[HelloWorld]-PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"HelloWorld": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_nterm_double_mod():
    # Test the case where there is a single N-Term mod
    proforma = "[Hello][World]-PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"Hello": 1, "World": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_nterm_double_multi_mod():
    # Test the case where there is a single N-Term mod
    proforma = "[Hello]^2[World]^3-PEPTIDE"
    try:
        pt.ProFormaAnnotation.parse(proforma)
        raise AssertionError("Expected ValueError not raised")
    except ValueError:
        pass


def test_cterm_mod(basic_mod: pt.Mod):
    # Test the case where there is a single C-Term mod
    proforma = "PEPTIDE-[HelloWorld]"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", cterm_mods={"HelloWorld": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_cterm_double_mod():
    # Test the case where there is a single C-Term mod
    proforma = "PEPTIDE-[Hello][World]"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", cterm_mods={"Hello": 1, "World": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_cterm_double_multi_mod():
    # Test the case where there is a single C-Term mod
    proforma = "PEPTIDE-[Hello][Hello][World][World][World]"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", cterm_mods={"Hello": 2, "World": 3})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_multiple_mods(basic_mod: pt.Mod):
    # Test the case with multiple modifications at different positions
    proforma = "P[HelloWorld]EPT[HelloWorld]IDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", internal_mods={0: {"HelloWorld": 1}, 3: {"HelloWorld": 1}})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_interval_mods(basic_mod: pt.Mod):
    # Test the case where there is an interval with modifications
    proforma = "PEP(TI)[HelloWorld]DE"
    expected_interval = pt.Interval(start=3, end=5, ambiguous=False, mods={"HelloWorld": 1})
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", intervals=[expected_interval])
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_interval_mods_double():
    proforma = "PEP(TI)[Hello][World]DE"
    expected_interval = pt.Interval(
        start=3,
        end=5,
        ambiguous=False,
        mods={"Hello": 1, "World": 1},
    )
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", intervals=[expected_interval])
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_ambiguous_interval():
    # Test an ambiguous interval
    proforma = "PEP(?TI)DE"
    expected_interval = pt.Interval(start=3, end=5, ambiguous=True, mods=None)
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", intervals=[expected_interval])
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_ambiguous_interval_with_mod(basic_mod: pt.Mod):
    # Test an ambiguous interval with a modification
    proforma = "PEP(?TI)[HelloWorld]DE"
    expected_interval = pt.Interval(start=3, end=5, ambiguous=True, mods={"HelloWorld": 1})
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", intervals=[expected_interval])
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_charge_state():
    # Test the case with a charge state at the end of the _sequence
    proforma = "PEPTIDE/2"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", charge=2)
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_negative_charge_state():
    # Test the case with a negative charge state at the end of the _sequence
    proforma = "PEPTIDE/-2"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", charge=-2)
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_charge_state_with_mod(basic_mod: pt.Mod):
    # Test the case with a charge state at the end of the _sequence
    proforma = "PEPTIDE[HelloWorld]/2"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", internal_mods={6: {"HelloWorld": 1}}, charge=2)
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_charge_state_with_adduct():
    # Test the case with a charge state at the end of the _sequence
    proforma = "PEPTIDE/[H:z+1]"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", charge="H:z+1")
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_labile_mod(basic_mod: pt.Mod):
    # Test the case where there is a labile modification
    proforma = "{HelloWorld}PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", labile_mods={"HelloWorld": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_labile_mod_double():
    # Test the case where there is a labile modification
    proforma = "{Hello}{World}PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", labile_mods={"Hello": 1, "World": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_labile_mod_double_multi():
    # Test the case where there is a labile modification
    proforma = "{Hello}{World}{World}PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", labile_mods={"Hello": 1, "World": 2})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_unknown_mod(basic_mod: pt.Mod):
    # Test the case with an unknown modification
    proforma = "[HelloWorld]?PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", unknown_mods={"HelloWorld": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_unknown_mod_double():
    # Test the case with an unknown modification
    proforma = "[Hello][World]?PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", unknown_mods={"Hello": 1, "World": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_unknown_mod_double_multi():
    # Test the case with an unknown modification
    proforma = "[Hello]^2[World]^3?PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", unknown_mods={"Hello": 2, "World": 3})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_isotope_mod():
    # Test the case with an isotope modification
    proforma = "<13C>PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", isotope_mods={"13C": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_isotope_mod_double():
    # Test the case with an isotope modification
    proforma = "<13C><15N>PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", isotope_mods={"13C": 1, "15N": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_isotope_mod_double_multi():
    # Test the case with an isotope modification
    proforma = "<13C>^2<15N>^3PEPTIDE"
    # expect ValueError
    with pytest.raises(ValueError):
        pt.ProFormaAnnotation.parse(proforma)


def test_static_mod_with_target():
    # Test static modification with a specified target
    proforma = "<[Phospho]@S>PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[Phospho]@S": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_static_mod_with_target_double():
    # Test static modification with a specified target
    proforma = "<[Phospho]@S><[Acetyl]@K>PEPTIDE"
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[Phospho]@S": 1, "[Acetyl]@K": 1})
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_static_mod_with_target_double_multi():
    # Test static modification with a specified target
    proforma = "<Phospho@S>^2<Acetyl@K>^3PEPTIDE"
    try:
        pt.ProFormaAnnotation.parse(proforma)
        raise AssertionError("Expected ValueError not raised")
    except ValueError:
        pass


def test_multiple_types_of_mods():
    # Test a _sequence with multiple types of modifications
    proforma = "[Acetyl]-PEP[T[Phospho]]TIDE-[Amide]"
    expected_output = pt.ProFormaAnnotation(
        sequence="PEPTIDE",
        nterm_mods={"Acetyl": 1},
        cterm_mods={"Amide": 1},
        internal_mods={2: {"T[Phospho]": 1}},
    )
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_mod_with_multiplier():
    # Test modification with a multiplier
    proforma = "PEPTIDE[HelloWorld]^3"
    try:
        pt.ProFormaAnnotation.parse(proforma)
        raise AssertionError("Expected ValueError not raised")
    except ValueError:
        pass


def test_sequence_with_multiple_intervals():
    # Test a _sequence with multiple intervals

    proforma = "PEP(TI)[Mod1](DE)[Mod2]"
    expected_intervals = [
        pt.Interval(start=3, end=5, ambiguous=False, mods={"Mod1": 1}),
        pt.Interval(start=5, end=7, ambiguous=False, mods={"Mod2": 1}),
    ]
    expected_output = pt.ProFormaAnnotation(sequence="PEPTIDE", intervals=expected_intervals)
    assert pt.ProFormaAnnotation.parse(proforma) == expected_output
    assert pt.serialize(expected_output) == proforma


def test_slice_annotation():
    # Test slicing a ProForma annotation
    proforma = "P[Phospho]EPTIDE"
    annotation = pt.ProFormaAnnotation.parse(proforma)
    sliced = annotation.slice(0, 4, inplace=False)
    assert sliced.sequence == "PEPT"
    assert sliced.serialize() == "P[Phospho]EPT"


def test_reverse_annotation():
    # Test reversing a ProForma annotation
    annotation = pt.ProFormaAnnotation.parse("P[Phospho]EPTIDE")
    reversed_annotation = annotation.reverse()
    assert reversed_annotation.serialize() == "EDITPEP[Phospho]"


def test_shuffle_annotation():
    # Test shuffling a ProForma annotation with a fixed seed
    annotation = pt.ProFormaAnnotation.parse("P[Phospho]EPTIDE")
    shuffled = annotation.shuffle(seed=42)
    assert len(shuffled.serialize()) == len(annotation.serialize())
    assert shuffled.serialize() == "ETIPEP[Phospho]D"


def test_shift_annotation():
    # Test shifting a ProForma annotation
    proforma = "P[Phospho]EPTIDE"
    annotation = pt.ProFormaAnnotation.parse(proforma)
    shifted = annotation.shift(2)
    assert shifted.serialize() == "PTIDEP[Phospho]E"


def test_is_subsequence():
    # Test subsequence detection
    parent = pt.ProFormaAnnotation.parse("PEPTIDEPEPTIDE")
    child = pt.ProFormaAnnotation.parse("TIDE")
    assert child.is_subsequence(parent)


def test_find_indices():
    # Test finding indices of a subsequence
    parent = pt.ProFormaAnnotation.parse("PEPTIDEPEPTIDE")
    child = pt.ProFormaAnnotation.parse("PEP")
    indices = child.find_indices(parent)
    assert indices == [0, 7]


def test_count_residues():
    # Test counting residues
    annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
    counts = annotation.count_residues()
    assert counts["P"] == 2
    assert counts["E"] == 2
    assert counts["T"] == 1
    assert counts["I"] == 1
    assert counts["D"] == 1


def test_count_residues_with_modifications():
    # Test counting residues with modifications
    proforma = "P[Phospho]EPTIDE"
    annotation = pt.ProFormaAnnotation.parse(proforma)
    counts = annotation.count_residues()
    assert counts["P"] == 1
    assert counts["P[Phospho]"] == 1
    assert counts["E"] == 2
    assert counts["T"] == 1
    assert counts["I"] == 1
    assert counts["D"] == 1
