"""append_mods takes one or several mods per key; extend_mods treats a bare string as one mod."""

import pytest

import peptacular as pt


@pytest.mark.parametrize(
    "mods, expected",
    [
        ({1: ["Oxidation"]}, "PE[Oxidation]MTIDE"),
        ({1: ("Oxidation", "Phospho", "Acetyl")}, "PE[Oxidation][Phospho][Acetyl]MTIDE"),
        ({"nterm": ["Acetyl", "Formyl"]}, "[Acetyl][Formyl]-PEMTIDE"),
        ({"cterm": ["Amidated"]}, "PEMTIDE-[Amidated]"),
        ({"labile": ["Hex"]}, "{Hex}PEMTIDE"),
        ({"unknown": ["Phospho"]}, "[Phospho]?PEMTIDE"),
        ({1: "Oxidation"}, "PE[Oxidation]MTIDE"),
        ({1: ("Oxidation", 2)}, "PE[Oxidation][Oxidation]MTIDE"),
    ],
)
def test_append_mods_accepts_single_and_many(mods, expected):
    assert pt.append_mods("PEMTIDE", mods) == expected
    assert pt.parse("PEMTIDE").append_mods(mods).serialize() == expected


@pytest.mark.parametrize(
    "mods, expected",
    [
        ({1: "Oxidation"}, "PE[Oxidation]MTIDE"),
        ({"nterm": "Acetyl"}, "[Acetyl]-PEMTIDE"),
        ({"cterm": "Amidated"}, "PEMTIDE-[Amidated]"),
        ({"labile": "Hex"}, "{Hex}PEMTIDE"),
        ({"unknown": "Phospho"}, "[Phospho]?PEMTIDE"),
        ({1: ["Oxidation", "Phospho"]}, "PE[Oxidation][Phospho]MTIDE"),
        ({1: 15.995}, "PE[+15.995]MTIDE"),
    ],
)
def test_extend_mods_treats_str_as_one_mod(mods, expected):
    assert pt.extend_mods("PEMTIDE", mods) == expected
    assert pt.parse("PEMTIDE").extend_mods(mods).serialize() == expected
