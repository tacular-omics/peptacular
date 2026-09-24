"""5.0: every user-input error on a public path is a PeptacularError."""

import pickle

import pytest

import peptacular as pt

USER_INPUT_ERRORS = {
    "parallel method name": lambda: pt.mass(["PEPTIDE"], method="threads"),  # ty: ignore[invalid-argument-type]
    "digest parallel method name": lambda: pt.digest(["PEPTIDEK"], "trypsin", method="threads"),  # ty: ignore[invalid-argument-type]
    "batch parallel method name": lambda: list(pt.iter_batch("mass", ["PEPTIDE"], method="threads")),  # ty: ignore[invalid-argument-type]
    "batch start method": lambda: list(pt.iter_batch("mass", ["PEPTIDE"], method="process", start_method="bogus")),
    "set_start_method": lambda: pt.set_start_method("bogus"),  # ty: ignore[invalid-argument-type]
    "secondary structure scale": lambda: pt.secondary_structure("PEPTIDE", scale="nope"),
    "permutations size": lambda: pt.permutations("PEPTIDE", -1),
    "combinations r": lambda: list(pt.parse("PEPTIDE").combinations(-1)),
    "product repeat": lambda: list(pt.parse("PEP").product(-2)),
    "ms2pip location": lambda: pt.from_ms2_pip(("PEPTIDE", "x|Oxidation")),
    "proforma json syntax": lambda: pt.from_proforma_json("{bad"),
    "join nothing": lambda: pt.join([]),
    "chem_comp element": lambda: pt.chem_comp({"Qq": 1}),
    "chem_mass element": lambda: pt.chem_mass({"Qq": 1}),
    "chem_formula element": lambda: pt.chem_formula({"Qq": 1}),
    "brain element": lambda: pt.brain_isotopic_distribution({"Xx": 3}),
    "averagine ion type": lambda: pt.estimate_averagine_comp(500, ion_type="q"),  # ty: ignore[invalid-argument-type]
    "negative missed cleavages": lambda: pt.digest("PEPTIDEKAAR", "trypsin", missed_cleavages=-1),
    "negative missed cleavages (span method)": lambda: list(pt.parse("PEPTIDEKAAR").digest_spans("trypsin", missed_cleavages=-1)),
    "annotation integer index": lambda: pt.parse("PEPTIDE")[2],  # ty: ignore[invalid-argument-type]
}


@pytest.mark.parametrize("call", USER_INPUT_ERRORS.values(), ids=USER_INPUT_ERRORS.keys())
def test_user_input_errors_are_peptacular_errors(call):
    with pytest.raises(pt.PeptacularError):
        call()


@pytest.mark.parametrize("call", [lambda: pt.chem_comp({"Qq": 1}), lambda: pt.brain_isotopic_distribution({"Xx": 3})])
def test_unknown_element_is_a_readable_key_error(call):
    with pytest.raises(pt.UnknownElementError) as info:
        call()
    assert isinstance(info.value, KeyError)
    assert isinstance(info.value, pt.PeptacularKeyError)
    assert not str(info.value).startswith('"')


def test_unknown_enzyme_is_a_peptacular_key_error():
    assert issubclass(pt.UnknownEnzymeError, pt.PeptacularKeyError)


def test_zero_missed_cleavages_is_still_allowed():
    assert pt.digest("PEPTIDEKAAR", "trypsin", missed_cleavages=0) == pt.digest("PEPTIDEKAAR", "trypsin")


def test_json_error_keeps_the_decoder_error_as_cause():
    import json

    with pytest.raises(pt.PeptacularError) as info:
        pt.from_proforma_json("{bad")
    assert isinstance(info.value.__cause__, json.JSONDecodeError)


class TestModsSnapshot:
    def test_mods_view_does_not_follow_later_edits(self):
        annot = pt.parse("[Acetyl]-PEPTIDE")
        mods = annot.nterm_mods
        before_hash, before_mass, before_str = hash(mods), mods.get_mass(), str(mods)
        annot.append_mods({"nterm": "Formyl"})
        assert str(mods) == before_str
        assert mods.get_mass() == before_mass
        assert hash(mods) == before_hash
        assert str(annot.nterm_mods) == "Mods@Nterm([Acetyl, Formyl])"

    @pytest.mark.parametrize("attr", ["internal_mods", "labile_mods", "unknown_mods", "static_mods", "isotope_mods", "cterm_mods"])
    def test_every_mod_view_is_a_snapshot(self, attr):
        annot = pt.parse("<13C><[Oxidation]@M>{Glycan:Hex}[Phospho]?PEM[Deamidated]TIDE-[Amidated]")
        view = getattr(annot, attr)
        before = str(view)
        annot.clear_mods(inplace=True)
        assert str(view) == before

    def test_interval_mods_snapshot(self):
        annot = pt.parse("P(EP)[Phospho]TIDE")
        interval = annot.intervals[0]
        view = interval.mods
        interval.set_mods({"Oxidation": 1})
        assert str(view) == "Mods@Interval([Phospho])"

    def test_snapshot_is_read_only_and_hashable(self):
        mods = pt.parse("[Acetyl]-PEPTIDE").nterm_mods
        with pytest.raises(TypeError):
            mods._mods["Formyl"] = 1  # ty: ignore[invalid-assignment]
        assert hash(mods) == hash(pt.parse("[Acetyl]-PEPTIDE").nterm_mods)
        assert mods == pt.parse("[Acetyl]-PEPTIDE").nterm_mods

    def test_snapshot_round_trips(self):
        mods = pt.parse("[Acetyl]-PEPTIDE").nterm_mods
        assert pt.parse("PEPTIDE").set_nterm_mods(mods, inplace=False).serialize() == "[Acetyl]-PEPTIDE"
        assert mods.copy() == mods
        assert pickle.loads(pickle.dumps(pt.parse("[Acetyl]-PEPTIDE"))).nterm_mods == mods


def test_integer_index_hint_and_bad_key_type():
    annot = pt.parse("PEPTIDE")
    with pytest.raises(pt.UnsupportedOperationError, match=r"annot\[2:3\]"):
        annot[2]  # ty: ignore[invalid-argument-type]
    with pytest.raises(TypeError):
        annot["x"]  # ty: ignore[invalid-argument-type]
    assert annot[2:3].serialize() == "P"


def test_batch_result_equality_and_hashability():
    [result] = pt.batch("mass", ["PEPTIDE"])
    assert result == pt.batch("mass", ["PEPTIDE"])[0]
    assert isinstance(hash(result), int)
    [annotated] = pt.batch("mass", [pt.parse("PEPTIDE")])
    with pytest.raises(TypeError):
        hash(annotated)
