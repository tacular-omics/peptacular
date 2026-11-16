"""Activation-method presets: the keyword arguments fragment() already takes."""

import pytest

import peptacular as pt

PEPTIDE = "ELVISLIVESK"  # 11 residues

# Measured with pt.fragment(PEPTIDE, **preset): total ions, then ions with no neutral loss,
# for a 2+ precursor (fragment charges (1,)) and a 3+ precursor (fragment charges (1, 2)).
EXPECTED = {
    "cid": {"series": {"b": 23, "y": 32}, "total": 55, "no_loss": 22},
    "hcd": {"series": {"b": 23, "y": 32}, "total": 55, "no_loss": 22},
    "etd": {"series": {"c": 23, "z.": 32}, "total": 55, "no_loss": 22},
    "ecd": {"series": {"c": 23, "z.": 32}, "total": 55, "no_loss": 22},
    "ethcd": {"series": {"b": 23, "y": 32, "c": 23, "z.": 32}, "total": 110, "no_loss": 44},
    "ead": {"series": {"a": 23, "b": 23, "c": 23, "x": 32, "y": 32, "z.": 32}, "total": 165, "no_loss": 66},
    "uvpd": {"series": {"a": 22, "b": 22, "c": 22, "x": 21, "y": 21, "z": 21}, "total": 129, "no_loss": 66},
}

PRESET_NAMES = sorted(pt.FRAGMENTATION_PRESETS)


def _counts(fragments):
    series: dict[str, int] = {}
    no_loss = 0
    for fragment in fragments:
        series[fragment.ion_type.value] = series.get(fragment.ion_type.value, 0) + 1
        if not fragment.deltas:
            no_loss += 1
    return series, no_loss


class TestPresetTable:
    def test_every_name_is_covered_by_the_expected_counts(self):
        assert PRESET_NAMES == sorted(EXPECTED)

    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_every_preset_has_a_citation(self, name: str):
        assert pt.FRAGMENTATION_PRESET_CITATIONS[name]

    def test_citations_cover_exactly_the_presets(self):
        assert sorted(pt.FRAGMENTATION_PRESET_CITATIONS) == PRESET_NAMES

    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_every_preset_sets_the_same_three_keys(self, name: str):
        assert set(pt.FRAGMENTATION_PRESETS[name]) == {"ion_types", "neutral_deltas", "max_ndeltas"}

    @pytest.mark.parametrize(("alias", "target"), sorted(pt.FRAGMENTATION_PRESET_ALIASES.items()))
    def test_every_alias_resolves_to_its_preset(self, alias: str, target: str):
        assert pt.fragmentation_preset(alias) == pt.fragmentation_preset(target)

    def test_etd_holds_no_collisional_series(self):
        # The b and y ions of an electron-transfer experiment come from supplemental collisional
        # activation of the charge-reduced survivor, which is the ethcd preset.
        assert pt.FRAGMENTATION_PRESETS["etd"]["ion_types"] == ("c", "z.")
        assert set(pt.FRAGMENTATION_PRESETS["ethcd"]["ion_types"]) > {"b", "y"}


class TestPresetLookup:
    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_name_is_case_and_space_insensitive(self, name: str):
        assert pt.fragmentation_preset(f"  {name.upper()}  ") == pt.fragmentation_preset(name)

    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_no_charge_state_sets_no_charges(self, name: str):
        assert "charges" not in pt.fragmentation_preset(name)

    @pytest.mark.parametrize(("charge_state", "charges"), [(1, (1,)), (2, (1,)), (3, (1, 2)), (4, (1, 2, 3)), (-3, (1, 2))])
    def test_charge_rule_is_one_to_precursor_minus_one(self, charge_state: int, charges: tuple[int, ...]):
        assert pt.fragmentation_preset("hcd", charge_state=charge_state)["charges"] == charges

    def test_each_call_returns_a_new_dict(self):
        first = pt.fragmentation_preset("hcd")
        first["ion_types"] = ("b",)
        assert pt.fragmentation_preset("hcd")["ion_types"] == ("b", "y")

    def test_unknown_method_raises_and_names_the_valid_ones(self):
        with pytest.raises(pt.PeptacularError) as excinfo:
            pt.fragmentation_preset("mystery")
        message = str(excinfo.value)
        assert "mystery" in message
        for name in PRESET_NAMES:
            assert name in message
        for alias in pt.FRAGMENTATION_PRESET_ALIASES:
            assert alias in message


class TestPresetCounts:
    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_counts_for_a_doubly_charged_precursor(self, name: str):
        fragments = pt.fragment(PEPTIDE, **pt.fragmentation_preset(name, charge_state=2))
        series, no_loss = _counts(fragments)
        assert series == EXPECTED[name]["series"]
        assert len(fragments) == EXPECTED[name]["total"]
        assert no_loss == EXPECTED[name]["no_loss"]

    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_a_triply_charged_precursor_doubles_every_count(self, name: str):
        fragments = pt.fragment(PEPTIDE, **pt.fragmentation_preset(name, charge_state=3))
        series, no_loss = _counts(fragments)
        assert series == {ion: 2 * count for ion, count in EXPECTED[name]["series"].items()}
        assert len(fragments) == 2 * EXPECTED[name]["total"]
        assert no_loss == 2 * EXPECTED[name]["no_loss"]

    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_one_no_loss_ion_per_residue_per_series(self, name: str):
        preset = pt.fragmentation_preset(name, charge_state=2)
        _, no_loss = _counts(pt.fragment(PEPTIDE, **preset))
        assert no_loss == pt.sequence_length(PEPTIDE) * len(preset["ion_types"])


class TestOverride:
    def test_the_callers_key_wins(self):
        hcd = pt.fragmentation_preset("hcd", charge_state=2)
        overridden = {**hcd, "neutral_deltas": ("H3PO4",)}
        assert overridden["ion_types"] == hcd["ion_types"]
        assert overridden["charges"] == hcd["charges"]
        phospho = pt.fragment("ELVIS[Phospho]LIVESK", **overridden)
        assert {fragment.ion_type.value for fragment in phospho} == {"b", "y"}
        assert len(phospho) < len(pt.fragment("ELVIS[Phospho]LIVESK", **hcd))

    def test_overriding_ion_types_changes_the_series(self):
        hcd = pt.fragmentation_preset("hcd", charge_state=2)
        fragments = pt.fragment(PEPTIDE, **{**hcd, "ion_types": ("b",)})
        series, _ = _counts(fragments)
        assert series == {"b": EXPECTED["hcd"]["series"]["b"]}

    def test_writing_a_new_method_into_the_table(self):
        pt.FRAGMENTATION_PRESETS["spike"] = {"ion_types": ("b",), "neutral_deltas": (), "max_ndeltas": 1}
        try:
            assert pt.fragmentation_preset("spike", charge_state=2)["charges"] == (1,)
            assert len(pt.fragment(PEPTIDE, **pt.fragmentation_preset("spike", charge_state=2))) == 11
        finally:
            del pt.FRAGMENTATION_PRESETS["spike"]


class TestIonsAreUsable:
    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_every_ion_round_trips_through_mzpaf_on_a_modified_peptide(self, name: str):
        modified = "[Acetyl]-ELVIS[Phospho]LIVEM[Oxidation]K/2"
        fragments = pt.fragment(modified, **pt.fragmentation_preset(name, charge_state=2))
        assert fragments
        for fragment in fragments:
            assert fragment.mz > 0
            assert fragment.to_mzpaf()
            assert fragment.sequence is not None

    @pytest.mark.parametrize("name", PRESET_NAMES)
    def test_the_filtered_preset_works_with_fast_fragment(self, name: str):
        preset = pt.fragmentation_preset(name, charge_state=2)
        fast_keys = {key: value for key, value in preset.items() if key in ("ion_types", "charges", "monoisotopic")}
        fast_ions = [ion for ion in preset["ion_types"] if ion in ("a", "b", "c", "x", "y", "z")]
        if len(fast_ions) != len(preset["ion_types"]):
            # fast_fragment covers the six main series only, so a radical series is rejected.
            with pytest.raises(pt.UnsupportedOperationError):
                pt.fast_fragment(PEPTIDE, **fast_keys)
            return
        assert len(pt.fast_fragment(PEPTIDE, **fast_keys)) == len(fast_ions)

    def test_splatting_a_whole_preset_into_fast_fragment_is_a_type_error(self):
        # Documented in the module docstring: fast_fragment takes three of the preset's keys.
        with pytest.raises(TypeError):
            pt.fast_fragment(PEPTIDE, **pt.fragmentation_preset("hcd", charge_state=2))
