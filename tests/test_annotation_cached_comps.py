"""Coverage for `peptacular.annotation.cached_comps`.

`IsotopeInfo` and `DeltaInfo` are exercised through the public API
(`ProFormaAnnotation.mass/comp`) wherever possible. `ChargeCarrierInfo`,
`handle_charge_input`, and `get_charge_adducts` are not currently wired into
any public code path (the only call site in annotation.py is commented out),
so they are exercised directly against `cached_comps` — matching the existing
pattern in tests/test_cache_and_frag.py, which already imports `ChargeCarrierInfo`
and `DeltaInfo` directly from this module.
"""

from collections import Counter

import pytest
from tacular import ELEMENT_LOOKUP

import peptacular as pt
from peptacular.annotation.cached_comps import (
    ChargeCarrierInfo,
    DeltaInfo,
    IsotopeInfo,
    _handle_delta_input,
    get_isotopes,
    get_losses,
    handle_charge_input,
)
from peptacular.proforma_components import ChargedFormula


class TestIsotopeInfoComposition:
    """IsotopeInfo.composition and average_mass_delta, reached via ProFormaAnnotation."""

    def test_composition_matches_mass_adjustment(self):
        a = pt.parse("PEPTIDE")
        base_comp = a.comp()
        iso_comp = a.comp(isotopes={"13C": 2})
        c12 = ELEMENT_LOOKUP["C"]
        c13 = ELEMENT_LOOKUP["13C"]
        assert iso_comp[c12] == base_comp[c12] - 2
        assert iso_comp[c13] == 2

    def test_average_mass_delta_differs_from_monoisotopic(self):
        a = pt.parse("PEPTIDE")
        mono = a.mass(isotopes={"13C": 2}, monoisotopic=True)
        avg = a.mass(isotopes={"13C": 2}, monoisotopic=False)
        assert mono != avg

    def test_average_mass_isotope_shift_is_positive(self):
        a = pt.parse("PEPTIDE")
        base_avg = a.mass(monoisotopic=False)
        shifted_avg = a.mass(isotopes={"13C": 2}, monoisotopic=False)
        assert shifted_avg > base_avg


class TestIsotopeInfoDirect:
    """Direct checks of IsotopeInfo/get_isotopes semantics not reachable via .mass()/.comp() alone."""

    def test_int_and_equivalent_dict_produce_same_composition(self):
        int_info = get_isotopes(2)
        dict_info = get_isotopes({"13C": 2})
        assert int_info.composition == dict_info.composition

    def test_zero_isotopes_has_empty_composition(self):
        info = get_isotopes(0)
        assert info.composition == Counter()

    def test_none_normalizes_to_zero(self):
        assert get_isotopes(None).composition == get_isotopes(0).composition

    def test_negative_isotope_count_raises(self):
        with pytest.raises(ValueError, match="cannot be negative"):
            get_isotopes(-1)

    def test_invalid_isotope_type_raises(self):
        with pytest.raises(TypeError, match="Invalid isotope type"):
            get_isotopes(3.5)  # type: ignore[arg-type]

    def test_adjust_composition_raises_on_over_removal(self):
        info = get_isotopes({"13C": 100})
        comp: Counter = Counter({ELEMENT_LOOKUP["C"]: 1})
        with pytest.raises(ValueError, match="negative element counts"):
            info.adjust_composition(comp)

    def test_to_fragment_mapping_single_13c_returns_bare_count(self):
        info = get_isotopes(3)
        assert info.to_fragment_mapping == 3

    def test_to_fragment_mapping_dict_isotope_returns_dict(self):
        info = get_isotopes({"15N": 1})
        assert info.to_fragment_mapping == {"15N": 1}

    def test_to_fragment_mapping_empty_returns_empty_dict(self):
        assert get_isotopes(0).to_fragment_mapping == {}

    def test_from_input_matches_get_isotopes(self):
        assert IsotopeInfo.from_input({"13C": 1}).composition == get_isotopes({"13C": 1}).composition

    def test_dict_with_element_info_key_matches_string_key(self):
        # dict keys may be ElementInfo objects directly, not just symbol strings.
        c13 = ELEMENT_LOOKUP["13C"]
        by_object = get_isotopes({c13: 2})
        by_string = get_isotopes({"13C": 2})
        assert by_object.composition == by_string.composition


class TestChargeCarrierInfo:
    """ChargeCarrierInfo has no current public call site (annotation.py:3303 is commented
    out) — exercised directly against cached_comps, following tests/test_cache_and_frag.py."""

    def test_from_input_int_charge(self):
        info = ChargeCarrierInfo.from_input(2)
        assert info.charge == 2

    def test_get_mass_monoisotopic_vs_average_differ(self):
        info = ChargeCarrierInfo.from_input(2)
        assert info.get_mass(monoisotopic=True) == info.monoisotopic_mass
        assert info.get_mass(monoisotopic=False) == info.average_mass
        assert info.monoisotopic_mass != info.average_mass

    def test_composition_counts_protons_as_hydrogen(self):
        info = ChargeCarrierInfo.from_input(2)
        assert info.composition == Counter({ELEMENT_LOOKUP["H"]: 2})

    def test_adjust_composition_adds_counts(self):
        info = ChargeCarrierInfo.from_input(2)
        comp: Counter = Counter()
        info.adjust_composition(comp)
        assert comp[ELEMENT_LOOKUP["H"]] == 2

    def test_adjust_composition_raises_on_negative_result(self):
        info = ChargeCarrierInfo.from_input("H-1:z+1")
        comp: Counter = Counter()
        with pytest.raises(ValueError, match="negative element counts"):
            info.adjust_composition(comp)

    def test_to_fragment_mapping_bare_proton_is_none(self):
        info = ChargeCarrierInfo.from_input(2)
        assert info.to_fragment_mapping is None

    def test_to_fragment_mapping_non_proton_returns_dict(self):
        info = ChargeCarrierInfo.from_input("Na:z+1")
        assert info.to_fragment_mapping == {"Na": 1}

    def test_to_fragment_mapping_empty_is_none(self):
        info = ChargeCarrierInfo.from_input(None)
        assert info.to_fragment_mapping is None

    def test_to_explicit_fragment_mapping_includes_protons(self):
        info = ChargeCarrierInfo.from_input(2)
        assert info.to_explicit_fragment_mapping == {"H": 2}

    def test_to_explicit_fragment_mapping_empty_is_empty_dict(self):
        info = ChargeCarrierInfo.from_input(None)
        assert info.to_explicit_fragment_mapping == {}

    def test_to_proforma_charge_bare_proton_returns_int(self):
        info = ChargeCarrierInfo.from_input(2)
        assert info.to_proforma_charge == 2

    def test_to_proforma_charge_non_proton_returns_dict(self):
        info = ChargeCarrierInfo.from_input("Na:z+1")
        assert info.to_proforma_charge == {"Na": 1}

    def test_to_proforma_charge_empty_is_none(self):
        assert ChargeCarrierInfo.from_input(None).to_proforma_charge is None

    def test_to_tuple_or_int_bare_proton_returns_int(self):
        info = ChargeCarrierInfo.from_input(2)
        assert info.to_tuple_or_int == 2

    def test_to_tuple_or_int_non_proton_returns_tuple(self):
        info = ChargeCarrierInfo.from_input("Na:z+1")
        assert info.to_tuple_or_int == ("Na",)

    def test_to_tuple_or_int_empty_is_none(self):
        assert ChargeCarrierInfo.from_input(None).to_tuple_or_int is None

    def test_from_input_is_cached_singleton(self):
        # Sorted-by-serialization normalization means these two inputs should hit the
        # same cached instance.
        a = ChargeCarrierInfo.from_input(2)
        b = ChargeCarrierInfo.from_input(2)
        assert a is b


class TestHandleChargeInput:
    """handle_charge_input is only reachable via ChargeCarrierInfo.from_input, which has
    no public caller — exercised directly, following the ChargeCarrierInfo precedent."""

    def test_int_produces_proton_adduct(self):
        adducts = handle_charge_input(2)
        assert len(adducts) == 1
        assert adducts[0].get_charge() == 2

    def test_str_parses_single_adduct(self):
        adducts = handle_charge_input("H:z+1")
        assert len(adducts) == 1

    def test_global_charge_carrier_passthrough(self):
        (carrier,) = handle_charge_input(2)
        assert handle_charge_input(carrier) == (carrier,)

    def test_tuple_flattens_mixed_entries(self):
        adducts = handle_charge_input((2, "H:z+1"))
        assert len(adducts) == 2

    def test_none_returns_empty_tuple(self):
        assert handle_charge_input(None) == ()

    def test_list_of_strings_builds_adducts(self):
        adducts = handle_charge_input(["H:z+1", "Na:z+1"])
        assert len(adducts) == 2

    def test_invalid_type_raises_type_error(self):
        with pytest.raises(TypeError, match="Invalid charge type"):
            handle_charge_input(3.5)  # type: ignore[arg-type]


class TestDeltaInfoComposition:
    """DeltaInfo.composition/adjust_composition, reached through ProFormaAnnotation.comp."""

    def test_composition_reflects_named_delta(self):
        a = pt.parse("PEPTIDE")
        base_comp = a.comp()
        shifted_comp = a.comp(deltas="H2O")
        # "H2O" in NEUTRAL_DELTA_LOOKUP represents *loss* of water.
        assert shifted_comp[ELEMENT_LOOKUP["H"]] == base_comp[ELEMENT_LOOKUP["H"]] - 2
        assert shifted_comp[ELEMENT_LOOKUP["O"]] == base_comp[ELEMENT_LOOKUP["O"]] - 1

    def test_composition_raises_when_deltas_include_floats(self):
        info = DeltaInfo.from_input(1.5)
        with pytest.raises(ValueError, match="Cannot get composition"):
            info.composition

    def test_adjust_composition_raises_when_deltas_include_floats(self):
        info = DeltaInfo.from_input(1.5)
        with pytest.raises(ValueError, match="Cannot adjust composition"):
            info.adjust_composition(Counter())

    def test_adjust_composition_raises_on_negative_result(self):
        # "H2O" is a *loss* of water (dict_composition={"H": -2, "O": -1}), so applying it
        # once to an empty (all-zero) composition drives H/O negative.
        info = DeltaInfo.from_input({"H2O": 1})
        with pytest.raises(ValueError, match="negative element counts"):
            info.adjust_composition(Counter())

    def test_average_mass_delta_matches_public_mass_call(self):
        a = pt.parse("PEPTIDE")
        avg = a.mass(deltas="H2O", monoisotopic=False)
        base_avg = a.mass(monoisotopic=False)
        info = DeltaInfo.from_input("H2O")
        assert avg == pytest.approx(base_avg + info.average_mass_delta)

    def test_composition_property_sums_charged_formula_keys(self):
        info = DeltaInfo.from_input({"H2O": 1})
        comp = info.composition
        assert comp[ELEMENT_LOOKUP["H"]] == -2
        assert comp[ELEMENT_LOOKUP["O"]] == -1

    def test_mass_delta_with_float_key_monoisotopic(self):
        info = DeltaInfo.from_input(2.5)
        assert info.monoisotopic_mass_delta == pytest.approx(2.5)

    def test_mass_delta_with_float_key_average(self):
        info = DeltaInfo.from_input(2.5)
        assert info.average_mass_delta == pytest.approx(2.5)

    def test_get_mass_delta_dispatches_on_monoisotopic_flag(self):
        info = DeltaInfo.from_input("H2O")
        assert info.get_mass_delta(monoisotopic=True) == info.monoisotopic_mass_delta
        assert info.get_mass_delta(monoisotopic=False) == info.average_mass_delta


class TestDeltaInfoFromInput:
    def test_int_delta_key_is_normalized_to_float(self):
        info = DeltaInfo.from_input(5)
        assert info.deltas == {5.0: 1}

    def test_raw_formula_string_not_in_named_lookup(self):
        # "C2H4" isn't a NEUTRAL_DELTA_LOOKUP name, so from_input falls back to parsing it
        # as a bare chemical formula.
        info = DeltaInfo.from_input("C2H4")
        assert info.composition[ELEMENT_LOOKUP["C"]] == 2
        assert info.composition[ELEMENT_LOOKUP["H"]] == 4

    def test_charged_formula_delta_raises(self):
        charged = ChargedFormula.from_string("H:z+1", require_formula_prefix=False)
        with pytest.raises(ValueError, match="must be neutral"):
            DeltaInfo.from_input(charged)


class TestDeltaInfoAdd:
    def test_add_combines_matching_keys(self):
        one = DeltaInfo.from_input({"H2O": 1})
        two = DeltaInfo.from_input({"H2O": 1})
        combined = one + two
        assert combined.deltas == DeltaInfo.from_input({"H2O": 2}).deltas

    def test_add_drops_keys_that_cancel_to_zero(self):
        one = DeltaInfo.from_input({"H2O": 1})
        minus_one = DeltaInfo.from_input({"H2O": -1})
        assert (one + minus_one).deltas == {}


class TestGetLosses:
    def test_nonempty_counter_returns_single_charged_formula(self):
        h = ELEMENT_LOOKUP["H"]
        result = get_losses(Counter({h: 2}))
        assert result is not None
        assert list(result.values()) == [1]

    def test_empty_counter_returns_none(self):
        assert get_losses(Counter()) is None


class TestDeltaInfoArithmetic:
    """__add__ is covered elsewhere; __sub__/__neg__/__mul__/__rmul__/__str__ are not."""

    def test_sub_removes_matching_counts(self):
        two = DeltaInfo.from_input({"H2O": 2})
        one = DeltaInfo.from_input({"H2O": 1})
        result = two - one
        assert result.deltas == one.deltas

    def test_sub_drops_zeroed_keys(self):
        one = DeltaInfo.from_input({"H2O": 1})
        result = one - one
        assert result.deltas == {}

    def test_neg_flips_sign(self):
        one = DeltaInfo.from_input({"H2O": 1})
        negated = -one
        assert negated.monoisotopic_mass_delta == pytest.approx(-one.monoisotopic_mass_delta)

    def test_mul_scales_counts(self):
        one = DeltaInfo.from_input({"H2O": 1})
        tripled = one * 3
        assert tripled.monoisotopic_mass_delta == pytest.approx(3 * one.monoisotopic_mass_delta)

    def test_rmul_matches_mul(self):
        one = DeltaInfo.from_input({"H2O": 1})
        assert (3 * one).deltas == (one * 3).deltas

    def test_mul_by_zero_drops_key(self):
        one = DeltaInfo.from_input({"H2O": 1})
        assert (one * 0).deltas == {}

    def test_str_includes_formula_prefix_and_count(self):
        info = DeltaInfo.from_input({"H2O": 2})
        assert str(info).startswith("Formula:")
        assert "x2" in str(info)

    def test_str_float_delta_has_no_multiplier_suffix(self):
        info = DeltaInfo.from_input(1.5)
        assert str(info) == "1.5"

    def test_str_empty_delta_is_no_delta(self):
        assert str(DeltaInfo.from_input(None)) == "No Delta"


class TestDeltaInfoFragmentMapping:
    def test_empty_deltas_returns_none(self):
        assert DeltaInfo.from_input(None).to_fragment_mapping is None

    def test_float_key_preserved_as_is(self):
        mapping = DeltaInfo.from_input(1.5).to_fragment_mapping
        assert mapping == {1.5: 1}

    def test_formula_key_serialized_without_prefix(self):
        mapping = DeltaInfo.from_input({"H2O": 1}).to_fragment_mapping
        assert mapping is not None
        assert "H2O" not in mapping  # "H2O" is stored as its (negative) elemental formula
        assert any(v == 1 for v in mapping.values())


class TestHandleDeltaInputTypeError:
    def test_invalid_type_raises(self):
        with pytest.raises(TypeError, match="Invalid delta type"):
            _handle_delta_input([1, 2])  # type: ignore[arg-type]

    def test_none_key_in_dict_normalizes_to_empty(self):
        # DeltaInfo.from_input({None: 1}) recurses _handle_delta_input(None) for the key,
        # which is the only public-reachable path that hits the deltas-is-None base case.
        info = DeltaInfo.from_input({None: 1})  # type: ignore[dict-item]
        assert info.deltas == {}
