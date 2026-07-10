"""Tests for peptacular.annotation.utils - low-level fragment mass/composition helpers."""

from collections import Counter

import pytest
from tacular import ELEMENT_LOOKUP, IonType

import peptacular as pt
from peptacular.annotation.cached_comps import DeltaInfo, IsotopeInfo
from peptacular.annotation.mod import Mods
from peptacular.annotation.utils import (
    adjust_comp,
    adjust_mass_mz,
    can_fragment_sequence,
    comp_frag,
    cumsum,
    process_losses,
)
from peptacular.proforma_components.comps import ChargedFormula, GlobalChargeCarrier

H_ELEMENT = ELEMENT_LOOKUP["H"]
C_ELEMENT = ELEMENT_LOOKUP["C"]
O_ELEMENT = ELEMENT_LOOKUP["O"]
NA_ELEMENT = ELEMENT_LOOKUP["Na"]


@pytest.fixture
def no_isotope() -> IsotopeInfo:
    return IsotopeInfo.from_input(None)


@pytest.fixture
def no_delta() -> DeltaInfo:
    return DeltaInfo.from_input(None)


@pytest.fixture
def protonated_charge() -> Mods[GlobalChargeCarrier]:
    """Default +2 charge, protonated (no explicit adduct)."""
    return pt.parse("PEPTIDE").set_charge(2).charge_adducts


@pytest.fixture
def zero_charge() -> Mods[GlobalChargeCarrier]:
    """No charge carriers at all."""
    return pt.parse("PEPTIDE").set_charge(0).charge_adducts


@pytest.fixture
def sodium_adduct_charge() -> Mods[GlobalChargeCarrier]:
    """A non-protonated (sodium) charge carrier."""
    return pt.parse("PEPTIDE/[Na:z+1]").charge_adducts


@pytest.fixture
def depleting_adduct_charge() -> Mods[GlobalChargeCarrier]:
    """A charge carrier whose formula removes more of an element than a small base
    composition can supply, so it should - in principle - drive an element count
    negative once merged into the base composition."""
    return pt.parse("PEPTIDE/[H-2Na1:z+1]").charge_adducts


class TestAdjustMassMz:
    """TESTS FOR: adjust_mass_mz"""

    def test_counter_base_sums_element_masses(self, protonated_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo) -> None:
        comp = Counter({H_ELEMENT: 10, C_ELEMENT: 5})
        expected_base_mass = H_ELEMENT.get_mass(monoisotopic=True) * 10 + C_ELEMENT.get_mass(monoisotopic=True) * 5

        fragment = adjust_mass_mz(
            base=comp,
            charge=protonated_charge,
            ion_type=IonType.PRECURSOR,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            position=None,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
        )

        direct = adjust_mass_mz(
            base=expected_base_mass,
            charge=protonated_charge,
            ion_type=IonType.PRECURSOR,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            position=None,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
        )

        assert fragment.mass == pytest.approx(direct.mass)

    def test_float_base_used_directly(self, zero_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo) -> None:
        fragment = adjust_mass_mz(
            base=100.0,
            charge=zero_charge,
            ion_type=IonType.B,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            position=3,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
        )
        assert fragment.charge_state == 0
        assert fragment.position == 3

    def test_protonated_charge_has_no_explicit_adducts(
        self, protonated_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo
    ) -> None:
        fragment = adjust_mass_mz(
            base=100.0,
            charge=protonated_charge,
            ion_type=IonType.PRECURSOR,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            position=None,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
        )
        assert fragment.is_protonated is True
        assert fragment._charge_adducts is None

    def test_non_protonated_charge_reports_explicit_adducts(
        self, sodium_adduct_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo
    ) -> None:
        fragment = adjust_mass_mz(
            base=100.0,
            charge=sodium_adduct_charge,
            ion_type=IonType.PRECURSOR,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            position=None,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
        )
        assert fragment.is_protonated is False
        assert fragment._charge_adducts is not None


class TestAdjustComp:
    """TESTS FOR: adjust_comp"""

    def test_inplace_true_mutates_input(self, zero_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo) -> None:
        comp = Counter({C_ELEMENT: 5})
        adjust_comp(
            base_comp=comp,
            charge=zero_charge,
            ion_type=IonType.B,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
            position=None,
            inplace=True,
        )
        # ion B composition is empty, no charge mods, so the input counter is mutated
        # in-place (same object identity used internally) - it should still reflect
        # the original element counts since nothing was added or removed.
        assert comp[C_ELEMENT] == 5

    def test_inplace_false_leaves_original_input_untouched(
        self, sodium_adduct_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo
    ) -> None:
        comp = Counter({C_ELEMENT: 5})
        original = comp.copy()

        fragment = adjust_comp(
            base_comp=comp,
            charge=sodium_adduct_charge,
            ion_type=IonType.B,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
            position=None,
            inplace=False,
        )

        assert comp == original
        assert fragment.composition is not None
        assert fragment.composition[NA_ELEMENT] == 1

    def test_protonated_charge_has_no_explicit_adducts(
        self, protonated_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo
    ) -> None:
        fragment = adjust_comp(
            base_comp=Counter({C_ELEMENT: 5}),
            charge=protonated_charge,
            ion_type=IonType.B,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
            position=None,
            inplace=False,
        )
        assert fragment.is_protonated is True
        assert fragment._charge_adducts is None

    def test_non_protonated_charge_reports_explicit_adducts(
        self, sodium_adduct_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo
    ) -> None:
        fragment = adjust_comp(
            base_comp=Counter({C_ELEMENT: 5}),
            charge=sodium_adduct_charge,
            ion_type=IonType.B,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
            position=None,
            inplace=False,
        )
        assert fragment.is_protonated is False
        assert fragment._charge_adducts is not None

    def test_isotope_map_replaces_original_element_with_substitute(
        self, zero_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo
    ) -> None:
        """Global isotope substitutions (e.g. <13C>) swap every occurrence of the
        original element for the replacement element in the final composition.
        Uses two mapped elements so the isotope_map loop iterates more than once."""
        carbon_12 = ELEMENT_LOOKUP["C"]
        carbon_13 = ELEMENT_LOOKUP["13C"]
        nitrogen_14 = ELEMENT_LOOKUP["N"]
        nitrogen_15 = ELEMENT_LOOKUP["15N"]
        oxygen_16 = ELEMENT_LOOKUP["O"]
        oxygen_18 = ELEMENT_LOOKUP["18O"]

        # oxygen is absent from base_comp, so its entry in isotope_map exercises the
        # "original element not present" branch, while carbon/nitrogen exercise the
        # "original element present and gets replaced" branch.
        fragment = adjust_comp(
            base_comp=Counter({carbon_12: 5, nitrogen_14: 2}),
            charge=zero_charge,
            ion_type=IonType.B,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
            position=None,
            inplace=False,
            isotope_map={carbon_12: carbon_13, nitrogen_14: nitrogen_15, oxygen_16: oxygen_18},
        )

        assert fragment.composition is not None
        assert carbon_12 not in fragment.composition
        assert nitrogen_14 not in fragment.composition
        assert oxygen_18 not in fragment.composition
        assert fragment.composition[carbon_13] == 5
        assert fragment.composition[nitrogen_15] == 2

    def test_depleting_adduct_raises_when_it_removes_more_atoms_than_present(
        self, depleting_adduct_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo
    ) -> None:
        """A charge adduct formula of 'H-2Na1' removes 2 H atoms and adds 1 Na atom.

        Regression test: ``adjust_comp`` used to merge element counts via
        ``base_comp += mod.get_composition()``, which relies on ``Counter.__iadd__``
        silently discarding any element whose *resulting* count is <= 0 instead of
        keeping it negative. That meant a charge adduct depleting more atoms of an
        element than the base composition has was never caught by the "negative
        element counts" guard -- the atoms were just dropped. Composition is now
        merged element-by-element, so this case correctly raises.
        """
        comp = Counter({C_ELEMENT: 5})  # no H atoms present at all
        with pytest.raises(ValueError, match="Negative element counts"):
            adjust_comp(
                base_comp=comp,
                charge=depleting_adduct_charge,
                ion_type=IonType.B,
                monoisotopic=True,
                isotope=no_isotope,
                delta=no_delta,
                parent_sequence="PEPTIDE",
                parent_sequence_length=7,
                position=None,
                inplace=False,
            )


class TestCompFrag:
    """TESTS FOR: comp_frag"""

    def test_matches_manual_mass_via_adjust_mass_mz(self, protonated_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo) -> None:
        comp = Counter({H_ELEMENT: 10, C_ELEMENT: 5, O_ELEMENT: 2})

        via_comp_frag = comp_frag(
            comp=comp,
            charge=protonated_charge,
            ion_type=IonType.PRECURSOR,
            monoisotopic=True,
            isotopes=no_isotope,
            deltas=no_delta,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
            position=None,
        )

        manual_mass = sum(elem.get_mass(monoisotopic=True) * count for elem, count in comp.items())
        via_adjust_mass_mz = adjust_mass_mz(
            base=manual_mass,
            charge=protonated_charge,
            ion_type=IonType.PRECURSOR,
            monoisotopic=True,
            isotope=no_isotope,
            delta=no_delta,
            position=None,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
        )

        assert via_comp_frag.mass == pytest.approx(via_adjust_mass_mz.mass)

    def test_average_mass_uses_average_element_masses(self, zero_charge: Mods[GlobalChargeCarrier], no_isotope: IsotopeInfo, no_delta: DeltaInfo) -> None:
        comp = Counter({C_ELEMENT: 1})
        fragment = comp_frag(
            comp=comp,
            charge=zero_charge,
            ion_type=IonType.B,
            monoisotopic=False,
            isotopes=no_isotope,
            deltas=no_delta,
            parent_sequence="PEPTIDE",
            parent_sequence_length=7,
            position=None,
        )
        assert fragment.mass == pytest.approx(C_ELEMENT.get_mass(monoisotopic=False))


class TestProcessLosses:
    """TESTS FOR: process_losses"""

    def test_known_loss_name_resolves_to_formula(self) -> None:
        result = process_losses("Water")
        assert len(result) == 1
        (loss, count) = next(iter(result.items()))
        assert isinstance(loss, ChargedFormula)
        assert count == 1

    def test_charged_formula_input_returned_with_count_one(self) -> None:
        formula = ChargedFormula.from_string("H2O", require_formula_prefix=False)
        result = process_losses(formula)
        assert result == {formula: 1}

    def test_float_input_returned_with_count_one(self) -> None:
        result = process_losses(18.0105)
        assert result == {18.0105: 1}

    def test_non_dict_non_supported_type_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="Invalid losses type"):
            process_losses(object())  # type: ignore[arg-type]

    def test_dict_with_known_loss_name_key(self) -> None:
        result = process_losses({"Water": 2})
        assert sum(result.values()) == 2

    def test_dict_with_unknown_str_key_falls_back_to_formula_parsing(self) -> None:
        # "CH4" is a valid neutral formula but is not a registered named/formula
        # loss in NEUTRAL_DELTA_LOOKUP, so it exercises the KeyError fallback path.
        result = process_losses({"CH4": 3})
        (loss, count) = next(iter(result.items()))
        assert isinstance(loss, ChargedFormula)
        assert count == 3

    def test_dict_with_unknown_str_key_and_charge_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Loss formula cannot have charge"):
            process_losses({"CH4:z+1": 1})

    def test_dict_with_charged_formula_key(self) -> None:
        formula = ChargedFormula.from_string("H2O", require_formula_prefix=False)
        result = process_losses({formula: 4})
        assert result[formula] == 4

    def test_dict_with_float_key(self) -> None:
        result = process_losses({18.0: 2})
        assert result[18.0] == 2

    def test_dict_accumulates_counts_for_repeated_equivalent_keys(self) -> None:
        formula = ChargedFormula.from_string("H2O", require_formula_prefix=False)
        result = process_losses({formula: 1, "Water": 1})
        # "Water" resolves to the same neutral-loss formula as an explicit H2O ChargedFormula
        assert sum(result.values()) == 2

    def test_dict_with_unsupported_key_type_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="Invalid key type for loss"):
            process_losses({5: 1})  # type: ignore[dict-item]


class TestCumsum:
    """TESTS FOR: cumsum"""

    def test_empty_sequence_returns_empty_list(self) -> None:
        assert cumsum([]) == []

    def test_numeric_forward_cumulative_sum(self) -> None:
        assert cumsum([1.0, 2.0, 3.0]) == [1.0, 3.0, 6.0]

    def test_numeric_reverse_cumulative_sum(self) -> None:
        assert cumsum([1.0, 2.0, 3.0], reverse=True) == [3.0, 5.0, 6.0]

    def test_int_sequence_is_treated_as_numeric(self) -> None:
        assert cumsum([1, 2, 3]) == [1.0, 3.0, 6.0]

    def test_counter_forward_cumulative_sum(self) -> None:
        first = Counter({H_ELEMENT: 1})
        second = Counter({C_ELEMENT: 1})

        result = cumsum([first, second])

        assert result[0] == Counter({H_ELEMENT: 1})
        assert result[1] == Counter({H_ELEMENT: 1, C_ELEMENT: 1})

    def test_counter_reverse_cumulative_sum(self) -> None:
        first = Counter({H_ELEMENT: 1})
        second = Counter({C_ELEMENT: 1})

        result = cumsum([first, second], reverse=True)

        assert result[0] == Counter({C_ELEMENT: 1})
        assert result[1] == Counter({H_ELEMENT: 1, C_ELEMENT: 1})

    def test_invalid_element_type_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="cumsum expects sequence of float or Counter"):
            cumsum(["a", "b"])  # type: ignore[list-item]


class TestCanFragmentSequence:
    """TESTS FOR: can_fragment_sequence"""

    def test_string_ion_type_is_converted(self) -> None:
        assert can_fragment_sequence("PEPTIDE", "b") == IonType.B

    def test_ion_type_not_in_fragment_rules_is_returned_unchanged(self) -> None:
        assert can_fragment_sequence("PEPTIDE", IonType.Y) == IonType.Y

    # -- D ion: end position, excludes G/A/P/I/T, maps V -> D_VALINE --

    def test_d_ion_raises_on_excluded_terminal_residue(self) -> None:
        with pytest.raises(ValueError, match="D fragments cannot be produced"):
            can_fragment_sequence("PEPTIDG", IonType.D)

    def test_d_ion_unmapped_terminal_residue_returns_d(self) -> None:
        assert can_fragment_sequence("PEPTIDE", IonType.D) == IonType.D

    def test_d_ion_valine_terminus_maps_to_d_valine(self) -> None:
        assert can_fragment_sequence("PEPTIDV", IonType.D) == IonType.D_VALINE

    # -- DA ion: end position, requires I/T --

    def test_da_ion_raises_when_terminal_residue_not_required(self) -> None:
        with pytest.raises(ValueError, match="DA fragments can only be produced"):
            can_fragment_sequence("PEPTIDE", IonType.DA)

    def test_da_ion_threonine_terminus_maps_to_da_threonine(self) -> None:
        assert can_fragment_sequence("PEPTIDT", IonType.DA) == IonType.DA_THREONINE

    def test_da_ion_isoleucine_terminus_maps_to_da_isoleucine(self) -> None:
        assert can_fragment_sequence("PEPTIDI", IonType.DA) == IonType.DA_ISOLEUCINE

    # -- DB ion: start position, requires I/T --

    def test_db_ion_raises_when_initial_residue_not_required(self) -> None:
        with pytest.raises(ValueError, match="DB fragments can only be produced"):
            can_fragment_sequence("EPEPTIDE", IonType.DB)

    def test_db_ion_threonine_start_maps_to_db_threonine(self) -> None:
        assert can_fragment_sequence("TPEPTIDE", IonType.DB) == IonType.DB_THREONINE

    def test_db_ion_isoleucine_start_maps_to_db_isoleucine(self) -> None:
        assert can_fragment_sequence("IPEPTIDE", IonType.DB) == IonType.DB_ISOLEUCINE

    # -- Leaf ion type with a required set but no specific_map (return path at end) --

    def test_d_valine_ion_raises_when_terminus_is_not_valine(self) -> None:
        with pytest.raises(ValueError, match="D_VALINE fragments can only be produced"):
            can_fragment_sequence("PEPTIDE", IonType.D_VALINE)

    def test_d_valine_ion_returns_unchanged_when_terminus_is_valine(self) -> None:
        assert can_fragment_sequence("PEPTIDV", IonType.D_VALINE) == IonType.D_VALINE

    # -- W ion: start position, excludes G/A/P/I/T, maps V -> W_VALINE --

    def test_w_ion_raises_on_excluded_initial_residue(self) -> None:
        with pytest.raises(ValueError, match="W fragments cannot be produced"):
            can_fragment_sequence("GPEPTIDE", IonType.W)

    def test_w_ion_unmapped_initial_residue_returns_w(self) -> None:
        assert can_fragment_sequence("EPEPTIDE", IonType.W) == IonType.W

    def test_w_ion_valine_start_maps_to_w_valine(self) -> None:
        assert can_fragment_sequence("VPEPTIDE", IonType.W) == IonType.W_VALINE

    # -- WA ion: start position, requires I/T --

    def test_wa_ion_raises_when_initial_residue_not_required(self) -> None:
        with pytest.raises(ValueError, match="WA fragments can only be produced"):
            can_fragment_sequence("EPEPTIDE", IonType.WA)

    def test_wa_ion_threonine_start_maps_to_wa_threonine(self) -> None:
        assert can_fragment_sequence("TPEPTIDE", IonType.WA) == IonType.WA_THREONINE

    def test_wa_ion_isoleucine_start_maps_to_wa_isoleucine(self) -> None:
        assert can_fragment_sequence("IPEPTIDE", IonType.WA) == IonType.WA_ISOLEUCINE

    # -- WB ion: end position, requires I/T --

    def test_wb_ion_raises_when_terminal_residue_not_required(self) -> None:
        with pytest.raises(ValueError, match="WB fragments can only be produced"):
            can_fragment_sequence("PEPTIDE", IonType.WB)

    def test_wb_ion_threonine_terminus_maps_to_wb_threonine(self) -> None:
        assert can_fragment_sequence("PEPTIDT", IonType.WB) == IonType.WB_THREONINE

    def test_wb_ion_isoleucine_terminus_maps_to_wb_isoleucine(self) -> None:
        assert can_fragment_sequence("PEPTIDI", IonType.WB) == IonType.WB_ISOLEUCINE
