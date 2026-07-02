"""Regression tests for glycan-occurrence, cache-mutation, and neutral-mass bugs."""

from collections import Counter

import peptacular as pt
from peptacular.annotation.cached_comps import ChargeCarrierInfo, DeltaInfo, get_losses


def _comp_mass(seq: str) -> float:
    a = pt.parse(seq)
    return sum(e.get_mass(monoisotopic=True) * n for e, n in a.comp().items())


class TestGlycanOccurrence:
    """GlycanComponent.get_composition must multiply by occurrence (match get_mass)."""

    def test_composition_matches_mass_for_counts_gt_one(self):
        for seq in ["PEPT[Glycan:Hex3]IDE", "PEPT[Glycan:HexNAc2Hex3]IDE", "PEPT[Glycan:Hex1]IDE"]:
            assert abs(_comp_mass(seq) - pt.parse(seq).mass()) < 1e-3, seq

    def test_hex3_is_three_hex(self):
        one = pt.parse("PEPT[Glycan:Hex1]IDE").mass() - pt.mass("PEPTIDE")
        three = pt.parse("PEPT[Glycan:Hex3]IDE").mass() - pt.mass("PEPTIDE")
        assert abs(three - 3 * one) < 1e-6


class TestCacheMutationIsolation:
    """Cached singleton accessors must not hand out shared mutable containers."""

    def test_charge_carrier_composition_isolated(self):
        from tacular import ELEMENT_LOOKUP

        h = ELEMENT_LOOKUP["H"]
        c = ChargeCarrierInfo.from_input(2).composition
        c[h] = c.get(h, 0) + 100
        assert ChargeCarrierInfo.from_input(2).composition.get(h) == 2

    def test_delta_fragment_mapping_isolated(self):
        m = DeltaInfo.from_input({"H2O": 1}).to_fragment_mapping
        m["INJECTED"] = 999  # type: ignore[index]
        assert "INJECTED" not in DeltaInfo.from_input({"H2O": 1}).to_fragment_mapping

    def test_get_losses_isolated(self):
        from tacular import ELEMENT_LOOKUP

        h = ELEMENT_LOOKUP["H"]
        losses = get_losses(Counter({h: 2}))
        losses["INJECTED"] = 999  # type: ignore[index]
        assert "INJECTED" not in get_losses(Counter({h: 2}))


class TestMzPAFChargeSign:
    """mzPAF charge component must be a bare magnitude, never signed (spec section 4.8:
    "The charge state component in the peak annotation MUST NOT include the minus sign")."""

    def test_positive_charge_omits_one(self):
        a = pt.parse("PEPTIDE")
        y3 = [f for f in a.fragment(ion_types=["y"], charges=[1]) if f.position == 3][0]
        assert y3.serialize(format="mzpaf") == "y3{IDE}"

    def test_positive_multi_charge_bare(self):
        a = pt.parse("PEPTIDE")
        y3 = [f for f in a.fragment(ion_types=["y"], charges=[2]) if f.position == 3][0]
        assert y3.serialize(format="mzpaf") == "y3{IDE}^2"

    def test_negative_charge_no_minus_sign(self):
        a = pt.parse("PEPTIDE")
        for z, expected in ((-1, "y3{IDE}^1"), (-2, "y3{IDE}^2")):
            y3 = [f for f in a.fragment(ion_types=["y"], charges=[z]) if f.position == 3][0]
            assert y3.serialize(format="mzpaf") == expected


class TestFragmentNeutralMass:
    """Fragment.neutral_mass must be charge-invariant (electron correction undone)."""

    def test_neutral_mass_constant_across_charges(self):
        a = pt.parse("PEPTIDE")
        vals = []
        for z in (1, 2, 3, -1):
            y3 = [f for f in a.fragment(ion_types=["y"], charges=[z]) if f.position == 3][0]
            vals.append(y3.neutral_mass)
        assert max(vals) - min(vals) < 1e-9
