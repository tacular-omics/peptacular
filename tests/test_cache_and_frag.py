"""Regression tests for glycan-occurrence, cache-mutation, and neutral-mass bugs."""

from collections import Counter

import pytest

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
    """A negative mzPAF charge is written signed (``^-n``) so the label parses back to the same
    m/z, as paftacular 2.0 does. ``signed_charge=False`` writes the bare magnitude (mzPAF 1.0.1
    section 4.8)."""

    def test_positive_charge_omits_one(self):
        a = pt.parse("PEPTIDE")
        y3 = [f for f in a.fragment(ion_types=["y"], charges=[1]) if f.position == 3][0]
        assert y3.serialize(format="mzpaf") == "y3{IDE}"

    def test_positive_multi_charge_bare(self):
        a = pt.parse("PEPTIDE")
        y3 = [f for f in a.fragment(ion_types=["y"], charges=[2]) if f.position == 3][0]
        assert y3.serialize(format="mzpaf") == "y3{IDE}^2"

    def test_negative_charge_signed(self):
        a = pt.parse("PEPTIDE")
        for z, expected in ((-1, "y3{IDE}^-1"), (-2, "y3{IDE}^-2")):
            y3 = [f for f in a.fragment(ion_types=["y"], charges=[z]) if f.position == 3][0]
            assert y3.serialize(format="mzpaf") == expected
            assert y3.to_mzpaf() == expected

    def test_negative_charge_unsigned_opt_out(self):
        a = pt.parse("PEPTIDE")
        # the magnitude 1 is implicit, as in paftacular's serialize(signed_charge=False)
        for z, expected in ((-1, "y3{IDE}"), (-2, "y3{IDE}^2")):
            y3 = [f for f in a.fragment(ion_types=["y"], charges=[z]) if f.position == 3][0]
            assert y3.serialize(format="mzpaf", signed_charge=False) == expected
            assert y3.to_mzpaf(signed_charge=False) == expected

    def test_str_keeps_negative_sign(self):
        y3 = pt.parse("PEPTIDE").frag(ion_type="y", charge=-2, position=3)
        assert "charge=-2" in str(y3)


class TestFragmentNeutralMass:
    """Fragment.neutral_mass must be charge-invariant (electron correction undone)."""

    def test_neutral_mass_constant_across_charges(self):
        a = pt.parse("PEPTIDE")
        vals = []
        for z in (1, 2, 3, -1):
            y3 = [f for f in a.fragment(ion_types=["y"], charges=[z]) if f.position == 3][0]
            vals.append(y3.neutral_mass)
        assert max(vals) - min(vals) < 1e-9


class TestMzPAFLabelMass:
    """An mzPAF label must describe the same ion as the fragment it was written from."""

    CASES = (
        ("PEPTIDE", "y", -1, 3),
        ("PEPTIDE", "b", -3, 2),
        ("[Acetyl]-PEPTIDE", "i", 1, 1),
        ("[Acetyl]-PEPTIDE", "i", -1, 1),
        ("PEPTIDE-[Amidated]", "i", 1, 7),
        ("PEP[+10]TIDE", "i", 2, 3),
    )

    def test_immonium_terminal_mod_mass_matches_residue_mod(self):
        # IP[Acetyl] means P carrying Acetyl: same mass as the N-terminal acetyl immonium
        nterm = pt.parse("[Acetyl]-PEPTIDE").frag(ion_type="i", charge=1, position=1)
        residue = pt.parse("P[Acetyl]EPTIDE").frag(ion_type="i", charge=1, position=1)
        assert nterm.to_mzpaf() == residue.to_mzpaf() == "IP[Acetyl]"
        assert abs(nterm.mz - residue.mz) < 1e-9

    def test_label_mz_matches_paftacular(self):
        paf = pytest.importorskip("paftacular")
        for seq, ion, z, pos in self.CASES:
            frag = pt.parse(seq).frag(ion_type=ion, charge=z, position=pos)
            parsed = paf.parse(frag.to_mzpaf())
            parsed = parsed[0] if isinstance(parsed, list) else parsed
            assert abs(parsed.mz() - frag.mz) < 1e-6, (seq, ion, z, pos, frag.to_mzpaf())
