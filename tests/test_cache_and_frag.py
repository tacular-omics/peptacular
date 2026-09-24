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
        ("PEPTIDE", "y", "H:z-1", 3),
        ("PEPTIDE", "y", ["H:z-1", "H:z-1"], 3),
        ("PEPTIDE", "ax", 1, (2, 5)),
        ("PEPTIDE", "bx", 2, (2, 5)),
        ("[Acetyl]-TIV", "v", 1, 3),
        ("[Acetyl]-VIT", "w", 1, 3),
        ("TIV-[Amidated]", "d", 1, 3),
        ("<[Oxidation]@P>PEPTIDE", "i", 1, 1),
        ("<13C>PEPTIDE", "i", 1, 1),
        ("<13C><15N>PEPTIDE", "i", 2, 1),
        ("<D>PEPTIDE", "i", 1, 1),
        ("<13C>[Acetyl]-PEPTIDE", "i", 1, 1),
        ("<[Oxidation]@P><13C>PEPTIDE", "i", -1, 1),
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


# Every ion type whose full-length ion contains both termini, with a sequence it can be made from.
_FULL_LENGTH_CASES = (
    ("a", "TIV"),
    ("b", "TIV"),
    ("c", "TIV"),
    ("x", "TIV"),
    ("y", "TIV"),
    ("z", "TIV"),
    ("z.", "TIV"),
    ("z+H", "TIV"),
    ("c-H", "TIV"),
    ("d", "TIV"),
    ("d-valine", "TIV"),
    ("da-threonine", "VIT"),
    ("db-threonine", "VIT"),
    ("v", "TIV"),
    ("w", "VIT"),
    ("w-valine", "VIT"),
    ("wa", "TIV"),
    ("wb", "TIV"),
    ("wa-threonine", "TIV"),
    ("wb-threonine", "TIV"),
    ("p", "TIV"),
)


class TestFullLengthTerminalMods:
    """A full-length ion contains both termini, so both terminal mods add their mass. The
    satellite ions (d, v, w) dropped the terminal mod on the residue whose side chain is lost."""

    @pytest.mark.parametrize(("ion", "core"), _FULL_LENGTH_CASES)
    def test_full_length_ion_carries_both_terminal_mods(self, ion, core):
        pos = None if ion == "p" else len(core)

        def mass(seq):
            return pt.parse(seq).frag(ion_type=ion, charge=1, position=pos).mass

        base = mass(core)
        acetyl = pt.mass("[Acetyl]-G") - pt.mass("G")
        amidated = pt.mass("G-[Amidated]") - pt.mass("G")
        assert abs(mass(f"[Acetyl]-{core}") - base - acetyl) < 1e-9
        assert abs(mass(f"{core}-[Amidated]") - base - amidated) < 1e-9

    def test_partial_satellite_ions_unchanged(self):
        # v2 of TIV is IV: no N-terminus, so the N-terminal mod is not in it
        assert pt.parse("[Acetyl]-TIV").frag(ion_type="v", charge=1, position=2).mass == pt.parse("TIV").frag(ion_type="v", charge=1, position=2).mass
        # da2 of TIV is TI: no C-terminus
        assert pt.parse("TIV-[Amidated]").frag(ion_type="da", charge=1, position=2).mass == pt.parse("TIV").frag(ion_type="da", charge=1, position=2).mass

    def test_fragment_matches_frag(self):
        full = [f for f in pt.parse("[Acetyl]-PPA").fragment(ion_types=["v"], charges=[1]) if f.position == 3][0]
        assert full.mass == pt.parse("[Acetyl]-PPA").frag(ion_type="v", charge=1, position=3).mass
        assert full.mass > pt.parse("PPA").frag(ion_type="v", charge=1, position=3).mass + 42


class TestHydrideCarrier:
    """``H:z-1`` is a hydride adduct (H plus an electron), not a removed proton."""

    def test_hydride_is_not_protonated(self):
        assert pt.GlobalChargeCarrier.from_string("H:z+1").is_protonated
        assert not pt.GlobalChargeCarrier.from_string("H:z-1").is_protonated

    def test_hydride_fragment(self):
        y3 = pt.parse("PEPTIDE").frag(ion_type="y", charge="H:z-1", position=3)
        deprot = pt.parse("PEPTIDE").frag(ion_type="y", charge=-1, position=3)
        assert y3.charge_state == -1
        assert not y3.is_protonated
        assert "H:z-1" in str(y3)
        # hydride adds H + e-; deprotonation removes a proton (H - e- plus the H binding energy)
        assert abs(y3.mass - deprot.mass - 2 * 1.00782503223 - pt.constants.HYDROGEN_BINDING_MASS) < 1e-9
        assert y3.to_mzpaf() == "y3{IDE}[M+H]^-1"
        assert deprot.to_mzpaf() == "y3{IDE}^-1"

    def test_two_hydrides(self):
        for charge in (["H:z-1", "H:z-1"], "H:z-1^2"):
            y3 = pt.parse("PEPTIDE").frag(ion_type="y", charge=charge, position=3)
            assert y3.to_mzpaf() == "y3{IDE}[M+2H]^-2"

    def test_carriers_summing_to_zero_raise(self):
        y3 = pt.parse("PEPTIDE").frag(ion_type="y", charge=["H:z-1", "Na:z+1"], position=3)
        assert y3.charge_state == 0
        with pytest.raises(pt.PeptacularError, match="uncharged"):
            y3.to_mzpaf()
        with pytest.raises(pt.PeptacularError, match="uncharged"):
            y3.serialize(format="mzpaf")

    def test_fragment_series_hydride(self):
        y3 = [f for f in pt.parse("PEPTIDE").fragment(ion_types=["y"], charges=["H:z-1"]) if f.position == 3][0]
        assert y3.to_mzpaf() == "y3{IDE}[M+H]^-1"


class TestMzPAFHydrogenLoss:
    """A loss of H2 is written in Hill order (``-H2``), like every other formula delta."""

    def test_delta_h2(self):
        b3 = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"H2": 1})
        assert b3.to_mzpaf() == "b3{PEP}+H2"
        b3 = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"H2": -1})
        assert b3.to_mzpaf() == "b3{PEP}-H2"

    def test_internal_offsets(self):
        ax = pt.parse("PEPTIDE").frag(ion_type="ax", charge=1, position=(2, 5))
        bx = pt.parse("PEPTIDE").frag(ion_type="bx", charge=1, position=(2, 5))
        assert ax.to_mzpaf() == "m2:5{EPTI}-H2"
        assert bx.to_mzpaf() == "m2:5{EPTI}+CO-H2"


class TestImmoniumGlobalMods:
    """Global fixed mods and isotope labels on the immonium residue are part of its label."""

    @pytest.mark.parametrize(
        ("seq", "pos", "label"),
        [
            ("<[Oxidation]@P>PEPTIDE", 1, "IP[Oxidation]"),
            ("<[Oxidation]@E>PEPTIDE", 1, "IP"),
            ("<[Acetyl]@N-term>PEPTIDE", 1, "IP[Acetyl]"),
            ("<[Acetyl]@N-term>PEPTIDE", 3, "IP"),
            ("<[+15.995]@P>PEPTIDE", 1, "IP[+15.995]"),
            ("<13C>PEPTIDE", 1, "IP+4i13C"),
            ("<15N>PEPTIDE", 1, "IP+i15N"),
            ("<13C><15N>PEPTIDE", 1, "IP+4i13C+i15N"),
            ("<D>PEPTIDE", 1, "IP+7i2H"),
            ("<13C>[Acetyl]-PEPTIDE", 1, "IP[Acetyl]+6i13C"),
            ("<[Oxidation]@P><13C>PEPTIDE", 1, "IP[Oxidation]+4i13C"),
            ("P[U:Oxidation]EPTIDE", 1, "IP[Oxidation]"),
        ],
    )
    def test_label(self, seq, pos, label):
        assert pt.parse(seq).frag(ion_type="i", charge=1, position=pos).to_mzpaf() == label

    @pytest.mark.parametrize("seq", ["<[Oxidation]@P>P[Phospho]EPTIDE", "<[Oxidation]@P>[Acetyl]-PEPTIDE", "<13C>P[+15.995]EPTIDE"])
    def test_unwritable_raises(self, seq):
        frag = pt.parse(seq).frag(ion_type="i", charge=1, position=1)
        with pytest.raises(pt.PeptacularError):
            frag.to_mzpaf()


class TestLabelledDeprotonation:
    """Deprotonating a labelled ion removes the hydrogen isotope the ion holds."""

    @staticmethod
    def _mass(comp: Counter) -> float:
        return sum(e.get_mass(monoisotopic=True) * n for e, n in comp.items())

    @pytest.mark.parametrize("seq", ["<D>PEK", "<2H>PEK", "<13C>PEK", "PEK"])
    @pytest.mark.parametrize("charge", [-1, -2])
    def test_negative_charge_mass_matches_composition(self, seq, charge):
        from peptacular.constants import ELECTRON_MASS

        neutral = pt.parse(seq).frag("y", charge=0, calculate_with_composition=True).composition
        frag = pt.parse(seq).frag("y", charge=charge, calculate_with_composition=True)
        comp = frag.composition
        assert all(n > 0 for n in comp.values())
        h_neutral = sum(n for e, n in neutral.items() if e.symbol == "H")
        h_ion = sum(n for e, n in comp.items() if e.symbol == "H")
        assert h_ion == h_neutral + charge
        from peptacular.constants import HYDROGEN_BINDING_MASS

        expected = self._mass(comp) - charge * ELECTRON_MASS + charge * HYDROGEN_BINDING_MASS
        assert frag.mass == pytest.approx(expected, abs=1e-11)
        # The mass path (isotope as mass) agrees with the composition path.
        assert pt.parse(seq).frag("y", charge=charge).mass == pytest.approx(frag.mass, abs=1e-9)

    @pytest.mark.parametrize("seq", ["<D>PEK", "<2H>PEK"])
    def test_deuterated_ion_loses_a_deuteron(self, seq):
        comp = pt.parse(seq).frag("b", charge=-1, calculate_with_composition=True).composition
        assert {e.mass_number for e in comp if e.symbol == "H"} == {2}
        assert pt.fragment(seq, ion_types="b", charges=-1)

    def test_carbon13_ion_loses_a_protium(self):
        comp = pt.parse("<13C>PEK").frag("y", charge=-1, calculate_with_composition=True).composition
        assert {e.mass_number for e in comp if e.symbol == "H"} == {None}
        assert {e.mass_number for e in comp if e.symbol == "C"} == {13}


class TestChargedIsotopePeaks:
    """A charged envelope's monoisotopic peak sits on the ion's mass (PROTON_MASS per proton)."""

    @pytest.mark.parametrize("charge", [1, 2, 3, -1, -2, "Na:z+1", None])
    def test_monoisotopic_peak_matches_fragment(self, charge):
        annot = pt.parse("PEPTIDEK")
        peak = annot.isotopic_distribution(charge=charge)[0]
        exact = annot.frag(charge=charge, calculate_with_composition=True)
        assert peak.neutron_count == 0
        assert peak.mass == pytest.approx(exact.mass, abs=1e-11)
        assert peak.mass / max(abs(exact.charge_state), 1) == pytest.approx(annot.mz(charge=charge), abs=1e-9)

    def test_neutral_distribution_unchanged(self):
        annot = pt.parse("PEPTIDEK")
        neutral = annot.isotopic_distribution(charge=0)
        assert neutral == pt.brain_isotopic_distribution(annot.comp(), charge=0)


class TestImmoniumLabelOnFinalIon:
    """An immonium isotope label counts the labelled atoms of the final ion."""

    _SHIFT = {"2H": 2.01410177812 - 1.00782503223, "15N": 15.00010889888 - 14.00307400443, "13C": 1.00335483507}

    @pytest.mark.parametrize(
        ("labelled", "plain", "kwargs", "label", "iso", "n"),
        [
            ("<2H>P", "P", {"charge": -1}, "IP+6i2H^-1", "2H", 6),
            ("<D>P", "P", {"charge": -1}, "IP+6i2H^-1", "2H", 6),
            ("<2H>P", "P", {"charge": 1}, "IP+7i2H", "2H", 7),
            ("<15N>K", "K", {"charge": 1, "deltas": {pt.NeutralDelta.AMMONIA: 1}}, "IK-NH3+i15N", "15N", 1),
            ("<15N>K", "K", {"charge": 1, "deltas": {"NH3": -1}}, "IK+NH3+3i15N", "15N", 3),
            ("<13C>P", "P", {"charge": 2}, "IP+4i13C^2", "13C", 4),
            ("<15N>K", "K", {"charge": 1, "deltas": {17.0: -1}}, "IK-17.0+2i15N", "15N", 2),
        ],
    )
    def test_label_matches_mz(self, labelled, plain, kwargs, label, iso, n):
        frag = pt.parse(labelled).frag(ion_type="i", **kwargs)
        assert frag.to_mzpaf() == label
        unlabelled = pt.parse(plain).frag(ion_type="i", **kwargs)
        z = abs(frag.charge_state)
        assert frag.mz == pytest.approx(unlabelled.mz + n * self._SHIFT[iso] / z, abs=1e-9)


class TestImmoniumLabelWithCarriers:
    """A non-proton charge carrier is written as an adduct and is not counted in the label."""

    _SHIFT = TestImmoniumLabelOnFinalIon._SHIFT

    @pytest.mark.parametrize(
        ("labelled", "plain", "charge", "plain_charge", "label", "iso", "n"),
        [
            ("<2H>P", "P", "D:z+1", "D:z+1", "IP+7i2H[M+[2H]]", "2H", 7),
            # An unlabelled P has no 2H for D-1:z-1 to remove, so compare with plain
            # deprotonation (-H): 7 labelled atoms minus the deuteron the adduct takes = 6.
            ("<2H>P", "P", "D-1:z-1", -1, "IP+7i2H[M-[2H]]^-1", "2H", 6),
            ("<13C>P", "P", "[13C]:z+1", "[13C]:z+1", "IP+4i13C[M+[13C]]", "13C", 4),
            ("<2H>P", "P", "Na:z+1", "Na:z+1", "IP+7i2H[M+Na]", "2H", 7),
            ("<15N>K", "K", "NH4:z+1", "NH4:z+1", "IK+2i15N[M+NH4]", "15N", 2),
        ],
    )
    def test_label_reads_back_to_mz(self, labelled, plain, charge, plain_charge, label, iso, n):
        frag = pt.parse(labelled).frag(ion_type="i", charge=charge)
        assert frag.to_mzpaf() == label
        unlabelled = pt.parse(plain).frag(ion_type="i", charge=plain_charge)
        z = abs(frag.charge_state)
        # A plain integer charge is a proton carrier, which adds HYDROGEN_BINDING_MASS per proton.
        binding = -plain_charge * pt.HYDROGEN_BINDING_MASS / z if isinstance(plain_charge, int) else 0.0
        assert frag.mz == pytest.approx(unlabelled.mz + n * self._SHIFT[iso] / z + binding, abs=1e-10)
        try:
            import paftacular
        except ImportError:
            return
        assert paftacular.parse(label).mz() == pytest.approx(frag.mz, abs=1e-9)


class TestIsotopeCarrierNeedsIsotope:
    """A carrier that removes 2H needs a 2H in the ion; it never takes a light H instead."""

    def test_mass_raises(self):
        with pytest.raises(pt.PeptacularError):
            pt.parse("PEK").mass(charge="D-1:z-1")

    def test_fragment_raises(self):
        with pytest.raises(pt.PeptacularError):
            pt.parse("PEK").frag(ion_type="y", position=2, charge="D-1:z-1")


class TestCarrierOrder:
    """Charge carriers are summed before they touch the composition, so their order does not matter."""

    @pytest.mark.parametrize("carriers", [["H:z+1", "H-1:z-1", "H-1:z-1"], ["H-1:z-1", "H-1:z-1", "H:z+1"], ["H-1:z-1", "H:z+1", "H-1:z-1"]])
    def test_order_independent(self, carriers):
        frag = pt.parse("<D>PEK").frag(ion_type="p", charge=carriers)
        h = {(e.mass_number, n) for e, n in frag.composition.items() if e.symbol == "H"}
        assert h == {(2, 27)}
        assert frag.mass == pytest.approx(398.3631, abs=1e-4)

    @pytest.mark.parametrize(("ion_type", "position"), [("i", 1), ("b", 1)])
    def test_existing_deficit_still_raises(self, ion_type, position):
        annot = pt.parse("G[cysteinyl mercury][Label:2H(6)15N(1)]IVK")
        with pytest.raises(pt.InvalidAdjustmentError):
            annot.frag(ion_type=ion_type, position=position, charge=-1)


class TestSatelliteIons:
    """d and w ions need an unmodified cleaved residue and the atoms their offset removes."""

    @pytest.mark.parametrize(
        ("seq", "ion_type", "position"),
        [
            ("PEPV[Oxidation]K", "d", 4),
            ("<[Oxidation]@V>PEPVK", "d", 4),
            ("PEPTV[Oxidation]K", "w", 2),
            ("PET[Phospho]VK", "wa", 3),
            ("<[Oxidation]@T>PETVK", "wa", 3),
        ],
    )
    def test_modified_cleaved_residue_raises(self, seq, ion_type, position):
        with pytest.raises(pt.PeptacularError, match="carries a modification"):
            pt.parse(seq).frag(ion_type=ion_type, position=position, charge=1)

    def test_v_ion_drops_the_modification(self):
        assert pt.parse("PEPTV[Oxidation]K").frag(ion_type="v", position=2, charge=1).to_mzpaf() == "v2{V[Oxidation]K}"

    def test_series_skips_modified_residue(self):
        d = [f.to_mzpaf() for f in pt.fragment("PEPV[Oxidation]KV", ion_types="d", charges=1)]
        assert "d4{PEPV[Oxidation]}" not in d
        assert "d5{PEPV[Oxidation]K}" in d
        w = [f.to_mzpaf() for f in pt.fragment("PEPV[Oxidation]KV", ion_types="w", charges=1)]
        assert "w3{V[Oxidation]KV}" not in w
        assert "w2{KV}" in w

    def test_series_skips_impossible_composition(self):
        assert pt.fragment("V-[Amidated]", ion_types="d", charges=-1) == []
        with pytest.raises(pt.InvalidAdjustmentError):
            pt.parse("V-[Amidated]").frag(ion_type="d", position=1, charge=-1)

    def test_impossible_caller_isotopes_raise(self):
        with pytest.raises(pt.PeptacularError, match="do not fit any requested ion"):
            pt.fragment("PEPTIDE", ion_types="b", charges=1, isotopes=[{"13C": 100}])
        assert len(pt.fragment("PEPTIDE", ion_types="b", charges=1, isotopes=[{"13C": 1}])) == 7

    @pytest.mark.parametrize("composition", [False, True])
    def test_caller_isotopes_skip_ions_they_do_not_fit(self, composition):
        # b1, b2, y1, y2 have fewer than 3 N; the other ten ions of b and y take the swap.
        frags = pt.fragment("PEPTIDE", ion_types=["b", "y"], charges=1, isotopes=[{"15N": 3}], calculate_with_composition=composition)
        assert sorted((f.ion_type.value, f.position) for f in frags) == sorted((t, i) for t in "by" for i in range(3, 8))

    def test_caller_deltas_still_raise(self):
        with pytest.raises(pt.InvalidAdjustmentError):
            pt.fragment("PEPTIDE", ion_types="b", charges=1, deltas=[{"H3PO4": 1}])
        with pytest.raises(pt.InvalidAdjustmentError):
            pt.fragment("PEPTIDE", ion_types="b", charges=1, deltas=[{"H3PO4": 1}], isotopes=[{"15N": 3}])

    def test_any_series_skips_impossible_ion(self):
        assert pt.fragment("G-[Amidated]", ion_types="a", charges=-1) == []
        assert [f.to_mzpaf() for f in pt.fragment("GK-[Amidated]", ion_types="a", charges=-1)] == ["a1{G}^-1", "a2{GK-[Amidated]}^-1"]
        with pytest.raises(pt.InvalidAdjustmentError):
            pt.parse("G-[Amidated]").frag(ion_type="a", position=1, charge=-1)


class TestMapIsotopes:
    def test_map_isotopes(self):
        annot = pt.parse("<13C><15N>PEP")
        mapping = {str(k.symbol): v.mass_number for k, v in annot.map_isotopes().items()}
        assert mapping == {"C": 13, "N": 15}
        assert not hasattr(annot, "_map_isotopes")
        assert pt.parse("PEP").map_isotopes() == {}


class TestUnchargedMzPAF:
    """mzPAF has no uncharged ion: a label with no charge reads as +1."""

    @pytest.mark.parametrize("kwargs", [{"ion_type": "b", "position": 3}, {"ion_type": "p"}, {"ion_type": "b", "position": 3, "charge": 0}])
    def test_raises(self, kwargs):
        frag = pt.parse("PEPK").frag(**kwargs)
        with pytest.raises(pt.PeptacularError, match="uncharged"):
            frag.to_mzpaf()
        with pytest.raises(pt.PeptacularError, match="uncharged"):
            frag.serialize(format="mzpaf")

    def test_records_give_none(self):
        assert pt.fragment_records([pt.parse("PEPK").frag(ion_type="b", position=3)])[0]["mzpaf"] is None
