"""Fragment immutability, the prefix-sum fragment() fast path and fast_fragment() parity (5.0)."""

import copy
import pickle
from dataclasses import FrozenInstanceError

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from tacular import NEUTRAL_DELTA_LOOKUP, IonType

import peptacular as pt
from peptacular.annotation.cached_comps import DeltaInfo, IsotopeInfo
from peptacular.annotation.frag import Fragment

PEPTIDES = [
    "PEPTIDE",
    "[Acetyl]-PEPM[Oxidation]TIDEKC[Carbamidomethyl]LLSGR/2",
    "AS[Phospho]DFGHIKLMNPQR/3",
    "LLLSEEPQRVVK-[Amidated]/2",
    "PEPT[+79.966]IDE",
    "C[Carbamidomethyl]DEFGHK/[Na:z+1,H:z+1]",
    "{Glycan:Hex}PEPTIDE/2",
    "M",
    "DTVWAK",
]


def _fields(frag: Fragment) -> tuple:
    return (
        frag.ion_type,
        frag.position,
        frag.monoisotopic,
        frag.charge_state,
        frag._charge_adducts,
        frag.external_charge,
        frag._isotopes,
        frag._deltas,
        frag.parent_sequence,
        frag.parent_sequence_length,
    )


def _series_both(
    sequence: str,
    ion_type: IonType,
    charge: int,
    *,
    isotopes=(0,),
    deltas=(None,),
    neutral_deltas=(),
    monoisotopic: bool = True,
    max_deltas: int = 1,
    min_length: int | None = None,
    max_length: int | None = None,
) -> tuple[list[Fragment], list[Fragment]]:
    annot = pt.parse(sequence).set_charge(charge, inplace=False)
    kwargs = {
        "forward": pt.IonType(ion_type) in (IonType.A, IonType.B, IonType.C, IonType.D, IonType.DA, IonType.DB),
        "monoisotopic": monoisotopic,
        "isotopes": [IsotopeInfo.from_input(i) for i in isotopes],
        "deltas": [DeltaInfo.from_input(d) for d in deltas],
        "neutral_deltas": [NEUTRAL_DELTA_LOOKUP[n] for n in neutral_deltas],
        "calculate_with_composition": False,
        "parent_sequence": annot.serialize(),
        "parent_sequence_length": len(annot),
        "max_deltas": max_deltas,
        "min_length": min_length,
        "max_length": max_length,
    }
    fast = list(annot._fragment_series(ion_type, **kwargs))
    slow = list(annot._fragment_series(ion_type, _fast=False, **kwargs))
    return fast, slow


def _assert_same(fast: list[Fragment], slow: list[Fragment]) -> None:
    assert len(fast) == len(slow)
    for f, s in zip(fast, slow, strict=True):
        assert _fields(f) == _fields(s)
        assert f.mass == pytest.approx(s.mass, abs=1e-9)


class TestFragmentFrozen:
    def test_setattr_raises(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2)
        with pytest.raises(FrozenInstanceError):
            frag.mass = 1.0
        with pytest.raises(FrozenInstanceError):
            frag.new_attribute = 1

    def test_delattr_raises(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2)
        with pytest.raises(FrozenInstanceError):
            del frag.mass

    def test_slots_no_dict(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="y", charge=2, position=3)
        assert not hasattr(frag, "__dict__")

    def test_pickle_and_copy_round_trip(self):
        frag = pt.parse("PEM[Oxidation]TIDE").frag(ion_type="y", charge=2, position=3, deltas={"H2O": -1}, isotopes=1)
        for clone in (pickle.loads(pickle.dumps(frag)), copy.copy(frag), copy.deepcopy(frag)):
            assert _fields(clone) == _fields(frag)
            assert clone.mass == frag.mass
            assert clone.to_mzpaf() == frag.to_mzpaf()
            with pytest.raises(FrozenInstanceError):
                clone.position = 1

    def test_replace_returns_new_object(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2)
        other = frag._replace(mass=10.0)
        assert other is not frag
        assert other.mass == 10.0
        assert frag.mass != 10.0
        assert other.position == frag.position

    def test_value_equality_and_hash(self):
        a = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2, deltas={"H2O": -1})
        b = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2, deltas={"H2O": -1})
        assert a == b and hash(a) == hash(b)
        assert len({a, b}) == 1
        assert a != a.replace(mass=a.mass + 1)
        assert a != pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2)
        assert a != "b2"
        # the composition cache is not part of the value
        assert a == a.replace(composition=a.composition)

    def test_public_replace(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"H2O": -1})
        other = frag.replace(mass=10.0, deltas={"H-3N-1": 1})
        assert other.mass == 10.0 and frag.mass != 10.0
        assert other._deltas == {"H-3N-1": 1} and other.to_mzpaf(include_sequence=False) == "b3-NH3"
        with pytest.raises(TypeError, match="losses"):
            frag.replace(losses={})
        with pytest.raises(TypeError, match="_deltas"):
            frag.replace(_deltas={})

    def test_replace_charge_moves_external_charge(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3)
        two = frag.replace(charge_state=2)
        assert two.external_charge == 2 and two.charge_adducts.serialize() == frag.replace(charge_state=2, external_charge=None).charge_adducts.serialize()
        assert frag.replace(charge_state=2, external_charge=1).external_charge == 1

    def test_replace_drops_stale_composition(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, calculate_with_composition=True)
        assert frag._composition is not None
        assert frag.replace(mass=1.0)._composition is frag._composition
        assert frag.replace(position=2)._composition is None
        assert frag.replace(position=2).composition == pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2).composition

    def test_constructor_options_are_keyword_only(self):
        with pytest.raises(TypeError):
            pt.Fragment(pt.IonType.B, 1, 100.0, True, 1, None)  # type: ignore[misc]

    def test_composition_path_still_sets_deltas(self):
        # The composition branch of _frag_impl rebuilds the fragment instead of mutating it.
        mass_mode = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"H2O": -1})
        assert mass_mode._deltas
        assert mass_mode._composition is None
        comp_mode = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"H2O": -1}, calculate_with_composition=True)
        assert comp_mode._deltas == mass_mode._deltas
        assert comp_mode.mass == pytest.approx(mass_mode.mass, abs=1e-9)


class TestFragmentFastPath:
    @pytest.mark.parametrize("sequence", PEPTIDES)
    @pytest.mark.parametrize("ion_type", [IonType.A, IonType.B, IonType.C, IonType.X, IonType.Y, IonType.Z, IonType.Z_RADICAL])
    @pytest.mark.parametrize("charge", [1, 2, -1])
    def test_matches_slicing(self, sequence, ion_type, charge):
        _assert_same(*_series_both(sequence, ion_type, charge))

    @pytest.mark.parametrize("sequence", PEPTIDES)
    def test_average_masses(self, sequence):
        _assert_same(*_series_both(sequence, IonType.B, 2, monoisotopic=False))
        _assert_same(*_series_both(sequence, IonType.Y, 1, monoisotopic=False))

    @pytest.mark.parametrize("sequence", PEPTIDES)
    def test_isotopes_deltas_and_neutral_deltas(self, sequence):
        # Mixed: mass-only items take the fast path, formula deltas and isotope swaps slice.
        fast, slow = _series_both(
            sequence,
            IonType.B,
            2,
            isotopes=(0, 1, 2, {"15N": 1}),
            deltas=(None, 1.5, {"H2O": -1}, {-17.0: 2}),
            neutral_deltas=("H2O", "NH3"),
            max_deltas=2,
        )
        _assert_same(fast, slow)
        _assert_same(*_series_both(sequence, IonType.Y, 1, isotopes=(0, 1), deltas=(None, -18.0), neutral_deltas=("NH3",)))

    @pytest.mark.parametrize("ion_type", [IonType.D, IonType.DA, IonType.DB, IonType.W, IonType.WA, IonType.WB, IonType.V])
    def test_satellite_ions_fall_back(self, ion_type):
        for sequence in ("VTIDELK/2", "[Acetyl]-TIVM[Oxidation]DEIK"):
            _assert_same(*_series_both(sequence, ion_type, 1))

    def test_length_filters(self):
        _assert_same(*_series_both("PEPTIDEKLLSGR/2", IonType.B, 1, min_length=3, max_length=6))
        fast, _ = _series_both("PEPTIDEKLLSGR/2", IonType.Y, 2, min_length=3, max_length=6)
        assert [f.position for f in fast] == [3, 4, 5, 6]

    @pytest.mark.parametrize(
        "sequence",
        [
            "<13C>PEPTIDE",  # global isotope label
            "<[Carbamidomethyl]@C>PEPCTIDE",  # static mod
            "PEPT[Formula:Na:z+1]IDE",  # charged mod
            "[Formula:Na:z+1]-PEPTIDE",
            "PEPTIDE-[Formula:Na:z+1]",
        ],
    )
    def test_unsupported_annotations_slice(self, sequence):
        annot = pt.parse(sequence)
        assert annot._series_mass_vector(True, False) is None
        _assert_same(*_series_both(sequence, IonType.B, 1))
        _assert_same(*_series_both(sequence, IonType.Y, 2))

    def test_composition_mode_and_atom_removing_carrier_slice(self):
        annot = pt.parse("PEPTIDE")
        assert annot._series_mass_vector(True, True) is None
        assert pt.parse("PEPTIDE/[H-1:z-1]")._series_mass_vector(True, False) is None
        _assert_same(*_series_both("PEPTIDE/[H-1:z-1]", IonType.B, -1))

    def test_unknown_residue_mass_slices(self):
        # B has no defined mass, so the vector is unavailable and the fallback raises as before.
        assert pt.parse("PEPBIDE")._series_mass_vector(True, False) is None

    def test_public_fragment_uses_fast_path_values(self):
        sequence = "[Acetyl]-PEPM[Oxidation]TIDEKC[Carbamidomethyl]LLSGR/2"
        frags = pt.fragment(sequence, ion_types=("b", "y"), charges=(1, 2))
        annot = pt.parse(sequence)
        for frag in frags:
            expected = annot.frag(ion_type=frag.ion_type, charge=frag.charge_state, position=frag.position)
            assert frag.mass == pytest.approx(expected.mass, abs=1e-9)
            assert frag.to_mzpaf() == expected.to_mzpaf()

    @settings(max_examples=60, deadline=None)
    @given(
        residues=st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=12),
        mod_site=st.integers(min_value=0, max_value=11),
        mod=st.sampled_from(["Oxidation", "Phospho", "+15.995", "Formula:C2H2O", "Carbamidomethyl"]),
        nterm=st.booleans(),
        cterm=st.booleans(),
        charge=st.integers(min_value=1, max_value=3),
    )
    def test_property_matches_slicing(self, residues, mod_site, mod, nterm, cterm, charge):
        site = mod_site % len(residues)
        sequence = residues[:site] + residues[site] + f"[{mod}]" + residues[site + 1 :]
        if nterm:
            sequence = "[Acetyl]-" + sequence
        if cterm:
            sequence = sequence + "-[Amidated]"
        for ion_type in (IonType.B, IonType.Y, IonType.C, IonType.Z):
            _assert_same(*_series_both(sequence, ion_type, charge, deltas=(None, {"H2O": -1}, 2.0)))


class TestFastFragmentParity:
    @pytest.mark.parametrize("sequence", [p for p in PEPTIDES if "Na:" not in p and "Glycan" not in p])
    @pytest.mark.parametrize("ion_types", [("b", "y"), ("a", "c", "x", "z")])
    def test_same_ion_set_and_mz_as_fragment(self, sequence, ion_types):
        annot = pt.parse(sequence)
        charges = (1, 2)
        fast = annot.fast_fragment(ion_types=ion_types, charges=charges)
        full = {(f.ion_type, f.charge_state, f.position): f.mz for f in annot.fragment(ion_types=ion_types, charges=charges)}
        flat = {(ion, z, i + 1): mz for (ion, z), values in fast.items() for i, mz in enumerate(values)}
        assert set(flat) == set(full)
        for key, mz in full.items():
            assert flat[key] == pytest.approx(mz, abs=1e-9)

    def test_includes_full_length_ions(self):
        annot = pt.parse("PEPTIDE")
        fast = annot.fast_fragment(ion_types=("b", "y"), charges=(1,))
        assert len(fast[(IonType.B, 1)]) == len(annot)
        assert fast[(IonType.B, 1)][-1] == pytest.approx(annot.frag(ion_type="b", charge=1, position=7).mz, abs=1e-9)
        assert fast[(IonType.Y, 1)][-1] == pytest.approx(annot.frag(ion_type="y", charge=1, position=7).mz, abs=1e-9)

    def test_proton_is_hydrogen_minus_electron(self):
        # b1 of G at z=1 and z=2 differ by exactly one (H - e) carrier, divided out.
        fast = pt.parse("GG").fast_fragment(ion_types=("b",), charges=(1, 2))
        b1_z1 = fast[(IonType.B, 1)][0]
        b1_z2 = fast[(IonType.B, 2)][0]
        carrier = pt.parse("H").frag(ion_type="p", charge=1).mass - pt.parse("H").frag(ion_type="p", charge=0).mass
        assert 2 * b1_z2 - b1_z1 == pytest.approx(carrier, abs=1e-9)
        assert carrier == pytest.approx(pt.PROTON_MASS, abs=1e-6)

    @pytest.mark.parametrize("sequence", ["PEPTIDE", "GG", "[Acetyl]-S[Phospho]AMPLEK/2"])
    def test_b1_agrees_with_fragment(self, sequence):
        # Both paths add the same (H - e) carrier, so b1 and y1 agree far below 1e-9 Da.
        annot = pt.parse(sequence)
        fast = annot.fast_fragment(ion_types=("b", "y"), charges=(1, 2, 3))
        for frag in annot.fragment(ion_types=("b", "y"), charges=(1, 2, 3)):
            if frag.position == 1:
                assert abs(fast[(frag.ion_type, frag.charge_state)][0] - frag.mz) < 1e-12


class TestDeltaInfoAdd:
    def test_adding_empty_returns_the_other_operand(self):
        empty = DeltaInfo.from_input(None)
        water = DeltaInfo.from_input({"H2O": -1})
        assert water + empty is water
        assert empty + water is water
        assert empty + empty is empty

    def test_adding_non_empty_combines_counts(self):
        water = DeltaInfo.from_input({"H2O": -1})
        combined = water + DeltaInfo.from_input({"H2O": -1, 1.5: 1})
        assert sorted(combined.to_fragment_mapping.values()) == [-2, 1]
