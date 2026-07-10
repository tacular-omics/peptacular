"""Tests for cross-linked / multi-chain peptidoform ions (ProForma 2.1 §9.2.2, §9.3).

Covers :class:`peptacular.MultiProFormaAnnotation` and the cross-link label validation
that backs it, plus the routing in ``ProFormaAnnotation.parse`` / ``pt.parse`` that returns
a container for ``//``-joined chains.
"""

import math

import pytest

import peptacular as pt
from peptacular import MultiProFormaAnnotation
from peptacular.annotation.annotation import validate_crosslink_labels

# A well-formed inter-chain cross-link: the DSBU linker (XLMOD:02001, +138.068 Da) defined
# on chain A and referenced from chain B.
INTERCHAIN = "EVTK[XLMOD:02001#XL1]LE//AK[#XL1]ENLYFQ/3"
# A branched peptide (ProForma 2.1 §9.3).
BRANCH = "ED[MOD:00093#BRANCH]//D[#BRANCH]ATR/1"

LINKER_MASS = 138.06808  # XLMOD:02001 delta


class TestParseRouting:
    def test_parse_returns_multi_for_crosslink(self):
        ion = pt.parse(INTERCHAIN)
        assert isinstance(ion, MultiProFormaAnnotation)
        assert len(ion) == 2

    def test_parse_returns_annotation_for_single_chain(self):
        annot = pt.parse("PEPTIDE")
        assert isinstance(annot, pt.ProFormaAnnotation)

    def test_intrachain_crosslink_stays_single_annotation(self):
        # Both ends on one chain -> still a single ProFormaAnnotation.
        annot = pt.parse("EVTK[XLMOD:02001#XL1]LEK[#XL1]SEFD")
        assert isinstance(annot, pt.ProFormaAnnotation)

    def test_chimeric_rejected_with_helpful_message(self):
        with pytest.raises(ValueError, match="parse_chimeric"):
            pt.parse("PEPTIDE+SEQUENCE")

    def test_combined_chimeric_and_crosslink_rejected(self):
        with pytest.raises(ValueError, match="not supported"):
            MultiProFormaAnnotation.parse("PEK[XLMOD:02001#XL1]//K[#XL1]AR+SEQUENCE")


class TestRoundTrip:
    @pytest.mark.parametrize(
        "seq",
        [
            INTERCHAIN,
            BRANCH,
            "EVTK[XLMOD:02001#XL1]LE//AK[#XL1]ENLYFQ",  # no charge
            "<13C>EVTK[XLMOD:02001#XL1]LE//AK[#XL1]ENLYFQ/2",  # global isotope mod
            "PEPK[XLMOD:02001#XL1]TIDE//SEQK[#XL1]//AAAK[#XL1]R/2",  # three chains
        ],
    )
    def test_round_trip(self, seq):
        assert pt.parse(seq).serialize() == seq

    def test_str_matches_serialize(self):
        ion = pt.parse(INTERCHAIN)
        assert str(ion) == ion.serialize()

    def test_global_mod_not_duplicated_across_chains(self):
        ion = pt.parse("<13C>PEPK[XLMOD:02001#XL1]//SEQK[#XL1]/2")
        # The <13C> prefix appears exactly once, not once per chain.
        assert ion.serialize().count("<13C>") == 1


class TestMass:
    def test_interchain_neutral_mass_sums_chains_plus_one_linker(self):
        ion = pt.parse(INTERCHAIN)
        backbone_a = pt.parse("EVTKLE").mass(charge=0)
        backbone_b = pt.parse("AKENLYFQ").mass(charge=0)
        expected = backbone_a + backbone_b + LINKER_MASS
        assert math.isclose(ion.neutral_mass(), expected, abs_tol=1e-3)

    def test_charged_mass_adds_protons_once(self):
        ion = pt.parse(INTERCHAIN)  # charge 3
        diff = ion.mass() - ion.neutral_mass()
        proton = 1.007276
        assert math.isclose(diff, 3 * proton, abs_tol=1e-3)

    def test_mz_is_mass_over_charge(self):
        ion = pt.parse(INTERCHAIN)
        assert math.isclose(ion.mz(), ion.mass() / 3, rel_tol=1e-9)

    def test_mz_neutral_raises(self):
        ion = MultiProFormaAnnotation.parse("EVTK[XLMOD:02001#XL1]LE//AK[#XL1]ENLYFQ")
        with pytest.raises(ValueError, match="neutral"):
            ion.mz()

    def test_average_mass_differs_from_monoisotopic(self):
        ion = pt.parse(INTERCHAIN)
        assert ion.mass(monoisotopic=False) > ion.mass(monoisotopic=True)

    def test_charge_override(self):
        ion = pt.parse(INTERCHAIN)  # native charge 3
        assert math.isclose(ion.mz(charge=1) * 1, ion.mass(charge=1), rel_tol=1e-9)
        assert not math.isclose(ion.mass(charge=1), ion.mass(charge=3), rel_tol=1e-6)

    def test_comp_mass_consistency(self):
        ion = pt.parse(INTERCHAIN)
        comp = ion.comp(charge=0)
        comp_mass = sum(el.get_mass() * n for el, n in comp.items())
        assert math.isclose(comp_mass, ion.neutral_mass(), abs_tol=1e-2)

    def test_adduct_charge_carrier_shared_across_chains(self):
        # A mixed adduct charge (list form) is applied once for the whole ion.
        seq = "EVTK[XLMOD:02001#XL1]LE//AK[#XL1]ENLYFQ/[Na:z+1,H:z+1]"
        ion = pt.parse(seq)
        assert ion.serialize() == seq
        assert ion.charge_state == 2
        assert math.isclose(ion.mz() * ion.charge_state, ion.mass(), rel_tol=1e-9)


class TestValidation:
    def test_valid_interchain_passes(self):
        # Should not raise.
        pt.parse(INTERCHAIN, validate=True)
        pt.parse(BRANCH, validate=True)

    def test_dangling_secondary_raises(self):
        # '#XL1' referenced but never defined.
        with pytest.raises(ValueError, match="referenced but never defined"):
            pt.parse("PEK[#XL1]TIDE//AKENLYFQ", validate=True)

    def test_dangling_primary_raises(self):
        # '#XL1' defined but never referenced.
        with pytest.raises(ValueError, match="never referenced"):
            pt.parse("PEK[XLMOD:02001#XL1]TIDE//AKENLYFQ", validate=True)

    def test_duplicate_primary_raises(self):
        with pytest.raises(ValueError, match="defined 2 times"):
            MultiProFormaAnnotation.parse("K[XLMOD:02001#XL1]A//K[XLMOD:02001#XL1]A//K[#XL1]A", validate=True)

    def test_validate_crosslink_labels_returns_empty_when_valid(self):
        chains = list(pt.parse(INTERCHAIN))
        assert validate_crosslink_labels(chains) == []

    def test_ambiguity_hash_not_treated_as_crosslink(self):
        # '#1' position grouping (§7.6) must not be mistaken for a cross-link.
        annot = pt.parse("PEP[Oxidation#1]M[#1]AT", validate=True)
        assert isinstance(annot, pt.ProFormaAnnotation)

    def test_single_chain_validate_crosslinks_method(self):
        # Intrachain: primary + secondary both present -> valid.
        pt.parse("EVTK[XLMOD:02001#XL1]LEK[#XL1]SEFD").validate_crosslinks()
        # Intrachain dangling secondary -> raises via the single-chain method.
        with pytest.raises(ValueError, match="referenced but never defined"):
            pt.parse("EVTK[#XL1]LESEFD").validate_crosslinks()

    def test_multi_validate_crosslinks_method(self):
        ion = MultiProFormaAnnotation.parse("PEK[#XL1]TIDE//AKENLYFQ")
        with pytest.raises(ValueError):
            ion.validate_crosslinks()


class TestFragmentation:
    def test_fragment_not_supported(self):
        ion = pt.parse(INTERCHAIN)
        with pytest.raises(NotImplementedError, match="not yet supported"):
            ion.fragment()


class TestContainerBehavior:
    def test_len(self):
        assert len(pt.parse(INTERCHAIN)) == 2

    def test_iter_and_index(self):
        ion = pt.parse(INTERCHAIN)
        chains = list(ion)
        assert len(chains) == 2
        assert ion[0] is chains[0]
        assert all(isinstance(c, pt.ProFormaAnnotation) for c in ion)

    def test_chains_are_charge_free(self):
        # The shared charge lives on the ion, not on the individual chains.
        ion = pt.parse(INTERCHAIN)
        assert all(chain.charge is None for chain in ion)
        assert ion.charge == 3
        assert ion.charge_state == 3

    def test_copy_is_independent_and_equal(self):
        ion = pt.parse(INTERCHAIN)
        clone = ion.copy()
        assert clone == ion
        assert clone is not ion
        assert clone.serialize() == ion.serialize()

    def test_equality_and_hash(self):
        a = pt.parse(INTERCHAIN)
        b = pt.parse(INTERCHAIN)
        assert a == b
        assert hash(a) == hash(b)

    def test_inequality_with_other_type(self):
        assert (pt.parse(INTERCHAIN) == "not an ion") is False

    def test_repr_mentions_chain_count(self):
        assert "chains=2" in repr(pt.parse(INTERCHAIN))

    def test_empty_chains_rejected(self):
        with pytest.raises(ValueError, match="at least one chain"):
            MultiProFormaAnnotation([])
