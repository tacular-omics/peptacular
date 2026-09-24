"""Localisation isomers, candidate sites and site-determining ions."""

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import peptacular as pt


def _ser(isomers):
    return [a.serialize() for a in isomers]


# --- localization_isomers: hand-checked cases -------------------------------------------


def test_phospho_range_over_s_and_t():
    assert _ser(pt.localization_isomers("PEP(ST)[Phospho]IDE")) == ["PEPS[Phospho]TIDE", "PEPST[Phospho]IDE"]


def test_unknown_position_goes_on_any_residue():
    assert _ser(pt.localization_isomers("[Phospho]?PEST")) == ["P[Phospho]EST", "PE[Phospho]ST", "PES[Phospho]T", "PEST[Phospho]"]


def test_unknown_position_with_count_uses_distinct_residues():
    assert _ser(pt.localization_isomers("[Phospho]^2?STY")) == ["S[Phospho]T[Phospho]Y", "S[Phospho]TY[Phospho]", "ST[Phospho]Y[Phospho]"]


def test_group_score_carried_onto_each_isomer():
    isomers = _ser(pt.localization_isomers("PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE"))
    assert isomers == ["PEPS[Phospho#g1(0.8)]TIDE", "PEPST[Phospho#g1(0.2)]IDE"]


def test_group_without_scores_keeps_label_only():
    assert _ser(pt.localization_isomers("PEPS[Phospho#g1]T[#g1]IY[#g1]DE")) == [
        "PEPS[Phospho#g1]TIYDE",
        "PEPST[Phospho#g1]IYDE",
        "PEPSTIY[Phospho#g1]DE",
    ]


def test_group_with_partial_scores():
    assert _ser(pt.localization_isomers("S[Phospho#g1]T[#g1(0.4)]")) == ["S[Phospho#g1]T", "ST[Phospho#g1(0.4)]"]


def test_group_member_keeps_its_other_mods():
    assert _ser(pt.localization_isomers("M[Oxidation][#g1(0.3)]S[Phospho#g1(0.7)]")) == [
        "M[Oxidation][Phospho#g1(0.3)]S",
        "M[Oxidation]S[Phospho#g1(0.7)]",
    ]


def test_group_inside_range_expands_over_the_group():
    assert _ser(pt.localization_isomers("PEP(S[Phospho#g1(0.3)]T[#g1(0.7)])IDE")) == ["PEPS[Phospho#g1(0.3)]TIDE", "PEPST[Phospho#g1(0.7)]IDE"]


def test_several_ambiguities_multiply_in_fixed_order():
    isomers = _ser(pt.localization_isomers("[Oxidation]?M(ST)[Phospho]M"))
    assert isomers == [
        "M[Oxidation]S[Phospho]TM",
        "MS[Phospho][Oxidation]TM",
        "MS[Phospho]T[Oxidation]M",
        "MS[Phospho]TM[Oxidation]",
        "M[Oxidation]ST[Phospho]M",
        "MS[Oxidation]T[Phospho]M",
        "MST[Phospho][Oxidation]M",
        "MST[Phospho]M[Oxidation]",
    ]


def test_identical_mods_are_deduplicated():
    # Swapping the two Phospho copies between the unknown mod and the range gives the same peptide.
    isomers = _ser(pt.localization_isomers("[Phospho]?(ST)[Phospho]"))
    assert isomers == ["S[Phospho][Phospho]T", "S[Phospho]T[Phospho]", "ST[Phospho][Phospho]"]
    assert len(isomers) == len(set(isomers))


def test_two_ranges_with_the_same_mod_dedup():
    isomers = _ser(pt.localization_isomers("(ST)[Phospho](ST)[Phospho]"))
    assert len(isomers) == 4
    assert len(set(isomers)) == 4


def test_no_ambiguity_returns_a_copy():
    annot = pt.parse("PEPS[Phospho]TIDE")
    isomers = pt.localization_isomers(annot)
    assert _ser(isomers) == ["PEPS[Phospho]TIDE"]
    assert isomers[0] is not annot


def test_input_is_not_modified():
    annot = pt.parse("[Phospho]?PEP(ST)[Oxidation]IDE")
    before = annot.serialize()
    annot.localization_isomers()
    assert annot.serialize() == before


def test_other_fields_are_kept():
    isomers = _ser(pt.localization_isomers("<13C>[Acetyl]-PEP(ST)[Phospho]IDE-[Amidated]/2"))
    assert isomers == ["<13C>[Acetyl]-PEPS[Phospho]TIDE-[Amidated]/2", "<13C>[Acetyl]-PEPST[Phospho]IDE-[Amidated]/2"]


def test_ambiguous_order_range_keeps_its_range():
    isomers = pt.localization_isomers("P(?ST)[Phospho]K")
    assert len(isomers) == 2
    assert all(len(a.intervals) == 1 and a.intervals[0].ambiguous and not a.intervals[0].has_mods for a in isomers)


def test_cross_link_labels_are_left_alone():
    assert _ser(pt.localization_isomers("K[XLMOD:02001#XL1]PEPK[#XL1]")) == ["K[XLMOD:02001#XL1]PEPK[#XL1]"]


def test_method_matches_function():
    annot = pt.parse("[Phospho]?PEPTIDE")
    assert _ser(annot.localization_isomers()) == _ser(pt.localization_isomers(annot))


def test_accepts_has_sequence_objects():
    class Entry:
        sequence = "PEP(ST)[Phospho]IDE"

    assert len(pt.localization_isomers(Entry())) == 2


# --- max_isomers ------------------------------------------------------------------------


def test_max_isomers_allows_exact_count():
    assert len(pt.localization_isomers("[Phospho]?PEST", max_isomers=4)) == 4


def test_max_isomers_raises_when_exceeded():
    with pytest.raises(pt.PeptacularError, match="max_isomers=3"):
        pt.localization_isomers("[Phospho]?PEST", max_isomers=3)


def test_max_isomers_stops_early_on_huge_expansions():
    with pytest.raises(pt.PeptacularError):
        pt.localization_isomers("[Phospho]^5?" + "S" * 60, max_isomers=10)


def test_max_isomers_counts_after_dedup():
    assert len(pt.localization_isomers("[Phospho]?(ST)[Phospho]", max_isomers=3)) == 3


@pytest.mark.parametrize("bad", [0, -1, 1.5, True, "3"])
def test_max_isomers_must_be_positive_int(bad):
    with pytest.raises(pt.PeptacularError):
        pt.localization_isomers("[Phospho]?PEST", max_isomers=bad)


def test_max_isomers_is_keyword_only():
    with pytest.raises(TypeError):
        pt.localization_isomers("[Phospho]?PEST", 3)  # type: ignore[misc]


# --- malformed or unsupported groups ------------------------------------------------------


@pytest.mark.parametrize(
    "sequence",
    [
        "[Phospho#g1]-S[#g1]T",
        "ST-[Phospho#g1]",
        "[Phospho#g1]?ST",
        "(ST)[Phospho#g1]",
        "{Phospho#g1}ST",
    ],
)
def test_group_off_a_residue_is_unsupported(sequence):
    with pytest.raises(pt.UnsupportedOperationError):
        pt.localization_isomers(sequence)


def test_group_with_only_references_raises():
    with pytest.raises(pt.PeptacularError, match="no modification"):
        pt.localization_isomers("S[#g1]T[#g1]")


def test_group_with_two_different_mods_raises():
    with pytest.raises(pt.PeptacularError, match="more than one"):
        pt.localization_isomers("S[Phospho#g1]T[Oxidation#g1]")


def test_more_copies_than_residues_raises():
    with pytest.raises(pt.PeptacularError, match="Cannot place 3"):
        pt.localization_isomers("[Phospho]^3?ST")


# --- round trips and invariants ---------------------------------------------------------


@pytest.mark.parametrize(
    "sequence",
    [
        "PEP(ST)[Phospho]IDE",
        "[Phospho]^2?PEPTIDE",
        "PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE",
        "[Oxidation]?M(ST)[Phospho]M",
        "PEPS[Phospho#g1]T[#g1]IY[#g1]DE",
    ],
)
def test_serialize_parse_round_trip(sequence):
    for isomer in pt.localization_isomers(sequence):
        text = isomer.serialize()
        assert pt.parse(text) == isomer
        assert pt.parse(text).serialize() == text
        # an isomer has nothing left to expand
        assert _ser(pt.localization_isomers(text)) == [text]


_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
_MODS = ["Phospho", "Oxidation", "Acetyl", "Methyl", "+15.995"]


@st.composite
def _ambiguous_peptides(draw):
    seq = draw(st.text(alphabet=_RESIDUES, min_size=2, max_size=12))
    kind = draw(st.sampled_from(["unknown", "range", "group", "mixed"]))
    mod = draw(st.sampled_from(_MODS))
    if kind in ("unknown", "mixed"):
        count = draw(st.integers(min_value=1, max_value=min(2, len(seq))))
        prefix = f"[{mod}]^{count}?" if count > 1 else f"[{mod}]?"
    else:
        prefix = ""
    if kind in ("range", "mixed"):
        start = draw(st.integers(min_value=0, max_value=len(seq) - 1))
        end = draw(st.integers(min_value=start + 1, max_value=len(seq)))
        body = f"{seq[:start]}({seq[start:end]})[{mod}]{seq[end:]}"
    elif kind == "group":
        positions = sorted(draw(st.sets(st.integers(min_value=0, max_value=len(seq) - 1), min_size=1, max_size=min(4, len(seq)))))
        scores = draw(st.booleans())
        parts = []
        for i, residue in enumerate(seq):
            if i == positions[0]:
                parts.append(f"{residue}[{mod}#g1{'(0.5)' if scores else ''}]")
            elif i in positions:
                parts.append(f"{residue}[#g1{'(0.25)' if scores else ''}]")
            else:
                parts.append(residue)
        body = "".join(parts)
    else:
        body = seq
    return prefix + body


@settings(max_examples=150, deadline=None)
@given(_ambiguous_peptides())
def test_property_isomers_keep_composition_and_mass(sequence):
    annot = pt.parse(sequence)
    isomers = pt.localization_isomers(annot)
    assert isomers
    texts = _ser(isomers)
    assert len(texts) == len(set(texts))
    mass = annot.mass()
    comp = None if "+15.995" in sequence else annot.comp()  # a mass-only mod has no composition
    for isomer in isomers:
        assert not isomer.has_unknown_mods
        assert all(not interval.has_mods for interval in isomer.intervals)
        assert isomer.mass() == pytest.approx(mass, abs=1e-6)
        if comp is not None:
            assert isomer.comp() == comp
        assert pt.parse(isomer.serialize()) == isomer
    # deterministic
    assert _ser(pt.localization_isomers(sequence)) == texts


# --- candidate_sites ----------------------------------------------------------------------


def test_candidate_sites_phospho():
    sites = pt.candidate_sites("PEPSTIDEYK", "Phospho", residues="STY")
    assert [(i, a.serialize()) for i, a in sites] == [
        (3, "PEPS[Phospho]TIDEYK"),
        (4, "PEPST[Phospho]IDEYK"),
        (8, "PEPSTIDEY[Phospho]K"),
    ]


def test_candidate_sites_skips_modified_residues():
    sites = pt.candidate_sites("PEPS[Phospho]TIDE", "Phospho", residues="ST")
    assert [i for i, _ in sites] == [4]


def test_candidate_sites_lowercase_residues_and_mass_mod():
    sites = pt.candidate_sites("MSM", 15.995, residues="m")
    assert [(i, a.serialize()) for i, a in sites] == [(0, "M[+15.995]SM"), (2, "MSM[+15.995]")]


def test_candidate_sites_no_match_is_empty():
    assert pt.candidate_sites("PEPTIDE", "Phospho", residues="Y") == []


def test_candidate_sites_method_and_input_untouched():
    annot = pt.parse("PEPST")
    sites = annot.candidate_sites("Phospho", residues="ST")
    assert [i for i, _ in sites] == [3, 4]
    assert annot.serialize() == "PEPST"


@pytest.mark.parametrize("bad", ["", "S1", None, 5])
def test_candidate_sites_rejects_bad_residues(bad):
    with pytest.raises(pt.PeptacularError):
        pt.candidate_sites("PEPST", "Phospho", residues=bad)


def test_candidate_sites_residues_is_required():
    with pytest.raises(TypeError):
        pt.candidate_sites("PEPST", "Phospho")  # type: ignore[call-arg]


# --- site_determining_ions ----------------------------------------------------------------


def _labels(per_isomer):
    return [[f"{f.ion_type}{f.position}+{f.charge_state}" for f in frags] for frags in per_isomer]


def test_site_determining_ions_for_adjacent_sites():
    isomers = pt.localization_isomers("PEP(ST)[Phospho]IDE")
    assert _labels(pt.site_determining_ions(isomers)) == [["b4+1", "y4+1"], ["b4+1", "y4+1"]]


def test_site_determining_ions_for_distant_sites():
    isomers = pt.localization_isomers("PES[Phospho#g1]AAT[#g1]K")
    ions = pt.site_determining_ions(isomers, ion_types=("b",))
    assert _labels(ions) == [["b3+1", "b4+1", "b5+1"], ["b3+1", "b4+1", "b5+1"]]


def test_site_determining_ions_are_fragment_objects_from_fragment():
    isomers = pt.localization_isomers("PEP(ST)[Phospho]IDE")
    ions = pt.site_determining_ions(isomers, charges=(1, 2))
    for isomer, frags in zip(isomers, ions, strict=True):
        all_frags = isomer.fragment(ion_types=("b", "y"), charges=(1, 2))
        keyed = {(f.ion_type, f.position, f.charge_state): f.mz for f in all_frags}
        for f in frags:
            assert isinstance(f, pt.Fragment)
            assert keyed[(f.ion_type, f.position, f.charge_state)] == f.mz
    assert {f.charge_state for frags in ions for f in frags} == {1, 2}


def test_site_determining_ions_checks_against_every_ion_of_other_isomers():
    # Isomer A's b-ion may coincide with a y-ion of isomer B; it is then not site-determining.
    isomers = pt.localization_isomers("PEP(ST)[Phospho]IDE")
    ions = pt.site_determining_ions(isomers)
    other_mz = [f.mz for f in isomers[1].fragment(ion_types=("b", "y"), charges=(1,))]
    for f in ions[0]:
        assert all(abs(f.mz - mz) > 1e-6 for mz in other_mz)


def test_site_determining_ions_tolerance_da_and_ppm():
    isomers = pt.localization_isomers("PEP(ST)[Phospho]IDE")
    exact = _labels(pt.site_determining_ions(isomers))
    assert _labels(pt.site_determining_ions(isomers, tolerance=0.02)) == exact
    assert _labels(pt.site_determining_ions(isomers, tolerance=10, unit="ppm")) == exact
    # A tolerance wider than the whole spectrum leaves nothing site-determining.
    assert pt.site_determining_ions(isomers, tolerance=5000.0) == [[], []]
    assert pt.site_determining_ions(isomers, tolerance=1e7, unit="ppm") == [[], []]


def test_site_determining_ions_single_isomer_returns_everything():
    (isomer,) = pt.localization_isomers("PEPTIDE")
    (frags,) = pt.site_determining_ions([isomer])
    assert len(frags) == len(isomer.fragment(ion_types=("b", "y"), charges=(1,)))


def test_site_determining_ions_accepts_strings():
    assert _labels(pt.site_determining_ions(["PEPS[Phospho]TIDE", "PEPST[Phospho]IDE"])) == [["b4+1", "y4+1"], ["b4+1", "y4+1"]]


def test_site_determining_ions_empty_input():
    assert pt.site_determining_ions([]) == []


@pytest.mark.parametrize("kwargs", [{"unit": "mz"}, {"tolerance": -1.0}, {"tolerance": True}, {"tolerance": float("nan")}])
def test_site_determining_ions_rejects_bad_options(kwargs):
    with pytest.raises(pt.PeptacularError):
        pt.site_determining_ions(["PEPTIDE"], **kwargs)
