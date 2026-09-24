"""Property-based tests (Hypothesis) for ProForma round-trips, fragmentation and digestion.

The default profile (``ci``) keeps these fast. Run a longer sweep with::

    HYPOTHESIS_PROFILE=thorough uv run pytest tests/test_hypothesis_properties.py
"""

from __future__ import annotations

import os

import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

import peptacular as pt

settings.register_profile("ci", max_examples=100, deadline=None, suppress_health_check=[HealthCheck.too_slow])
settings.register_profile("thorough", max_examples=2000, deadline=None, suppress_health_check=[HealthCheck.too_slow])
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "ci"))

PROTON = 1.007276466621
TOL = 1e-6

AAS = "ACDEFGHIKLMNPQRSTVWY"
# Modification tags that parse and carry a mass, covering every tag form: Unimod name and
# accession, PSI-MOD accession, mass deltas, formula and glycan.
RESIDUE_MODS = ["Oxidation", "Phospho", "U:Deamidated", "UNIMOD:35", "MOD:00046", "+15.995", "-17.027", "Formula:C2H2O", "Glycan:Hex1"]
NTERM_MODS = ["Acetyl", "UNIMOD:1", "+42.011", "Formula:C2H2O"]
CTERM_MODS = ["Amidated", "-0.984", "Methyl"]
LABILE_MODS = ["Glycan:Hex1HexNAc1", "Phospho", "Formula:HPO3"]
STATIC_MODS = ["<[Carbamidomethyl]@C>", "<[Oxidation]@M>", "<[+1.0]@K,R>"]
ISOTOPE_MODS = ["<13C>", "<15N>", "<D>"]

residue_seq = st.text(alphabet=AAS, min_size=1, max_size=15)


@st.composite
def proforma(draw, *, fragmentable: bool = False) -> str:
    """Build a valid ProForma 2.0 string.

    ``fragmentable`` restricts to features every fragmentation path supports: no labile,
    unknown-position or interval/ambiguity mods and no charge adducts.
    """
    seq = draw(residue_seq)
    residues = []
    ambiguous_done = False
    i = 0
    while i < len(seq):
        # Ambiguity group or interval over the next few residues.
        if not fragmentable and not ambiguous_done and i + 2 <= len(seq) and draw(st.booleans()):
            width = draw(st.integers(2, min(4, len(seq) - i)))
            chunk = seq[i : i + width]
            kind = draw(st.sampled_from(["ambiguous", "interval"]))
            if kind == "ambiguous":
                residues.append(f"(?{chunk})")
            else:
                residues.append(f"({chunk})[{draw(st.sampled_from(RESIDUE_MODS))}]")
            ambiguous_done = True
            i += width
            continue
        aa = seq[i]
        if draw(st.integers(0, 4)) == 0:
            aa += f"[{draw(st.sampled_from(RESIDUE_MODS))}]"
        residues.append(aa)
        i += 1
    body = "".join(residues)

    prefix = ""
    if draw(st.booleans()):
        prefix += "".join(draw(st.lists(st.sampled_from(ISOTOPE_MODS), max_size=1, unique=True)))
    if draw(st.booleans()):
        prefix += "".join(draw(st.lists(st.sampled_from(STATIC_MODS), max_size=2, unique=True)))
    if not fragmentable and draw(st.booleans()):
        prefix += f"[{draw(st.sampled_from(RESIDUE_MODS))}]?"
    if not fragmentable and draw(st.booleans()):
        prefix += f"{{{draw(st.sampled_from(LABILE_MODS))}}}"
    if draw(st.booleans()):
        prefix += f"[{draw(st.sampled_from(NTERM_MODS))}]-"
    suffix = ""
    if draw(st.booleans()):
        suffix += f"-[{draw(st.sampled_from(CTERM_MODS))}]"
    charges = ["", "/1", "/2", "/3"] + ([] if fragmentable else ["/[Na:z+1^2]", "/[Na:z+1,H:z+1]", "/-1"])
    if "<D>" in prefix:
        # Spec gap: with every H replaced by D, a negative charge removes a light proton
        # that is no longer there (InvalidAdjustmentError). Whether it should remove a
        # deuteron instead is not defined by ProForma; reported, not tested.
        charges.remove("/-1") if "/-1" in charges else None
    return prefix + body + suffix + draw(st.sampled_from(charges))


# --------------------------------------------------------------------------- 1, 2: round-trip


@given(proforma())
def test_parse_serialize_roundtrip(s):
    annot = pt.parse(s)
    again = pt.parse(annot.serialize())
    assert again == annot
    assert again.mass() == pytest.approx(annot.mass(), abs=TOL)


@given(proforma())
def test_serialize_is_idempotent(s):
    once = pt.parse(s).serialize()
    assert pt.parse(once).serialize() == once


# --------------------------------------------------------------------------- 3: fast_fragment == frag

FAST_IONS = ["a", "b", "c", "x", "y", "z"]


@given(proforma(fragmentable=True), st.sets(st.sampled_from(FAST_IONS), min_size=1), st.sets(st.integers(1, 3), min_size=1))
def test_fast_fragment_matches_frag(s, ions, charges):
    annot = pt.parse(s)
    n = len(annot)
    ref: dict[tuple[str, int], list[float]] = {}
    for ion in sorted(ions):
        for charge in sorted(charges):
            try:
                ref[(ion, charge)] = [annot.frag(ion_type=ion, charge=charge, position=p).mz for p in range(1, n + 1)]
            except pt.InvalidAdjustmentError:
                # e.g. a1 of "<13C>A-[Amidated]": the full-length a ion keeps the
                # C-terminal mod and its composition goes negative. fast_fragment must
                # fail the same way rather than return a number.
                with pytest.raises(pt.InvalidAdjustmentError):
                    annot.fast_fragment(ion_types=[ion], charges=[charge])
                return
    fast = annot.fast_fragment(ion_types=sorted(ions), charges=sorted(charges))
    for (ion, charge), mzs in fast.items():
        expected = ref[(str(ion.value) if hasattr(ion, "value") else ion, charge)]
        assert len(mzs) == n
        assert mzs == pytest.approx(expected, abs=TOL), (ion, charge)


# --------------------------------------------------------------------------- 4: b + y = precursor


@given(proforma(fragmentable=True))
def test_complementary_b_y_sum_to_precursor(s):
    annot = pt.parse(s)
    n = len(annot)
    assume(n >= 2)
    # Neutral b_i + y_(n-i) = neutral precursor (the residue sum plus water, plus any
    # terminal, static and isotope modifications).
    precursor = annot.mass(charge=0)
    for i in range(1, n):
        b = annot.frag(ion_type="b", charge=1, position=i).mz - PROTON
        y = annot.frag(ion_type="y", charge=1, position=n - i).mz - PROTON
        assert b + y == pytest.approx(precursor, abs=TOL), i


# --------------------------------------------------------------------------- 5: digestion

ENZYMES = ["trypsin", "lys_c", "arg_c", "glu_c", "asp_n", "chymotrypsin", "lys_n"]


@given(st.text(alphabet=AAS, min_size=1, max_size=60), st.sampled_from(ENZYMES))
def test_digest_spans_reassemble_protein(protein, enzyme):
    annot = pt.parse(protein)
    spans = sorted(annot.digest_spans(enzyme, missed_cleavages=0))
    assert spans, "a digest with no length filter returns at least the whole protein"
    assert spans[0].start == 0
    assert spans[-1].end == len(protein)
    for left, right in zip(spans, spans[1:], strict=False):
        assert left.end == right.start
    assert "".join(protein[sp.start : sp.end] for sp in spans) == protein


@given(st.text(alphabet=AAS, min_size=1, max_size=40), st.sampled_from(ENZYMES), st.integers(0, 2))
def test_missed_cleavage_spans_join_zero_mc_spans(protein, enzyme, mc):
    annot = pt.parse(protein)
    base = sorted(annot.digest_spans(enzyme, missed_cleavages=0))
    cuts = {sp.start for sp in base} | {sp.end for sp in base}
    for sp in annot.digest_spans(enzyme, missed_cleavages=mc):
        assert sp.start in cuts and sp.end in cuts
        inner = sum(1 for c in cuts if sp.start < c < sp.end)
        assert inner == sp.missed_cleavages <= mc


# --------------------------------------------------------------------------- 6: invalid input

TYPED_ERRORS = (pt.ProFormaFormatError, pt.UnsupportedOperationError, pt.UnknownModificationError, pt.CompositionError, pt.InvalidAdjustmentError)
BREAKERS = ["[", "]", "{", "}", "(", ")", "<", ">", "@", "#", "^", "?", "%", "1", "a", " "]


@given(proforma(), st.data())
def test_mutated_proforma_raises_only_typed_errors(s, data):
    pos = data.draw(st.integers(0, len(s)))
    mutated = s[:pos] + data.draw(st.sampled_from(BREAKERS)) + s[pos:]
    # Parsing is partly lazy: a string can parse and fail only when its mass is needed.
    # Either way the error must be one of peptacular's typed errors.
    try:
        pt.parse(mutated).mass()
    except TYPED_ERRORS:
        pass


@pytest.mark.parametrize("bad", ["PEPTIDE/[+H]", "PEPTIDE/[Na]", "PEPTIDE/[Na:z+1^x]"])
def test_bad_charge_adduct_raises_typed_error(bad):
    with pytest.raises(TYPED_ERRORS):
        pt.parse(bad).mass()


@given(proforma(), st.sampled_from(["[", "{", "(", "<"]))
def test_unclosed_bracket_is_a_format_error(s, opener):
    with pytest.raises(pt.ProFormaFormatError):
        pt.parse(s + opener)
