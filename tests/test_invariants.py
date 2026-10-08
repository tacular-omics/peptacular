"""Invariant tests: relations that must hold for every input, not hand-picked examples.

Three groups, each extending ``test_hypothesis_properties.py`` rather than repeating it:

1. ProForma round trip (Hypothesis) over a wider grammar than the base strategy there:
   every modification vocabulary and tag form, cross-link and group labels with scores,
   unknown-position counts, multiple labile and terminal mods, fixed mods on residues and
   termini, isotope labels, charge carriers of both signs (mixed, with occurrence counts)
   and chimeric (multi-chain) strings.
2. Mass relations over a peptide x ion type x charge (-3..+3) x carrier grid.
3. Data-driven loops over every residue code in tacular's amino-acid table and every ion
   type tacular defines.

Example counts follow the Hypothesis profile in ``tests/conftest.py``.
"""

from __future__ import annotations

import math

import pytest
from hypothesis import given
from hypothesis import strategies as st
from tacular import AA_LOOKUP, IonType

import peptacular as pt

TOL = 1e-6
ELECTRON = 0.000548579909065
# Monoisotopic element masses (AME 2020, as in tacular) for the carrier formulas below.
H, N, NA, K, CL = 1.00782503223, 14.00307400443, 22.989769282, 38.963706486, 34.968852682
N15 = 15.00010889888

TYPED_ERRORS = (pt.PeptacularError, pt.ProFormaFormatError, pt.UnsupportedOperationError, pt.UnknownModificationError, pt.CompositionError)


def _outcome(fn):
    """A value, or the exception type it raised (so two paths can be compared either way)."""
    try:
        return fn()
    except TYPED_ERRORS as exc:
        return type(exc)


def _same(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return a == pytest.approx(b, abs=TOL)
    return a == b


# =========================================================================== 1. round trip

# Residue tags: Unimod / PSI-MOD names with and without prefix, accessions in every CV that
# resolves (Unimod, PSI-MOD, RESID, GNO, XL-MOD), delta masses with and without a CV or
# Obs prefix, formulas (with isotopes), glycans (with formula and mass components) and
# INFO tags as alternatives.
RESIDUE_TAGS = [
    "Oxidation",
    "U:Phospho",
    "M:O-phospho-L-serine",
    "monohydroxylated residue",
    "UNIMOD:35",
    "MOD:00046",
    "RESID:AA0037",
    "GNO:G59626AS",
    "G:G59626AS",
    "XLMOD:02001",
    "+15.995",
    "-17.0265",
    "U:+79.966",
    "M:+79.966",
    "Obs:+12.3",
    "Formula:HPO3",
    "Formula:[13C2]C-2H2",
    "Glycan:HexNAc1Hex2",
    "Glycan:Hex{H2O}{+204.068}",
    "Phospho|INFO:a note",
    "+79.966|Obs:+79.978",
    "Oxidation|UNIMOD:35",
]
NTERM_TAGS = ["Acetyl", "UNIMOD:1", "Formula:C2H2O", "+42.011", "TMT6plex"]
CTERM_TAGS = ["Amidated", "Methyl", "-0.984"]
LABILE_TAGS = ["Glycan:Hex1HexNAc1", "Glycan:Hex", "Phospho", "Formula:HPO3"]
UNKNOWN_TAGS = ["Phospho", "Oxidation", "+12.0", "Formula:C2H2O"]
FIXED_MODS = ["<[Carbamidomethyl]@C>", "<[Oxidation]@M>", "<[+1.0]@K,R>", "<[Acetyl]@N-term>", "<[Amidated]@C-term>", "<[TMT6plex]@K,N-term>"]
ISOTOPES = ["<13C>", "<15N>", "<18O>", "<D>"]

# Charges: bare integers of both signs and carrier lists (proton, sodium, potassium,
# ammonium with a 15N label, a zinc dication, deprotonation and chloride), with
# occurrence counts and mixed signs.
CHARGES = [
    "",
    "/1",
    "/3",
    "/-1",
    "/-3",
    "/[H:z+1^2]",
    "/[Na:z+1]",
    "/[Na:z+1,H:z+1]",
    "/[K:z+1^2,H:z+1]",
    "/[[15N1]H4:z+1]",
    "/[Zn:z+2]",
    "/[H-1:z-1^2]",
    "/[Cl:z-1]",
    "/[Cl:z-1,H-1:z-1]",
    "/[Na:z+1,H-1:z-1^2]",
]
NEGATIVE_CHARGES = {c for c in CHARGES if c.startswith("/-") or "H-1:" in c}

AAS = "ACDEFGHIKLMNPQRSTVWY"


@st.composite
def chain(draw, *, allow_globals: bool = True) -> str:
    """One peptidoform ion (a single chain) in ProForma 2.1 notation.

    ``allow_globals=False`` leaves out isotope labels and fixed mods, which in a chimeric
    string may only open the whole string and then apply to every chain.
    """
    seq = draw(st.text(alphabet=AAS, min_size=2, max_size=12))
    n = len(seq)
    tags: list[str] = [""] * n
    for i in range(n):
        if draw(st.integers(0, 4)) == 0:
            tags[i] = "".join(f"[{t}]" for t in draw(st.lists(st.sampled_from(RESIDUE_TAGS), min_size=1, max_size=2)))

    # At most one positional label group: an intrachain cross-link, or an ambiguity
    # group with optional localisation scores.
    group_kind = draw(st.sampled_from([None, "xl", "group", "scored"]))
    if group_kind is not None and n >= 3:
        sites = sorted(draw(st.lists(st.integers(0, n - 1), min_size=2, max_size=3, unique=True)))
        if group_kind == "xl":
            sites = sites[:2]
            tags[sites[0]] += "[X:DSS#XL1]"
            tags[sites[1]] += "[#XL1]"
        else:
            tag = draw(st.sampled_from(["Phospho", "Oxidation", "+15.995"]))
            scores = [""] * len(sites)
            if group_kind == "scored":
                scores = [f"({s})" for s in ("0.75", "0.2", "0.05")[: len(sites)]]
            tags[sites[0]] += f"[{tag}#g1{scores[0]}]"
            for site, score in zip(sites[1:], scores[1:], strict=True):
                tags[site] += f"[#g1{score}]"

    residues = [aa + tag for aa, tag in zip(seq, tags, strict=True)]

    # A range (with a mod) or an ambiguous span over residues that carry no tag.
    span_kind = draw(st.sampled_from([None, "range", "ambiguous"]))
    if span_kind is not None:
        start = draw(st.integers(0, n - 2))
        width = draw(st.integers(2, min(4, n - start)))
        if not any(tags[start : start + width]):
            inner = seq[start : start + width]
            if span_kind == "range":
                block = f"({inner})[{draw(st.sampled_from(['+19.0523', 'Oxidation', 'Formula:H2O']))}]"
            else:
                block = f"(?{inner})"
            residues[start : start + width] = [block] + [""] * (width - 1)
    body = "".join(residues)

    prefix = ""
    isotope = draw(st.lists(st.sampled_from(ISOTOPES), max_size=2, unique=True)) if allow_globals else []
    prefix += "".join(isotope)
    if allow_globals:
        prefix += "".join(draw(st.lists(st.sampled_from(FIXED_MODS), max_size=2, unique=True)))
    unknown = draw(st.lists(st.sampled_from(UNKNOWN_TAGS), max_size=2, unique=True))
    if unknown:
        prefix += "".join(f"[{t}]" + (f"^{draw(st.integers(2, 3))}" if draw(st.booleans()) else "") for t in unknown) + "?"
    prefix += "".join(f"{{{t}}}" for t in draw(st.lists(st.sampled_from(LABILE_TAGS), max_size=2)))
    nterm = draw(st.lists(st.sampled_from(NTERM_TAGS), max_size=2, unique=True))
    if nterm:
        prefix += "".join(f"[{t}]" for t in nterm) + "-"
    cterm = draw(st.lists(st.sampled_from(CTERM_TAGS), max_size=2, unique=True))
    suffix = "-" + "".join(f"[{t}]" for t in cterm) if cterm else ""

    options = list(CHARGES)
    if "<D>" in isotope:
        # Spec gap (see test_hypothesis_properties.py): with every H replaced by D a
        # deprotonation removes a light proton that is not there.
        options = [c for c in options if c not in NEGATIVE_CHARGES]
    return prefix + body + suffix + draw(st.sampled_from(options))


def _mass_profile(annot: pt.ProFormaAnnotation):
    return (
        _outcome(lambda: annot.mass()),
        _outcome(lambda: annot.mz()),
        _outcome(lambda: annot.mass(charge=0)),
        annot.charge_state,
    )


@given(chain())
def test_roundtrip_single_chain(s):
    annot = pt.parse(s)
    out = annot.serialize()
    again = pt.parse(out)
    assert again == annot, (s, out)
    assert again.serialize() == out, "serialization is not idempotent"
    for a, b in zip(_mass_profile(annot), _mass_profile(again), strict=True):
        assert _same(a, b), (s, out, a, b)


@given(chain())
def test_roundtrip_mass_is_finite_or_typed_error(s):
    # Every mass of a valid string is a finite number or one of peptacular's typed errors.
    for value in _mass_profile(pt.parse(s))[:3]:
        if isinstance(value, float):
            assert math.isfinite(value), s
        else:
            assert issubclass(value, TYPED_ERRORS), (s, value)


GLOBAL_PREFIXES = ["", "<13C>", "<15N>", "<[Carbamidomethyl]@C>", "<13C><[Oxidation]@M>"]


@given(st.sampled_from(GLOBAL_PREFIXES), st.lists(chain(allow_globals=False), min_size=2, max_size=3))
def test_roundtrip_chimeric(globals_, chains):
    s = globals_ + "+".join(chains)
    parts = list(pt.parse_chimeric(s))
    assert len(parts) == len(chains)
    out = pt.serialize_chimeric(parts)
    again = list(pt.parse_chimeric(out))
    assert again == parts
    assert pt.serialize_chimeric(again) == out
    for single, a, b in zip(chains, parts, again, strict=True):
        # Each chain of a chimeric string is the same peptidoform as on its own, with
        # the leading global mods applied to every chain.
        assert a == pt.parse(globals_ + single)
        for x, y in zip(_mass_profile(a), _mass_profile(b), strict=True):
            assert _same(x, y), (single, x, y)


@pytest.mark.parametrize(
    "s",
    [
        # The 2.0 adduct form "/-1[+e-]" is not ProForma 2.1 and an electron is not a
        # chemical formula: parsing is lazy, but a mass must be a typed error, never a number.
        "PEPTIDE/[e-:z-1]",
        "PEPTIDE/[H:z+1,e-:z-1]",
        # A CV prefix before an accession number is incorrect usage per the spec.
        "PEP[R:AA0037]TIDE",
    ],
)
def test_invalid_carriers_and_tags_fail_typed(s):
    with pytest.raises(TYPED_ERRORS):
        pt.parse(s).mass()


# =========================================================================== 2. mass grid

GRID_PEPTIDES = [
    "PEPTIDEK",
    "[Acetyl]-M[Oxidation]PEPS[Phospho]K",
    "<13C>ACDEFGHIK",
    "<[Carbamidomethyl]@C>CGCPEPC-[Amidated]",
    "GN[Formula:C8H13NO5]GTR",
    "SEQ[+15.995]UENCEO",
]
# (carrier for one positive charge, its mass delta) and the same for one negative charge.
POS_CARRIERS = {
    "H:z+1": H - ELECTRON,
    "Na:z+1": NA - ELECTRON,
    "K:z+1": K - ELECTRON,
    "[15N1]H4:z+1": N15 + 4 * H - ELECTRON,
}
NEG_CARRIERS = {
    "H-1:z-1": -H + ELECTRON,
    "Cl:z-1": CL + ELECTRON,
}
GRID_CHARGES = [-3, -2, -1, 1, 2, 3]
ALL_IONS = list(IonType)
FAST_IONS = ["a", "b", "c", "x", "y", "z", "p", "n"]


def _carrier_specs(z: int):
    """Carrier lists of total charge z: one carrier type repeated, and a mixed list."""
    pool = POS_CARRIERS if z > 0 else NEG_CARRIERS
    specs = [([name] * abs(z), delta * abs(z)) for name, delta in pool.items()]
    names = list(pool)
    mixed = [names[i % len(names)] for i in range(abs(z))]
    specs.append((mixed, sum(pool[c] for c in mixed)))
    return specs


@pytest.mark.parametrize("seq", GRID_PEPTIDES)
@pytest.mark.parametrize("z", GRID_CHARGES)
def test_precursor_mass_grid(seq, z):
    annot = pt.parse(seq)
    neutral = annot.mass(charge=0)
    assert annot.neutral_mass() == pytest.approx(neutral, abs=TOL)
    # An integer charge means protons (z > 0) or deprotonation (z < 0).
    proton_delta = (POS_CARRIERS["H:z+1"] if z > 0 else NEG_CARRIERS["H-1:z-1"]) * abs(z)
    assert annot.mass(charge=z) == pytest.approx(neutral + proton_delta, abs=TOL)
    for carriers, delta in _carrier_specs(z):
        charged = annot.mass(charge=carriers)
        assert charged == pytest.approx(neutral + delta, abs=TOL), carriers
        assert annot.mz(charge=carriers) * abs(z) == pytest.approx(charged, abs=TOL), carriers
        # Encoded in the string, the carriers give the same numbers and the same charge.
        encoded = pt.parse(seq + "/[" + ",".join(carriers) + "]")
        assert encoded.charge_state == z
        assert encoded.mass() == pytest.approx(charged, abs=TOL)


@pytest.mark.parametrize("seq", GRID_PEPTIDES)
def test_fragment_mass_grid(seq):
    annot = pt.parse(seq)
    neutral = {ion: {f.position: f for f in annot.fragment([ion], [0])} for ion in ALL_IONS}
    for z in GRID_CHARGES:
        for carriers, delta in [(z, None), *_carrier_specs(z)]:
            for ion in ALL_IONS:
                charged = {f.position: f for f in annot.fragment([ion], [carriers])}
                assert set(charged) <= set(neutral[ion]), (ion, carriers)
                for pos in set(neutral[ion]) - set(charged):
                    # Documented skip: an ion too small to lose the atoms a deprotonation
                    # removes (az of one D at -3 has two H). frag() raises for it.
                    assert z < 0, (ion, carriers, pos)
                    with pytest.raises(pt.InvalidAdjustmentError):
                        annot.frag(ion, carriers, position=pos)
                for pos, f in charged.items():
                    assert f.charge_state == z, (ion, carriers, pos)
                    assert math.isfinite(f.mz) and f.mz > 0, (ion, carriers, pos, f.mz)
                    assert f.mz * abs(z) == pytest.approx(f.mass, abs=TOL)
                    if delta is not None:
                        assert f.mass == pytest.approx(neutral[ion][pos].mass + delta, abs=TOL), (ion, carriers, pos)


@pytest.mark.parametrize("seq", GRID_PEPTIDES)
def test_fast_fragment_matches_fragment_all_charges(seq):
    # Extends test_fast_fragment_matches_frag (charges 1..3, a-z) to negative charges and
    # the p and n series.
    annot = pt.parse(seq)
    fast = annot.fast_fragment(FAST_IONS, GRID_CHARGES)
    for (ion, z), mzs in fast.items():
        series = annot.fragment([ion], [z])
        if ion in (IonType.PRECURSOR, IonType.NEUTRAL):
            # One whole-molecule ion, repeated per position by fast_fragment.
            assert len(series) == 1
            assert mzs == pytest.approx([series[0].mz] * len(annot), abs=TOL), (ion, z)
        else:
            by_pos = {f.position: f.mz for f in series}
            assert mzs == pytest.approx([by_pos[p] for p in range(1, len(annot) + 1)], abs=TOL), (ion, z)


@pytest.mark.parametrize("seq", GRID_PEPTIDES)
def test_complementary_ions_sum_is_constant_any_charge(seq):
    # Neutralised, b_i + y_(n-i) is the precursor; a_i + x_(n-i) and c_i + z_(n-i) are each
    # a constant (the precursor minus a fixed offset), whatever the cleavage site and charge.
    annot = pt.parse(seq)
    n = len(annot)
    precursor = annot.mass(charge=0)
    sums: dict[tuple[str, str], list[float]] = {}
    for z in GRID_CHARGES:
        proton = (POS_CARRIERS["H:z+1"] if z > 0 else NEG_CARRIERS["H-1:z-1"]) * abs(z)
        for left, right in (("b", "y"), ("a", "x"), ("c", "z")):
            fast = annot.fast_fragment([left, right], [z])
            lm, rm = fast[(IonType(left), z)], fast[(IonType(right), z)]
            sums.setdefault((left, right), []).extend((lm[i - 1] * abs(z) - proton) + (rm[n - i - 1] * abs(z) - proton) for i in range(1, n))
    assert sums[("b", "y")] == pytest.approx([precursor] * len(sums[("b", "y")]), abs=1e-5)
    for pair, values in sums.items():
        assert values == pytest.approx([values[0]] * len(values), abs=1e-5), pair


# =========================================================================== 3. data loops

RESIDUE_CODES = sorted(code for code in AA_LOOKUP.keys() if isinstance(code, str) and len(code) == 1)
# B (Asx), Z (Glx) and X (any) have no single mass; documented to raise
# "Mass not available for amino acid" (X since the bare-X fix in [Unreleased]).
NO_MASS_CODES = {code for code in RESIDUE_CODES if AA_LOOKUP[code].monoisotopic_mass is None} | {"X"}


def _placements(aa: str) -> list[str]:
    return [aa, aa + "GEK", "GE" + aa + "AK", "GEK" + aa, aa * 3]


def test_residue_table_includes_extended_codes():
    assert {"U", "O", "X", "B", "Z", "J"} <= set(RESIDUE_CODES)
    assert NO_MASS_CODES == {"B", "Z", "X"}


@pytest.mark.parametrize("aa", RESIDUE_CODES)
def test_every_residue_mass_paths(aa):
    for seq in _placements(aa):
        calls = {
            "mass": lambda s=seq: pt.mass(s),
            "mz": lambda s=seq: pt.mz(s, charge=2),
            "mz_neg": lambda s=seq: pt.mz(s, charge=-1),
            "comp": lambda s=seq: pt.comp(s),
            "fragment": lambda s=seq: pt.fragment(s, ["b", "y"], [1]),
            "fast_fragment": lambda s=seq: pt.fast_fragment(s, ["b", "y"], [1]),
        }
        for name, call in calls.items():
            if aa in NO_MASS_CODES:
                with pytest.raises((pt.PeptacularError, pt.CompositionError), match=f"not available for amino acid: {aa}"):
                    call()
                continue
            value = call()
            if isinstance(value, float):
                assert math.isfinite(value) and value > 0, (seq, name, value)
            elif isinstance(value, dict) and name == "fast_fragment":
                assert all(math.isfinite(m) for mzs in value.values() for m in mzs), (seq, name)
            elif name == "fragment":
                assert value and all(math.isfinite(f.mz) for f in value), (seq, name)


@pytest.mark.parametrize("aa", RESIDUE_CODES)
def test_every_residue_pi_and_charge(aa):
    # pI and charge-at-pH ignore masses, so every residue code (B, Z, X, J too) works.
    for seq in _placements(aa):
        charges = [pt.charge_at_ph(seq, pH=ph) for ph in (0.0, 2.0, 4.0, 7.0, 10.0, 12.0, 14.0)]
        assert all(math.isfinite(c) for c in charges), seq
        # Net charge never rises with pH, starts positive and ends negative.
        assert all(a >= b - 1e-12 for a, b in zip(charges, charges[1:], strict=False)), (seq, charges)
        assert charges[0] > 0 > charges[-1], (seq, charges)
        p = pt.pi(seq)
        assert math.isfinite(p) and 0.0 < p < 14.0, (seq, p)
        # pI is where the charge crosses zero.
        assert pt.charge_at_ph(seq, pH=p - 0.01) > 0 > pt.charge_at_ph(seq, pH=p + 0.01), (seq, p)


MASS_CODES = [code for code in RESIDUE_CODES if code not in NO_MASS_CODES]


@pytest.mark.parametrize("ion", ALL_IONS, ids=lambda i: i.value)
def test_every_ion_type_finite_mz(ion):
    # Every ion type tacular defines, on a peptide holding every mass-bearing residue code
    # (each at N-term, C-term and inside somewhere across the three sequences).
    body = "".join(MASS_CODES)
    for seq in (body, body[::-1], "K" + body[len(body) // 2 :] + body[: len(body) // 2] + "R"):
        for charge in (1, 2, -1, "Na:z+1"):
            frags = pt.fragment(seq, [ion], [charge])
            for f in frags:
                assert math.isfinite(f.mz) and f.mz > 0, (seq, ion, charge, f.position, f.mz)
                assert math.isfinite(f.mass), (seq, ion, charge, f.position)
