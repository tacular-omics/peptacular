"""Named modifications keep their listed mass on loss and isotope ions.

A loss or isotope peak must be exactly the plain ion plus its delta: switching to the
elemental composition for those ions would move a named mod off its listed database mass
(Oxidation is listed as 15.994915, its composition O gives 15.9949146).
"""

import pytest

import peptacular as pt

NAMED = [
    "PEM[Oxidation]TIDEK",
    "PEPS[Phospho]TIDEK",
    "[TMT6plex]-PEPTIDEK[TMT6plex]",
    "PEM[UNIMOD:35]TIDEK",
    "PEM[MOD:00719]TIDEK",
    "PEPN[Glycan:HexNAc]TIDEK",
    "<[Carbamidomethyl]@C>PEPCTIDEK",
]

H2O = pt.chem_mass("H2O")
NH3 = pt.chem_mass("NH3")


def _by_key(fragments):
    return {(f.ion_type, f.position, f.charge_state): f for f in fragments}


@pytest.mark.parametrize("sequence", NAMED)
@pytest.mark.parametrize("monoisotopic", [True, False])
@pytest.mark.parametrize(("loss", "loss_mass"), [("H2O", H2O), ("NH3", NH3), (-18.010565, None)])
def test_loss_is_plain_ion_minus_loss(sequence, monoisotopic, loss, loss_mass):
    annotation = pt.parse(sequence)
    kwargs = {"ion_types": ["b", "y", "p"], "charges": [1, 2], "monoisotopic": monoisotopic}
    plain = _by_key(annotation.fragment(**kwargs))
    lossy = annotation.fragment(deltas=[loss], **kwargs)
    assert lossy
    if loss_mass is None:
        shift = -loss
    elif monoisotopic:
        shift = loss_mass
    else:
        shift = pt.chem_mass(loss, monoisotopic=False)
    for fragment in lossy:
        reference = plain[(fragment.ion_type, fragment.position, fragment.charge_state)]
        assert fragment.mass + shift == pytest.approx(reference.mass, abs=1e-9, rel=0), fragment


@pytest.mark.parametrize("sequence", NAMED)
def test_neutral_loss_is_plain_ion_minus_loss(sequence):
    annotation = pt.parse(sequence)
    kwargs = {"ion_types": ["b", "y"], "charges": [1]}
    plain = _by_key(annotation.fragment(**kwargs))
    lossy = [f for f in annotation.fragment(neutral_deltas=["H2O"], **kwargs) if f.deltas]
    assert lossy
    for fragment in lossy:
        reference = plain[(fragment.ion_type, fragment.position, fragment.charge_state)]
        assert fragment.mass + H2O == pytest.approx(reference.mass, abs=1e-9, rel=0), fragment


@pytest.mark.parametrize("sequence", NAMED)
@pytest.mark.parametrize("n", [1, 2, 3])
def test_isotope_peak_is_plain_ion_plus_neutrons(sequence, n):
    annotation = pt.parse(sequence)
    kwargs = {"ion_types": ["b", "y", "p"], "charges": [1, 2]}
    plain = _by_key(annotation.fragment(**kwargs))
    heavy = annotation.fragment(isotopes=[n], **kwargs)
    assert len(heavy) == len(plain)
    for fragment in heavy:
        reference = plain[(fragment.ion_type, fragment.position, fragment.charge_state)]
        assert fragment.mass == pytest.approx(reference.mass + n * pt.C13_NEUTRON_MASS, abs=1e-9, rel=0), fragment


@pytest.mark.parametrize("sequence", NAMED)
def test_frag_matches_fragment_with_loss(sequence):
    annotation = pt.parse(sequence)
    series = _by_key(annotation.fragment(ion_types=["y"], charges=[1], deltas=["H2O"]))
    for (ion_type, position, charge), fragment in series.items():
        single = annotation[len(annotation) - position :].frag(ion_type, charge, deltas="H2O")
        assert single.mass == pytest.approx(fragment.mass, abs=1e-9, rel=0)
        assert annotation[len(annotation) - position :].mass(ion_type=ion_type, charge=charge, deltas="H2O") == pytest.approx(fragment.mass, abs=1e-9, rel=0)


@pytest.mark.parametrize("sequence", [s for s in NAMED if "Glycan" not in s])
def test_fast_fragment_matches_fragment(sequence):
    annotation = pt.parse(sequence)
    fast = annotation.fast_fragment(ion_types=["b", "y"], charges=[1, 2])
    for (ion_type, charge), mzs in fast.items():
        slow = [f.mz for f in annotation.fragment(ion_types=[ion_type], charges=[charge])]
        assert mzs == pytest.approx(slow, abs=1e-9, rel=0)


def test_composition_mode_still_uses_composition():
    annotation = pt.parse("PEM[Oxidation]TIDEK")
    by_comp = annotation.mass(deltas="H2O", calculate_with_composition=True)
    assert by_comp == pytest.approx(pt.parse("PEMTIDEK").mass(deltas="H2O", calculate_with_composition=True) + pt.chem_mass("O"), abs=1e-9, rel=0)


def test_labile_mods_leave_fragments_and_stay_on_precursor():
    labile = pt.parse("{Glycan:Hex}PEPTIDEK")
    bare = pt.parse("PEPTIDEK")
    hex_mass = pt.parse("{Glycan:Hex}PEPTIDEK").mass() - bare.mass()
    assert hex_mass == pytest.approx(162.0528234, abs=1e-6)
    for kwargs in ({}, {"deltas": ["H2O"]}, {"isotopes": [1]}):
        fragments = labile.fragment(ion_types=["b", "y"], charges=[1], **kwargs)
        expected = bare.fragment(ion_types=["b", "y"], charges=[1], **kwargs)
        assert [f.mass for f in fragments] == pytest.approx([f.mass for f in expected], abs=1e-9, rel=0)
    precursor = labile.fragment(ion_types=["p"], charges=[1])[0]
    assert precursor.mass == pytest.approx(bare.fragment(ion_types=["p"], charges=[1])[0].mass + hex_mass, abs=1e-9, rel=0)
    assert labile.mass(deltas="H2O") == pytest.approx(labile.mass() - H2O, abs=1e-9, rel=0)


def test_negative_charge_uses_listed_masses():
    # A deprotonating charge takes the composition branch; named mods must still use their
    # listed masses, so negative mode mirrors positive mode.
    from peptacular.constants import PROTON_MASS

    annotation = pt.parse("PEM[Oxidation]TIDEK")
    assert annotation.mass(charge=-2) == pytest.approx(975.4230103537, abs=1e-9, rel=0)
    assert annotation.mass(charge=-2) == pytest.approx(annotation.mass() - 2 * PROTON_MASS, abs=1e-9, rel=0)
    assert pt.parse("PEM[Oxidation]TIDEK/-2").mass() == pytest.approx(annotation.mass(charge=-2), abs=1e-9, rel=0)
    precursor = annotation.fragment(ion_types=["p"], charges=[-2])[0]
    assert precursor.mass == pytest.approx(annotation.mass(charge=-2), abs=1e-9, rel=0)
