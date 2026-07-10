"""
Cross-validation against pyteomics, an independent, widely-used proteomics
library (a peptacular dev dependency), to catch mass/digestion regressions
that peptacular's own tests might not surface.
"""

import pytest
from pyteomics import mass as pmass
from pyteomics import parser as pparser

import peptacular as pt

MASS_ABS_TOL = 1e-4

PEPTIDES = [
    "PEPTIDE",
    "ACDEFGHIKLMNPQRSTVWY",
    "MVIMSEFSADPAGQGQGQQK",
    "SAMPLER",
    "K",
    "GG",
]


class TestNeutralMass:
    @pytest.mark.parametrize("seq", PEPTIDES)
    def test_neutral_monoisotopic_mass_matches(self, seq):
        expected = pmass.calculate_mass(sequence=seq)
        actual = pt.mass(seq)
        assert actual == pytest.approx(expected, abs=MASS_ABS_TOL)

    @pytest.mark.parametrize("seq", PEPTIDES)
    @pytest.mark.parametrize("charge", [1, 2, 3])
    def test_charged_mz_matches(self, seq, charge):
        expected = pmass.calculate_mass(sequence=seq, charge=charge)
        actual = pt.mz(seq, charge=charge)
        assert actual == pytest.approx(expected, abs=MASS_ABS_TOL)

    def test_standard_residue_masses_match_std_aa_mass(self):
        water = pmass.calculate_mass(formula="H2O")
        for aa, residue_mass in pmass.std_aa_mass.items():
            if len(aa) != 1:
                continue
            actual_residue_mass = pt.mass(aa) - water
            assert actual_residue_mass == pytest.approx(residue_mass, abs=MASS_ABS_TOL), aa


class TestProFormaMass:
    @pytest.mark.parametrize("charge", [None, 1, 2])
    def test_delta_mass_modification_matches(self, charge):
        seq = "PEPTIC[+57.021464]IDE"
        expected = pmass.calculate_mass(proforma=seq, charge=charge)
        actual = pt.mz(seq, charge=charge) if charge else pt.mass(seq)
        assert actual == pytest.approx(expected, abs=MASS_ABS_TOL)

    @pytest.mark.parametrize(
        "seq",
        [
            "NEEYN[Glycan:Hex5HexNAc4NeuAc1]K",
            "NEEYN[GNO:G59626AS]K",
        ],
    )
    def test_glycan_mass_matches_pyteomics_gno_reference(self, seq):
        # Reference value asserted by pyteomics' own test suite
        # (tests/test_proforma.py::test_gnome / test_glycan).
        assert pt.mass(seq) == pytest.approx(2709.016, abs=1e-2)


class TestDigestion:
    PROTEIN = "MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTE"

    def test_pyteomics_own_trypsin_example(self):
        # tests/test_parser.py::test_cleave in pyteomics
        assert pparser.xcleave("PEPTIDEKS", pparser.expasy_rules["trypsin"]) == [(0, "PEPTIDEK"), (8, "S")]
        result = pt.digest("PEPTIDEKS", pt.Proteases.TRYPSIN, missed_cleavages=0)
        assert result == [("PEPTIDEK", pt.Span(0, 8, 0)), ("S", pt.Span(8, 9, 0))]

    def test_pyteomics_own_semi_tryptic_example(self):
        # tests/test_parser.py::test_cleave_semi in pyteomics
        expected = {
            "PEPTIDEK",
            "P",
            "PE",
            "PEP",
            "PEPT",
            "PEPTI",
            "PEPTID",
            "PEPTIDE",
            "EPTIDEK",
            "PTIDEK",
            "TIDEK",
            "IDEK",
            "DEK",
            "EK",
            "K",
            "S",
        }
        result = {p for p, _ in pt.digest("PEPTIDEKS", pt.Proteases.TRYPSIN, missed_cleavages=0, semi=True)}
        assert result == expected

    @pytest.mark.parametrize(
        "protease_id",
        [
            "trypsin",
            "trypsin_full",
            "arg_c",
            "asp_n",
            "lys_c",
            "lys_n",
            "glu_c",
            "chymotrypsin",
            "chymotrypsin_low",
            "proteinase_k",
            "thermolysin",
        ],
    )
    @pytest.mark.parametrize("missed_cleavages", [0, 1, 2])
    @pytest.mark.parametrize("semi", [False, True])
    def test_digest_agrees_with_pyteomics_for_same_regex(self, protease_id, missed_cleavages, semi):
        # Feed peptacular's own protease regex into pyteomics.parser.cleave so this
        # isolates agreement on cleavage-site/cut-position semantics, independent of
        # any differences between the two libraries' enzyme-name-to-regex tables.
        regex = pt.PROTEASE_LOOKUP.get(protease_id).regex
        pt_result = {p for p, _ in pt.digest(self.PROTEIN, regex, missed_cleavages=missed_cleavages, semi=semi)}
        py_result = set(pparser.cleave(self.PROTEIN, regex, missed_cleavages=missed_cleavages, semi=semi))
        assert pt_result == py_result

    @pytest.mark.parametrize("min_len,max_len", [(None, None), (5, 20), (7, 15)])
    def test_semi_digest_length_filters_agree_with_pyteomics(self, min_len, max_len):
        regex = pt.PROTEASE_LOOKUP.get("trypsin").regex
        pt_result = {p for p, _ in pt.digest(self.PROTEIN, regex, missed_cleavages=2, semi=True, min_len=min_len, max_len=max_len)}
        py_result = set(pparser.cleave(self.PROTEIN, regex, missed_cleavages=2, semi=True, min_length=min_len))
        if max_len is not None:
            py_result = {p for p in py_result if len(p) <= max_len}
        assert pt_result == py_result
