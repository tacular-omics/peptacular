import pytest

import peptacular as pt
from peptacular import (
    PeptacularError,
    convert_casanovo_sequence,
    convert_diann_sequence,
    convert_ip2_sequence,
)


class TestConvertIP2Sequence:
    """Tests for convert_ip2_sequence function."""

    def test_basic_conversion(self):
        """Test basic IP2 sequence conversion."""
        result = convert_ip2_sequence("K.PEP(phospho)TIDE.K")
        assert result == "PEP[phospho]TIDE"

    def test_nterm_modification(self):
        """Test N-terminal modification conversion."""
        result = convert_ip2_sequence("K.(-1)PEP(phospho)TIDE.K")
        assert result == "[-1]-PEP[phospho]TIDE"

    def test_numeric_modification(self):
        """Test numeric modification conversion."""
        result = convert_ip2_sequence("K.PEPTIDE(2).K")
        assert result == "PEPTIDE[2]"

    def test_multiple_cterm_modifications(self):
        """Test multiple C-terminal modifications."""
        result = convert_ip2_sequence("K.PEPTIDE(2)(3).K")
        assert result == "PEPTIDE[2]-[3]"

    def test_complex_modifications(self):
        """Test complex modification pattern."""
        result = convert_ip2_sequence("-.(1)PEP(phospho)TIDE(2)(3).-")
        assert result == "[1]-PEP[phospho]TIDE[2]-[3]"

    def test_single_amino_acid(self):
        """Test single amino acid."""
        result = convert_ip2_sequence("P")
        assert result == "P"

    def test_empty_string(self):
        """Test empty string."""
        result = convert_ip2_sequence("")
        assert result == ""

    def test_no_flanking_residues(self):
        """Test sequence without flanking residues."""
        result = convert_ip2_sequence("PEPTIDE")
        assert result == "PEPTIDE"

    def test_parallel_processing(self):
        """Test parallel processing with list of sequences."""
        sequences = [
            "K.PEP(phospho)TIDE.K",
            "K.PEPTIDE(2).K",
            "P",
            "PEPTIDE",
        ]
        expected = [
            "PEP[phospho]TIDE",
            "PEPTIDE[2]",
            "P",
            "PEPTIDE",
        ]
        results = convert_ip2_sequence(sequences, n_workers=2, method="sequential")
        assert results == expected

    def test_parallel_with_thread_method(self):
        """Test parallel processing with threading."""
        sequences = ["K.PEPTIDE.K", "R.SEQUENCE.R"]
        expected = ["PEPTIDE", "SEQUENCE"]
        results = convert_ip2_sequence(sequences, method="thread", n_workers=2)
        assert results == expected

    def test_parallel_with_process_method(self):
        """Test parallel processing with multiprocessing."""
        sequences = ["K.PEPTIDE.K", "R.SEQUENCE.R"]
        expected = ["PEPTIDE", "SEQUENCE"]
        results = convert_ip2_sequence(sequences, method="process", n_workers=2)
        assert results == expected


class TestConvertDiannSequence:
    """Tests for convert_diann_sequence function."""

    def test_basic_conversion(self):
        """Test basic DIANN sequence conversion."""
        result = convert_diann_sequence("_YMGTLRGC[Carbamidomethyl]LLRLYHD_")
        assert result == "YMGTLRGC[Carbamidomethyl]LLRLYHD"

    def test_nterm_modification(self):
        """Test N-terminal modification conversion."""
        result = convert_diann_sequence("_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD_")
        assert result == "[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD"

    def test_cterm_modification(self):
        """Test C-terminal modification conversion."""
        result = convert_diann_sequence("_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]_[Methyl]")
        assert result == "[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]-[Methyl]"

    def test_no_underscores(self):
        """Test sequence without underscores."""
        result = convert_diann_sequence("PEPTIDE")
        assert result == "PEPTIDE"

    def test_only_nterm_underscore(self):
        """Test sequence with only N-terminal underscore."""
        result = convert_diann_sequence("_PEPTIDE")
        assert result == "PEPTIDE"

    def test_only_cterm_underscore(self):
        """Test sequence with only C-terminal underscore."""
        result = convert_diann_sequence("PEPTIDE_")
        assert result == "PEPTIDE"

    def test_parallel_processing(self):
        """Test parallel processing with list of sequences."""
        sequences = [
            "_YMGTLRGC[Carbamidomethyl]LLRLYHD_",
            "_[Acytel]PEPTIDE_",
            "SEQUENCE",
        ]
        expected = [
            "YMGTLRGC[Carbamidomethyl]LLRLYHD",
            "[Acytel]-PEPTIDE",
            "SEQUENCE",
        ]
        results = convert_diann_sequence(sequences, n_workers=2)
        assert results == expected

    def test_parallel_with_thread_method(self):
        """Test parallel processing with threading."""
        sequences = ["_PEPTIDE_", "_SEQUENCE_"]
        expected = ["PEPTIDE", "SEQUENCE"]
        results = convert_diann_sequence(sequences, method="thread", n_workers=2)
        assert results == expected


class TestConvertCasanovoSequence:
    """Tests for convert_casanovo_sequence function."""

    def test_basic_conversion(self):
        """Test basic Casanovo sequence conversion."""
        result = convert_casanovo_sequence("+43.006P+100EPTIDE")
        assert result == "[+43.006]-P[+100]EPTIDE"

    def test_nterm_modification_only(self):
        """Test N-terminal modification only."""
        result = convert_casanovo_sequence("+42.010PEPTIDE")
        assert result == "[+42.010]-PEPTIDE"

    def test_internal_modifications(self):
        """Test internal modifications."""
        result = convert_casanovo_sequence("P+80EPTIDE")
        assert result == "P[+80]EPTIDE"

    def test_negative_modification(self):
        """Test negative mass modification."""
        result = convert_casanovo_sequence("P-18EPTIDE")
        assert result == "P[-18]EPTIDE"

    def test_no_modifications(self):
        """Test sequence without modifications."""
        result = convert_casanovo_sequence("PEPTIDE")
        assert result == "PEPTIDE"

    def test_multiple_modifications(self):
        """Test multiple modifications."""
        result = convert_casanovo_sequence("+43P+16E+80PTIDE")
        assert result == "[+43]-P[+16]E[+80]PTIDE"

    def test_modification_at_end(self):
        """Test modification at sequence end."""
        result = convert_casanovo_sequence("PEPTIDE+10")
        assert result == "PEPTIDE[+10]"

    def test_parallel_processing(self):
        """Test parallel processing with list of sequences."""
        sequences = [
            "+43.006P+100EPTIDE",
            "PEPTIDE",
            "P+80EPTIDE",
        ]
        expected = [
            "[+43.006]-P[+100]EPTIDE",
            "PEPTIDE",
            "P[+80]EPTIDE",
        ]
        results = convert_casanovo_sequence(sequences, n_workers=2)
        assert results == expected

    def test_parallel_with_thread_method(self):
        """Test parallel processing with threading."""
        sequences = ["+43PEPTIDE", "SEQUENCE"]
        expected = ["[+43]-PEPTIDE", "SEQUENCE"]
        results = convert_casanovo_sequence(sequences, method="thread", n_workers=2)
        assert results == expected

    def test_parallel_with_process_method(self):
        """Test parallel processing with multiprocessing."""
        sequences = ["+43PEPTIDE", "P+80EPTIDE"]
        expected = ["[+43]-PEPTIDE", "P[+80]EPTIDE"]
        results = convert_casanovo_sequence(sequences, method="process", n_workers=2)
        assert results == expected


class TestConvertersEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_list(self):
        """Test with empty list."""
        assert convert_ip2_sequence([]) == []
        assert convert_diann_sequence([]) == []
        assert convert_casanovo_sequence([]) == []

    def test_single_item_list(self):
        """Test with single item list."""
        result = convert_ip2_sequence(["K.PEPTIDE.K"])
        assert result == ["PEPTIDE"]

        result = convert_diann_sequence(["_PEPTIDE_"])
        assert result == ["PEPTIDE"]

        result = convert_casanovo_sequence(["+43PEPTIDE"])
        assert result == ["[+43]-PEPTIDE"]

    def test_large_list(self):
        """Test with larger list to ensure parallel processing works."""
        sequences = ["K.PEPTIDE.K"] * 100
        expected = ["PEPTIDE"] * 100
        results = convert_ip2_sequence(sequences, n_workers=4)
        assert results == expected
        assert len(results) == 100


class TestToIP2:
    """Tests for ProFormaAnnotation.to_ip2."""

    @pytest.mark.parametrize(
        ("proforma", "expected"),
        [
            ("PEP[Phospho]TIDE", "PEP(Phospho)TIDE"),
            ("PEP[+79.9663]TIDE", "PEP(+79.9663)TIDE"),
            ("PEP[UNIMOD:21]TIDE", "PEP(UNIMOD:21)TIDE"),
            ("[Acetyl]-PEPTIDE", "(Acetyl)PEPTIDE"),
            ("PEPTIDE[Phospho]-[Amidated]", "PEPTIDE(Phospho)(Amidated)"),
            ("PEPTIDE", "PEPTIDE"),
            ("[A][B]-PEPTIDE", "(A)(B)PEPTIDE"),
            ("PEP[A][B]TIDE", "PEP(A)(B)TIDE"),
            ("[A]-P[B]-[C][D]", "(A)P(B)(C)(D)"),
            ("<[Carbamidomethyl]@C>PECTIDE", "PEC(Carbamidomethyl)TIDE"),
        ],
    )
    def test_output(self, proforma: str, expected: str):
        assert pt.parse(proforma).to_ip2() == expected

    @pytest.mark.parametrize(
        "ip2",
        [
            "PEP(Phospho)TIDE",
            "(Acetyl)PEPTIDE",
            "(Acetyl)(Oxidation)(Phospho)PEPTIDE",
            "PEPTIDE(Phospho)(Amidated)",
            "PEPTIDE(Phospho)(Amidated)(Oxidation)",
            "PEP(Phospho)(Oxidation)TIDE",
            "(Acetyl)PEP(Phospho)TIDE(Oxidation)(Amidated)",
            "K.PEP(Phospho)TIDE.R",
            "-.(Acetyl)PEPTIDE(Phospho)(Amidated).-",
            "K.PEPTIDE.R",
        ],
    )
    def test_round_trip(self, ip2: str):
        original = pt.ProFormaAnnotation.from_ip2_sequence(ip2)
        written = original.to_ip2()
        assert "." not in written
        again = pt.ProFormaAnnotation.from_ip2_sequence(written)
        assert again == original
        assert again.mass() == pytest.approx(original.mass())

    @pytest.mark.parametrize(
        ("proforma", "feature"),
        [
            ("<13C>PEPTIDE", "isotope"),
            ("{Hex}PEPTIDE", "labile"),
            ("[Phospho]?PEPTIDE", "unknown"),
            ("PEP(TI)[Phospho]DE", "interval"),
            ("PEPTIDE/2", "charge"),
            ("PEPTIDE/[+2Na+]", "charge adducts"),
            ("PEPTIDE-[Amidated]", "C-terminal modification on a last residue"),
            ("PEPTIDEK[A][B]", "multiple modifications on the last residue"),
            ("PEP[Carbamidomethyl(C)]TIDE", "parentheses"),
        ],
    )
    def test_refused_features(self, proforma: str, feature: str):
        with pytest.raises(PeptacularError, match=feature):
            pt.parse(proforma).to_ip2()

    def test_modification_multiplier_refused(self):
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"Acetyl": 2})
        with pytest.raises(PeptacularError, match="multipliers"):
            annot.to_ip2()

    def test_does_not_modify_input(self):
        annot = pt.parse("<[Carbamidomethyl]@C>PECTIDE")
        annot.to_ip2()
        assert annot.serialize() == "<[Carbamidomethyl]@C>PECTIDE"
