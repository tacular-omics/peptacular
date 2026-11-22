import pytest

import peptacular as pt
from peptacular import (
    PeptacularError,
    ProFormaFormatError,
    convert_casanovo_sequence,
    convert_diann_sequence,
    convert_ip2_sequence,
    convert_maxquant_sequence,
    convert_skyline_sequence,
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
            ("PEPTIDE/[Na:z+1]", "charge adducts"),
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


class TestToCasanovo:
    """Tests for ProFormaAnnotation.to_casanovo."""

    @pytest.mark.parametrize(
        ("proforma", "expected"),
        [
            ("PEPTIDE", "PEPTIDE"),
            ("PEPM[+15.995]IDE", "PEPM+15.995IDE"),
            ("[+43.006]-PEPTIDE", "+43.006PEPTIDE"),
            ("[+42.011]-PEPM[+15.995]IDEC[+57.021]", "+42.011PEPM+15.995IDEC+57.021"),
            ("<[+57.021]@C>PECTIDE", "PEC+57.021TIDE"),
        ],
    )
    def test_output(self, proforma: str, expected: str):
        assert pt.parse(proforma).to_casanovo() == expected

    @pytest.mark.parametrize(
        "casanovo",
        [
            "PEPTIDE",
            "+43.006PEPTIDE",
            "PEPM+15.995IDE",
            "+42.011PEPM+15.995IDEC+57.021",
            "+42.011PEPM-17.027IDEC+57.021",
        ],
    )
    def test_round_trip(self, casanovo: str):
        original = pt.ProFormaAnnotation.from_casanovo(casanovo)
        written = original.to_casanovo()
        again = pt.ProFormaAnnotation.from_casanovo(written)
        assert again == original
        assert again.mass() == pytest.approx(original.mass())

    @pytest.mark.parametrize(
        ("proforma", "feature"),
        [
            ("PEP[Phospho]TIDE", "signed delta mass"),
            ("PEP[UNIMOD:21]TIDE", "signed delta mass"),
            ("PEP[Formula:HO3P]TIDE", "signed delta mass"),
            ("PEP[Glycan:HexNAc]TIDE", "signed delta mass"),
            ("PEP[15.995]TIDE", "signed delta mass"),
            ("[Acetyl]-PEPTIDE", "signed delta mass"),
            ("PEPM[+15.995][+1.0]IDE", "multiple modifications"),
            ("[+42.011][+1.0]-PEPTIDE", "multiple N-terminal modifications"),
            ("PEPTIDE-[-17.027]", "C-terminal"),
            ("<13C>PEPTIDE", "isotope"),
            ("{Hex}PEPTIDE", "labile"),
            ("[+79.966]?PEPTIDE", "unknown"),
            ("PEP(TI)[+79.966]DE", "interval"),
            ("PEPTIDE/2", "charge"),
            ("PEPTIDE/[Na:z+1]", "charge adducts"),
        ],
    )
    def test_refused_features(self, proforma: str, feature: str):
        with pytest.raises(PeptacularError, match=feature):
            pt.parse(proforma).to_casanovo()

    def test_modification_multiplier_refused(self):
        annot_nterm = pt.ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"+42.011": 2})
        with pytest.raises(PeptacularError, match="multipliers"):
            annot_nterm.to_casanovo()

        annot_internal = pt.ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"+15.995": 2}})
        with pytest.raises(PeptacularError, match="multipliers"):
            annot_internal.to_casanovo()

    def test_does_not_modify_input(self):
        annot = pt.parse("<[+57.021]@C>PECTIDE")
        annot.to_casanovo()
        assert annot.serialize() == "<[+57.021]@C>PECTIDE"


class TestToDiann:
    """Tests for ProFormaAnnotation.to_diann."""

    @pytest.mark.parametrize(
        ("proforma", "expected"),
        [
            ("PEPTIDE", "_PEPTIDE_"),
            ("PEP[Phospho]TIDE", "_PEP[Phospho]TIDE_"),
            ("PEP[+79.9663]TIDE", "_PEP[+79.9663]TIDE_"),
            ("[Acetyl]-PEPTIDE", "_[Acetyl]PEPTIDE_"),
            ("PEPTIDE-[Amidated]", "_PEPTIDE_[Amidated]"),
            ("[Acetyl]-PEP[Phospho]TIDE-[Amidated]", "_[Acetyl]PEP[Phospho]TIDE_[Amidated]"),
            ("<[Carbamidomethyl]@C>PECTIDE", "_PEC[Carbamidomethyl]TIDE_"),
            ("PEP[Phospho][Oxidation]TIDE", "_PEP[Phospho][Oxidation]TIDE_"),
            ("[Acetyl]-P[Phospho]-[Amidated]", "_[Acetyl]P[Phospho]_[Amidated]"),
        ],
    )
    def test_output(self, proforma: str, expected: str):
        assert pt.parse(proforma).to_diann() == expected

    @pytest.mark.parametrize(
        "diann",
        [
            "_PEPTIDE_",
            "_PEP[Phospho]TIDE_",
            "_PEP[+79.9663]TIDE_",
            "_[Acetyl]PEPTIDE_",
            "_PEPTIDE_[Amidated]",
            "_[Acetyl]PEP[Phospho]TIDE_[Amidated]",
            "_[Acetyl]PEP[Phospho]TIDE[Oxidation]_[Amidated]",
            "_PEP[Phospho][Oxidation]TIDE_",
            "_[Acetyl]P[Phospho]_[Amidated]",
        ],
    )
    def test_round_trip(self, diann: str):
        original = pt.ProFormaAnnotation.from_diann(diann)
        written = original.to_diann()
        again = pt.ProFormaAnnotation.from_diann(written)
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
            ("PEPTIDE/[Na:z+1]", "charge adducts"),
        ],
    )
    def test_refused_features(self, proforma: str, feature: str):
        with pytest.raises(PeptacularError, match=feature):
            pt.parse(proforma).to_diann()

    def test_modification_multiplier_refused(self):
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"Acetyl": 2})
        with pytest.raises(PeptacularError, match="multipliers"):
            annot.to_diann()

    @pytest.mark.parametrize(
        ("proforma", "unwritable"),
        [
            ("[Acetyl][Phospho]-PEPTIDE", "_[Acetyl][Phospho]PEPTIDE_"),
            ("[Acetyl][Methyl]-PEP[Phospho]TIDE", "_[Acetyl][Methyl]PEP[Phospho]TIDE_"),
        ],
    )
    def test_two_nterm_mods_refused(self, proforma: str, unwritable: str):
        """`from_diann` dashes only the first leading group, so the second one cannot be read back."""

        with pytest.raises(PeptacularError, match="more than one N-terminal modification"):
            pt.parse(proforma).to_diann()

        with pytest.raises(ProFormaFormatError):
            pt.ProFormaAnnotation.from_diann(unwritable)

    @pytest.mark.parametrize(
        ("proforma", "unwritable"),
        [
            ("PEPTIDE-[Amidated][Phospho]", "_PEPTIDE_[Amidated][Phospho]"),
            ("[Acetyl]-P[Phospho]-[Amidated][Phospho]", "_[Acetyl]P[Phospho]_[Amidated][Phospho]"),
        ],
    )
    def test_two_cterm_mods_refused(self, proforma: str, unwritable: str):
        """`from_diann` rewrites one trailing group, so with two the inner underscore reaches the parser."""

        with pytest.raises(PeptacularError, match="more than one C-terminal modification"):
            pt.parse(proforma).to_diann()

        with pytest.raises(ProFormaFormatError):
            pt.ProFormaAnnotation.from_diann(unwritable)

    def test_does_not_modify_input(self):
        annot = pt.parse("<[Carbamidomethyl]@C>PECTIDE")
        annot.to_diann()
        assert annot.serialize() == "<[Carbamidomethyl]@C>PECTIDE"


MAXQUANT_PAIRS = [
    ("PEPTIDE", "_PEPTIDE_"),
    ("PEPM[Oxidation (M)]TIDE", "_PEPM(Oxidation (M))TIDE_"),
    ("[ac]-SM[ox]K", "_(ac)SM(ox)K_"),
    ("PEPTIDE-[Amidated]", "_PEPTIDE(Amidated)_"),
    ("PEP[Phospho]TIDE-[Amidated]", "_PEP(Phospho)TIDE(Amidated)_"),
]


class TestToMaxquant:
    """Tests for ProFormaAnnotation.to_maxquant."""

    @pytest.mark.parametrize(("proforma", "expected"), MAXQUANT_PAIRS)
    def test_output(self, proforma: str, expected: str):
        assert pt.parse(proforma).to_maxquant() == expected

    def test_static_mod_condensed(self):
        assert pt.parse("<[Carbamidomethyl]@C>PECTIDE").to_maxquant() == "_PEC(Carbamidomethyl)TIDE_"

    @pytest.mark.parametrize(
        ("proforma", "feature"),
        [
            ("<13C>PEPTIDE", "isotope"),
            ("{Hex}PEPTIDE", "labile"),
            ("[Phospho]?PEPTIDE", "unknown"),
            ("PEP(TI)[Phospho]DE", "interval"),
            ("PEPTIDE/2", "charge"),
            ("PEPTIDE/[Na:z+1]", "charge adducts"),
        ],
    )
    def test_refused_features(self, proforma: str, feature: str):
        with pytest.raises(PeptacularError, match=feature):
            pt.parse(proforma).to_maxquant()

    def test_modification_multiplier_refused(self):
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"Acetyl": 2})
        with pytest.raises(PeptacularError, match="multipliers"):
            annot.to_maxquant()

    @pytest.mark.parametrize(
        ("proforma", "where"),
        [
            ("[Acetyl][Phospho]-PEPTIDE", "at the N-terminus"),
            ("PEPTIDE-[Amidated][Phospho]", "at the C-terminus"),
            ("PEP[Phospho][Oxidation]TIDE", "on one residue"),
        ],
    )
    def test_two_mods_at_one_position_refused(self, proforma: str, where: str):
        with pytest.raises(PeptacularError, match=f"more than one modification {where}"):
            pt.parse(proforma).to_maxquant()

    @pytest.mark.parametrize("proforma", ["PEPTIDEK[TMT6plex]", "PEPTIDE[Phospho]-[Amidated]"])
    def test_last_residue_mod_refused(self, proforma: str):
        with pytest.raises(PeptacularError, match="last residue"):
            pt.parse(proforma).to_maxquant()

    def test_last_residue_mod_is_read_as_cterm(self):
        assert pt.ProFormaAnnotation.from_maxquant("_PEPTIDEK(TMT6plex)_").serialize() == "PEPTIDEK-[TMT6plex]"

    def test_unbalanced_name_refused(self):
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"a)": 1}})
        with pytest.raises(PeptacularError, match="unbalanced"):
            annot.to_maxquant()

    def test_does_not_modify_input(self):
        annot = pt.parse("<[Carbamidomethyl]@C>PECTIDE")
        annot.to_maxquant()
        assert annot.serialize() == "<[Carbamidomethyl]@C>PECTIDE"


class TestFromMaxquant:
    """Tests for ProFormaAnnotation.from_maxquant and convert_maxquant_sequence."""

    @pytest.mark.parametrize(("proforma", "mq"), MAXQUANT_PAIRS)
    def test_round_trip(self, proforma: str, mq: str):
        annot = pt.ProFormaAnnotation.from_maxquant(mq)
        assert annot == pt.parse(proforma)
        assert annot.to_maxquant() == mq

    @pytest.mark.parametrize(
        ("mq", "expected"),
        [
            ("_ALTDM(Oxidation (M))PQM(Oxidation (M))R_", "ALTDM[Oxidation (M)]PQM[Oxidation (M)]R"),
            ("_P(Phospho)EPTIDE_", "P[Phospho]EPTIDE"),
            ("_PEPTIDE(Amidated)_", "PEPTIDE-[Amidated]"),
            ("_(ac)SM(ox)K_", "[ac]-SM[ox]K"),
            ("PEPTIDE", "PEPTIDE"),
        ],
    )
    def test_reads(self, mq: str, expected: str):
        assert pt.ProFormaAnnotation.from_maxquant(mq).serialize() == expected
        assert convert_maxquant_sequence(mq) == expected

    def test_first_residue_mod_is_not_nterminal(self):
        annot = pt.ProFormaAnnotation.from_maxquant("_P(Phospho)EPTIDE_")
        assert not annot.has_nterm_mods

    @pytest.mark.parametrize("bad", ["_PEP(Phospho_", "_PEP)TIDE_", "_PEP(Phospho))TIDE_"])
    def test_unbalanced_raises(self, bad: str):
        with pytest.raises(ProFormaFormatError):
            pt.ProFormaAnnotation.from_maxquant(bad)

    def test_convert_list(self):
        seqs = ["_PEPTIDE_", "_(ac)SM(ox)K_", "_PEPM(Oxidation (M))TIDE_"]
        expected = ["PEPTIDE", "[ac]-SM[ox]K", "PEPM[Oxidation (M)]TIDE"]
        assert convert_maxquant_sequence(seqs, n_workers=2) == expected
        assert convert_maxquant_sequence(seqs, method="thread", n_workers=2) == expected
        assert convert_maxquant_sequence([]) == []


class TestBracketedModificationNamesAtATerminus:
    """A name holding a square bracket breaks the readers' terminal detection.

    Both readers delimit a terminal group with a bracket pattern that stops at the first
    ``]``, so a terminal name such as ``Formula:[13C2]H3`` was written into a string that
    either does not parse or parses as something else. The writers refuse it instead. An
    internal name is read back unchanged and must keep working.
    """

    ISOTOPE_NAME = "Formula:[13C2]H3"

    @pytest.mark.parametrize(
        "proforma",
        [
            "[Formula:[13C2]H3]-PEPTIDE",
            "PEPTIDE-[Formula:[13C2]H3]",
        ],
    )
    def test_both_writers_refuse_a_bracketed_terminal_name(self, proforma: str):
        for writer in ("to_ip2", "to_diann"):
            with pytest.raises(PeptacularError, match="square brackets"):
                getattr(pt.parse(proforma), writer)()

    def test_ip2_refuses_a_cterm_modification_after_a_bracketed_last_residue_name(self):
        # Written out it was `EEN(Formula:[13C2]H3)(UNIMOD:35)`, which the reader gave back
        # as `EEN[Formula:[13C2]H3][UNIMOD:35]`: the C-terminal modification silently moved
        # onto the residue.
        with pytest.raises(PeptacularError, match="square brackets"):
            pt.parse("EEN[Formula:[13C2]H3]-[UNIMOD:35]").to_ip2()

    @pytest.mark.parametrize(
        "proforma",
        [
            "PEP[Formula:[13C2]H3]TIDE",
            "PEPTIDEK[Formula:[13C2]H3]",
            "[UNIMOD:1]-EEN[Formula:[13C2]H3]",
            "PEPTIDE[Phospho]-[Amidated]",
        ],
    )
    def test_ip2_still_round_trips_what_it_can_carry(self, proforma: str):
        written = pt.parse(proforma).to_ip2()
        assert pt.ProFormaAnnotation.from_ip2_sequence(written).serialize() == proforma

    @pytest.mark.parametrize(
        "proforma",
        [
            "PEP[Formula:[13C2]H3]TIDE",
            "PEPTIDEK[Formula:[13C2]H3]",
            "[UNIMOD:1]-EEN[Formula:[13C2]H3]",
            "EEN[Formula:[13C2]H3]-[UNIMOD:35]",
            "PEPTIDE[Phospho]-[Amidated]",
            "PEPTIDE-[UNIMOD:35]",
        ],
    )
    def test_diann_still_round_trips_what_it_can_carry(self, proforma: str):
        written = pt.parse(proforma).to_diann()
        assert pt.ProFormaAnnotation.from_diann(written).serialize() == proforma


SKYLINE_PAIRS = [
    ("PEPTIDE", "PEPTIDE"),
    ("PEPT[+79.9663]IDE", "PEPT[+79.9663]IDE"),
    ("PEPT[-18.0106]IDE", "PEPT[-18.0106]IDE"),
    ("PEPC[Carbamidomethyl]IDE", "PEPC[Carbamidomethyl]IDE"),
]


class TestToSkyline:
    """Tests for ProFormaAnnotation.to_skyline."""

    @pytest.mark.parametrize(("proforma", "expected"), SKYLINE_PAIRS)
    def test_output(self, proforma: str, expected: str):
        assert pt.parse(proforma).to_skyline() == expected

    def test_static_mod_condensed(self):
        assert pt.parse("<[Carbamidomethyl]@C>PECTIDE").to_skyline() == "PEC[Carbamidomethyl]TIDE"

    def test_last_residue_mod_written(self):
        assert pt.parse("GSMAFVPK[+8.0]").to_skyline() == "GSMAFVPK[+8.0]"

    @pytest.mark.parametrize(
        ("proforma", "feature"),
        [
            ("<13C>PEPTIDE", "isotope"),
            ("{Hex}PEPTIDE", "labile"),
            ("[Phospho]?PEPTIDE", "unknown"),
            ("PEP(TI)[Phospho]DE", "interval"),
            ("PEPTIDE/2", "charge"),
            ("PEPTIDE/[Na:z+1]", "charge adducts"),
        ],
    )
    def test_refused_features(self, proforma: str, feature: str):
        with pytest.raises(PeptacularError, match=feature):
            pt.parse(proforma).to_skyline()

    def test_modification_multiplier_refused(self):
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"Phospho": 2}})
        with pytest.raises(PeptacularError, match="multipliers"):
            annot.to_skyline()

    @pytest.mark.parametrize("proforma", ["[Acetyl]-PEPTIDE", "PEPTIDE-[Amidated]", "[Acetyl]-PEPTIDE-[Amidated]"])
    def test_terminal_mod_refused(self, proforma: str):
        with pytest.raises(PeptacularError, match="terminal modification"):
            pt.parse(proforma).to_skyline()

    def test_two_mods_on_one_residue_refused(self):
        with pytest.raises(PeptacularError, match="more than one modification on one residue"):
            pt.parse("PEP[+1][+2]TIDE").to_skyline()

    def test_unsigned_number_refused(self):
        with pytest.raises(PeptacularError, match="unsigned number"):
            pt.parse("PEPT[15.9949]IDE").to_skyline()

    def test_unbalanced_name_refused(self):
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"a]": 1}})
        with pytest.raises(PeptacularError, match="unbalanced"):
            annot.to_skyline()

    def test_does_not_modify_input(self):
        annot = pt.parse("<[Carbamidomethyl]@C>PECTIDE")
        annot.to_skyline()
        assert annot.serialize() == "<[Carbamidomethyl]@C>PECTIDE"


class TestFromSkyline:
    """Tests for ProFormaAnnotation.from_skyline and convert_skyline_sequence."""

    @pytest.mark.parametrize(("proforma", "sky"), SKYLINE_PAIRS)
    def test_round_trip(self, proforma: str, sky: str):
        annot = pt.ProFormaAnnotation.from_skyline(sky)
        assert annot == pt.parse(proforma)
        assert annot.to_skyline() == sky

    @pytest.mark.parametrize(
        ("sky", "expected"),
        [
            ("GVDC[+57]QEVSQEK", "GVDC[+57]QEVSQEK"),
            ("GSMAFVPK(unimod:259)", "GSMAFVPK[UNIMOD:259]"),
            ("GSMAFVPK(UniMod:259)", "GSMAFVPK[UNIMOD:259]"),
            ("GSMAFVPK[+8.0]", "GSMAFVPK[+8.0]"),
            ("PEPT[Phospho (ST)]IDER", "PEPT[Phospho (ST)]IDER"),
            ("AAM[15.9949]Q[-0.98]R", "AAM[+15.9949]Q[-0.98]R"),
        ],
    )
    def test_reads(self, sky: str, expected: str):
        assert pt.ProFormaAnnotation.from_skyline(sky).serialize() == expected
        assert convert_skyline_sequence(sky) == expected

    def test_masses(self):
        assert pt.ProFormaAnnotation.from_skyline("GSMAFVPK(unimod:259)").mass() == pytest.approx(843.4404, abs=1e-4)
        assert pt.ProFormaAnnotation.from_skyline("AAM[15.9949]Q[-0.98]R").mass() == pytest.approx(590.2999, abs=1e-4)

    def test_unsigned_number_gets_plus(self):
        annot = pt.ProFormaAnnotation.from_skyline("PEPT[15.9949]IDE")
        assert annot.serialize() == "PEPT[+15.9949]IDE"
        assert annot.mass() == pytest.approx(pt.parse("PEPTIDE").mass() + 15.9949, abs=1e-4)

    def test_balanced_square_brackets(self):
        # peptacular's own choice, not Skyline's behaviour. Skyline cannot emit this: a name
        # holding `]` is written in round or curly brackets instead
        # (pwiz_tools/Skyline/Model/ModifiedSequence.cs:324-336), and Skyline's reader closes
        # at the first `]` (Model/ModificationMatcher.cs:294). Reading the balanced pair is
        # strictly more permissive, so it costs nothing, but no Skyline string exercises it.
        assert pt.ProFormaAnnotation.from_skyline("PEPT[a[b]c]IDE").serialize() == "PEPT[a[b]c]IDE"

    def test_modification_is_never_terminal(self):
        annot = pt.ProFormaAnnotation.from_skyline("P[+42]EPTIDEK[+8]")
        assert not annot.has_nterm_mods
        assert not annot.has_cterm_mods

    # --- strings taken from Skyline's own repository, ProteoWizard/pwiz at 408a38ca --------
    # Each label is the file the string was read from, verified 2026-11-19. Base path:
    # https://github.com/ProteoWizard/pwiz/blob/408a38ca6e2042ba886a2cd495eaa989be4436f8/pwiz_tools/Skyline/
    @pytest.mark.parametrize(
        ("skyline", "proforma", "provenance"),
        [
            ("GVDC[+57]QEVSQEK", "GVDC[+57]QEVSQEK", "TestTutorial/TutorialAuditLogs/en/TestOrbiPrmTutorial.log:69, real Skyline output"),
            ("GVDC[+57.021464]QEVSQEK", "GVDC[+57.021464]QEVSQEK", "the same log at :198, the full-precision form of the same peptide"),
            ("C[+57.0]PEPTIDER[+10.0]", "C[+57.0]PEPTIDER[+10.0]", "Test/ModifiedSequenceTest.cs:72, modified_sequence of a heavy precursor"),
            ("A[+56.0]GLQFPVGR", "A[+56.0]GLQFPVGR", "Test/ModifiedSequenceTest.cs:191, an N-terminal PropionylNT written onto residue 0"),
            (
                "AFC[UniMod:4]AVPWQGTM[oXi]T(UniMod:21)LSK",
                "AFC[UniMod:4]AVPWQGTM[oXi]T[UNIMOD:21]LSK",
                "Test/ModificationMatcherTest.cs:392 STR_CYS_OXI_PHOS, square and round UniMod in one string",
            ),
        ],
    )
    def test_reads_real_skyline_strings(self, skyline: str, proforma: str, provenance: str):
        assert pt.ProFormaAnnotation.from_skyline(skyline).serialize() == proforma, provenance

    def test_reads_the_two_forms_the_skyline_tutorial_lists(self):
        # Documentation/Tutorials/ImportingIntegrationBoundaries/en/index.html: "peptide
        # modified sequences can be expressed in many different ways, so GSMAFVPK(UniMod:259),
        # GSMAFVPK[+8.0], and GSMAFVPK[Label (K) 6] are all acceptable." The first two are the
        # ones peptacular resolves to a mass.
        by_mass = pt.ProFormaAnnotation.from_skyline("GSMAFVPK[+8.0]")
        by_unimod = pt.ProFormaAnnotation.from_skyline("GSMAFVPK(UniMod:259)")
        assert by_mass.serialize() == "GSMAFVPK[+8.0]"
        assert by_unimod.serialize() == "GSMAFVPK[UNIMOD:259]"
        # UNIMOD:259 is Label:13C(6)15N(2), 8.014199 Da, so the tutorial's two spellings are
        # the same peptide to within the rounding of its own +8.0.
        assert by_unimod.mass() == pytest.approx(by_mass.mass(), abs=0.015)

    # --- fix: a bracket written with padding ----------------------------------------------
    @pytest.mark.parametrize("skyline", ["PEPT[ 15.99]IDE", "PEPT[15.99 ]IDE", "PEPT[  +15.99  ]IDE"])
    def test_bracket_text_is_stripped_before_the_number_test(self, skyline: str):
        # Unstripped, the bare-number test missed these and the mass became a modification
        # named "15.99", which mass() then refused.
        annot = pt.ProFormaAnnotation.from_skyline(skyline)
        assert annot.serialize() == "PEPT[+15.99]IDE"
        assert annot.mass() == pytest.approx(pt.parse("PEPTIDE").mass() + 15.99, abs=1e-4)

    def test_a_bracket_of_only_whitespace_is_empty(self):
        with pytest.raises(ProFormaFormatError, match="Empty modification"):
            pt.ProFormaAnnotation.from_skyline("PEPT[ ]IDE")

    # --- fix: a comma decimal is Skyline's culture-dependent mass form ---------------------
    @pytest.mark.parametrize("skyline", ["PEPM[15,9949]K", "PEPM[+15,9949]K", "PEPM[-0,98]K", "PEPM[1,000]K"])
    def test_comma_decimal_is_refused(self, skyline: str):
        # Model/ModifiedSequence.cs:308 writes the mass with CultureInfo.CurrentCulture and
        # Model/MassModification.cs:108 accepts `[.,]` as the separator, so the value depends
        # on the culture that wrote it: "1,000" is 1000 in one and 1.0 in another.
        with pytest.raises(ProFormaFormatError, match="Comma decimal"):
            pt.ProFormaAnnotation.from_skyline(skyline)

    def test_a_name_holding_a_comma_is_not_refused(self):
        # The refusal matches only a bare number, so a modification name keeps working.
        assert pt.ProFormaAnnotation.from_skyline("PEPT[Label:2H(4),13C(1)]IDE").serialize() == "PEPT[Label:2H(4),13C(1)]IDE"

    # --- fix: the guard against a modification before the first residue -------------------
    @pytest.mark.parametrize("skyline", ["[+42]-PEPTIDE", "[+1]?PEPTIDE", "[+42]PEPTIDE", "(unimod:35)PEPTIDE"])
    def test_modification_before_the_first_residue_is_refused(self, skyline: str):
        # Skyline's reader does accept this shape and folds it onto residue 0
        # (Model/PeptideGroup.cs:409 NormalizeNTerminalMod), but Skyline never writes it, and
        # reading it would put a modification where peptacular cannot write it back. Without
        # this guard the first two parse as ProForma terminal and unknown-position
        # modifications instead of being refused.
        with pytest.raises(ProFormaFormatError):
            pt.ProFormaAnnotation.from_skyline(skyline)

    @pytest.mark.parametrize(
        "bad",
        ["PEPT[+1IDE", "PEPT+1]IDE", "PEPT[]IDE", "[+1]PEPTIDE", "PEPT(phospho)IDE", "PEPT(unimod:35IDE", "(unimod:35)PEPTIDE"],
    )
    def test_malformed_raises(self, bad: str):
        with pytest.raises(ProFormaFormatError):
            pt.ProFormaAnnotation.from_skyline(bad)

    def test_convert_list(self):
        seqs = ["PEPTIDE", "AAM[15.9949]Q[-0.98]R", "GSMAFVPK(unimod:259)"]
        expected = ["PEPTIDE", "AAM[+15.9949]Q[-0.98]R", "GSMAFVPK[UNIMOD:259]"]
        assert convert_skyline_sequence(seqs, n_workers=2) == expected
        assert convert_skyline_sequence(seqs, method="thread", n_workers=2) == expected
        assert convert_skyline_sequence([]) == []


class TestMSFraggerConverter:
    """Tests for MSFragger format conversion (to_msfragger, from_msfragger)."""

    # Issue #116 listed pair 3 as "N-term(42.0106), 4M(15.9949)". Philosopher sorts the entry
    # strings (`sort.Strings(assL)`, lib/rep/psm.go:564), so the residue entry comes first:
    # "4M(15.9949), N-term(42.0106)", which is also what the real rows in
    # pFindStudio/pDeep3 `psm.tsv` show ("7M(15.9949), N-term(42.0106)").
    PAIRS = [
        ("PEPTIDE", ("PEPTIDE", "")),
        ("PEPM[+15.9949]IDE", ("PEPMIDE", "4M(15.9949)")),
        ("[+42.0106]-PEPM[+15.9949]IDE", ("PEPMIDE", "4M(15.9949), N-term(42.0106)")),
        ("PEPTIDE-[-0.9840]", ("PEPTIDE", "C-term(-0.9840)")),
        ("PEPC[+57.0215]M[+15.9949]IDE", ("PEPCMIDE", "4C(57.0215), 5M(15.9949)")),
    ]

    @pytest.mark.parametrize(("proforma", "expected_pair"), PAIRS)
    def test_to_msfragger_target_pairs(self, proforma: str, expected_pair: tuple[str, str]):
        assert pt.parse(proforma).to_msfragger() == expected_pair
        assert pt.to_msfragger(proforma) == expected_pair

    def test_from_msfragger_target_pairs(self):
        # 1. PEPTIDE -> ('PEPTIDE', '') -> PEPTIDE
        assert pt.ProFormaAnnotation.from_msfragger("PEPTIDE", "") == pt.parse("PEPTIDE")
        assert pt.from_msfragger(("PEPTIDE", "")) == "PEPTIDE"

        # 2. PEPM[+15.9949]IDE -> ('PEPMIDE', '4M(15.9949)') -> PEPM[+15.9949]IDE
        assert pt.ProFormaAnnotation.from_msfragger("PEPMIDE", "4M(15.9949)") == pt.parse("PEPM[+15.9949]IDE")
        assert pt.from_msfragger(("PEPMIDE", "4M(15.9949)")) == "PEPM[+15.9949]IDE"

        # 3. [+42.0106]-PEPM[+15.9949]IDE -> ('PEPMIDE', 'N-term(42.0106), 4M(15.9949)') -> [+42.0106]-PEPM[+15.9949]IDE
        assert pt.ProFormaAnnotation.from_msfragger("PEPMIDE", "N-term(42.0106), 4M(15.9949)") == pt.parse("[+42.0106]-PEPM[+15.9949]IDE")
        assert pt.from_msfragger(("PEPMIDE", "N-term(42.0106), 4M(15.9949)")) == "[+42.0106]-PEPM[+15.9949]IDE"

        # 4. PEPTIDE-[-0.9840] -> ('PEPTIDE', 'C-term(-0.9840)') -> PEPTIDE-[-0.984]
        # (psm_utils gives PEPTIDE-[-0.984] because -0.9840 parsed as float normalizes trailing zeros)
        assert pt.ProFormaAnnotation.from_msfragger("PEPTIDE", "C-term(-0.9840)") == pt.parse("PEPTIDE-[-0.984]")
        assert pt.from_msfragger(("PEPTIDE", "C-term(-0.9840)")) == "PEPTIDE-[-0.984]"

        # 5. PEPC[+57.0215]M[+15.9949]IDE -> ('PEPCMIDE', '4C(57.0215), 5M(15.9949)') -> PEPC[+57.0215]M[+15.9949]IDE
        assert pt.ProFormaAnnotation.from_msfragger("PEPCMIDE", "4C(57.0215), 5M(15.9949)") == pt.parse("PEPC[+57.0215]M[+15.9949]IDE")
        assert pt.from_msfragger(("PEPCMIDE", "4C(57.0215), 5M(15.9949)")) == "PEPC[+57.0215]M[+15.9949]IDE"

    # --- rows from real FragPipe psm.tsv files and from Philosopher's own test ------------
    # Fetched 2026-11-19. The third column of each tuple is where the row came from.
    REAL_ROWS = [
        (
            "FLLEAGADQEHK",
            "",
            "FLLEAGADQEHK",
            "compomics/psm_utils example_files/fragpipe.psm.tsv, the one row with an empty cell (its Modified Peptide is blank too)",
        ),
        ("SEDCFILDHGK", "4C(57.0214)", "SEDC[+57.0214]FILDHGK", "the same file, a fixed carbamidomethyl in Assigned Modifications"),
        ("REELSNVLAAMR", "11M(15.9949)", "REELSNVLAAM[+15.9949]R", "the same file"),
        (
            "LHISPSNMTNQNTPEYMEK",
            "8M(15.9949), C-term(-0.9840)",
            "LHISPSNM[+15.9949]TNQNTPEYMEK-[-0.984]",
            "the same file, a C-terminal entry",
        ),
        ("ANIAVQR", "N-term(42.0106)", "[+42.0106]-ANIAVQR", "the same file, an N-terminal entry"),
        (
            "CEMSEHPSQTVRK",
            "1C(57.0215), 3M(15.9949)",
            "C[+57.0215]EM[+15.9949]SEHPSQTVRK",
            "pFindStudio/pDeep3 pDeep/tmp/MSFragger-tsv-from-fc/psm.tsv",
        ),
        (
            "THTGERPNACADCGK",
            "10C(57.0215), 13C(57.0215)",
            "THTGERPNAC[+57.0215]ADC[+57.0215]GK",
            "the same pDeep3 file, two residues of one kind",
        ),
        (
            "AASGSGMAQK",
            "7M(15.9949), N-term(42.0106)",
            "[+42.0106]-AASGSGM[+15.9949]AQK",
            "the same pDeep3 file, the order sort.Strings gives",
        ),
    ]

    @pytest.mark.parametrize(("peptide", "modifications", "proforma", "provenance"), REAL_ROWS)
    def test_reads_real_fragpipe_rows(self, peptide: str, modifications: str, proforma: str, provenance: str):
        assert pt.ProFormaAnnotation.from_msfragger(peptide, modifications).serialize() == proforma, provenance

    @pytest.mark.parametrize(("peptide", "modifications", "proforma", "provenance"), REAL_ROWS)
    def test_real_fragpipe_rows_rewrite_byte_identically(self, peptide: str, modifications: str, proforma: str, provenance: str):
        # Every one of these rewrites to exactly the string FragPipe wrote, which needs both
        # the four-decimal format and the string sort.
        annot = pt.ProFormaAnnotation.from_msfragger(peptide, modifications)
        assert annot.to_msfragger() == (peptide, modifications), provenance

    @pytest.mark.parametrize(
        ("peptide", "modifications", "calculated_peptide_mass"),
        [
            ("CEMSEHPSQTVRK", "1C(57.0215), 3M(15.9949)", 1603.7081),
            ("THTGERPNACADCGK", "10C(57.0215), 13C(57.0215)", 1672.7046),
            ("AASGSGMAQK", "7M(15.9949), N-term(42.0106)", 964.4284),
        ],
    )
    def test_mass_matches_the_calculated_peptide_mass_column(self, peptide: str, modifications: str, calculated_peptide_mass: float):
        # The expected value is the row's own "Calculated Peptide Mass" cell in
        # pFindStudio/pDeep3 pDeep/tmp/MSFragger-tsv-from-fc/psm.tsv, so this compares against
        # the tool's own arithmetic rather than against ours.
        annot = pt.ProFormaAnnotation.from_msfragger(peptide, modifications)
        assert annot.neutral_mass() == pytest.approx(calculated_peptide_mass, abs=0.001)

    def test_reads_two_entries_on_one_residue(self):
        # Philosopher's own test expects exactly this, with the same peptide:
        # lib/id/pep_test.go:116 wants {"15M:15.9949:true", "16C:-32.0085:true",
        # "16C:57.0214:false", "N-term:42.0106:true"} for peptide "DDDIAALVVDNGSGMCK", and its
        # ModifiedPeptide is "n[43]DDDIAALVVDNGSGM[147]C[71]K". The `false` is a fixed
        # modification, which is why that column is not only variable ones.
        mods = "15M(15.9949), 16C(-32.0085), 16C(57.0214), N-term(42.0106)"
        annot = pt.ProFormaAnnotation.from_msfragger("DDDIAALVVDNGSGMCK", mods)
        assert annot.serialize() == "[+42.0106]-DDDIAALVVDNGSGM[+15.9949]C[-32.0085][+57.0214]K"
        assert annot.to_msfragger() == ("DDDIAALVVDNGSGMCK", mods)
        assert pt.to_msfragger(annot.serialize()) == ("DDDIAALVVDNGSGMCK", mods)

    def test_to_msfragger_writes_two_modifications_on_one_residue_roundtrip(self):
        # Example from issue #138 and Philosopher lib/id/pep_test.go line 116
        pair = ("PEPCIDE", "4C(-32.0085), 4C(57.0214)")
        proforma = pt.from_msfragger(pair)
        assert proforma == "PEPC[-32.0085][+57.0214]IDE"
        assert pt.to_msfragger(proforma) == pair

    def test_entries_are_sorted_as_strings_not_numerically(self):
        # `sort.Strings(assL)` (lib/rep/psm.go:564) is a plain string sort, so "10P" precedes
        # "4P", which a numeric sort would reverse. The real row is
        # "10P(15.9949),4P(15.9949)" in siwar003/Hard-Matrisome-Digestion-Project data/psm/C3_psm.tsv.
        annot = pt.ProFormaAnnotation.from_msfragger("PPPPPPPPPPPP", "10P(15.9949),4P(15.9949)")
        assert annot.to_msfragger() == ("PPPPPPPPPPPP", "10P(15.9949), 4P(15.9949)")
        # And a residue entry sorts before "N-term" because "4" < "N".
        assert pt.parse("[+42.0106]-PEPM[+15.9949]IDE").to_msfragger() == ("PEPMIDE", "4M(15.9949), N-term(42.0106)")

    def test_mass_is_always_written_with_four_decimals(self):
        # `fmt.Sprintf("%s%s(%.4f)", ...)` (lib/rep/rep.go:596).
        assert pt.parse("PEPM[+16]IDE").to_msfragger() == ("PEPMIDE", "4M(16.0000)")
        assert pt.parse("PEPM[+15.994915]IDE").to_msfragger() == ("PEPMIDE", "4M(15.9949)")
        assert pt.parse("PEPM[-0.98]IDE").to_msfragger() == ("PEPMIDE", "4M(-0.9800)")

    @pytest.mark.parametrize(
        "peptide",
        ["PEPM[+15.9949]IDE", "PEPTIDE/2", " PEPMIDE ", "PEP TIDE", "PEPT1DE", "pepmide", "PEPTIDE*", ""],
    )
    def test_peptide_must_be_a_bare_residue_sequence(self, peptide: str):
        # The FragPipe guide calls Peptide the "'stripped' peptide sequence". Reading anything
        # else doubled a modification already in the sequence, or carried a charge through.
        with pytest.raises(PeptacularError, match="bare upper-case residue sequence"):
            pt.ProFormaAnnotation.from_msfragger(peptide, "")

    def test_non_ascii_digits_are_not_positions(self):
        with pytest.raises(PeptacularError, match="Invalid MSFragger modification entry"):
            pt.ProFormaAnnotation.from_msfragger("PEPMIDE", "\uff17E(1.0000)")

    @pytest.mark.parametrize("empty", [None, float("nan")])
    def test_none_and_nan_mean_an_empty_cell(self, empty):
        # An empty Assigned Modifications cell is an unmodified peptide: the one such row in
        # compomics/psm_utils example_files/fragpipe.psm.tsv (FLLEAGADQEHK) has a blank
        # Modified Peptide, which the FragPipe guide defines as "blank if peptide is
        # unmodified". pandas hands an empty cell over as NaN.
        assert pt.ProFormaAnnotation.from_msfragger("PEPTIDE", empty) == pt.parse("PEPTIDE")

    @pytest.mark.parametrize("bad", [5, b"4M(15.9949)", ["4M(15.9949)"], {"4M": 15.9949}])
    def test_other_non_string_modifications_are_refused(self, bad):
        with pytest.raises(PeptacularError, match="Assigned Modifications must be a string"):
            pt.ProFormaAnnotation.from_msfragger("PEPMIDE", bad)

    def test_a_tuple_of_pairs_is_a_batch(self):
        # `Sequence[tuple[str, str]]` admits a tuple, so a tuple of pairs matched the batch
        # overload and then took the single path, raising AttributeError on `item[0].strip()`.
        pairs = (("PEPTIDE", ""), ("PEPMIDE", "4M(15.9949)"))
        assert pt.from_msfragger(pairs) == ["PEPTIDE", "PEPM[+15.9949]IDE"]
        assert pt.from_msfragger(pairs + (("ANIAVQR", "N-term(42.0106)"),))[2] == "[+42.0106]-ANIAVQR"
        # A pair of strings is still one pair.
        assert pt.from_msfragger(("PEPMIDE", "4M(15.9949)")) == "PEPM[+15.9949]IDE"

    def test_roundtrip_equivalence(self):
        for proforma, (pep, mods) in self.PAIRS:
            back = pt.ProFormaAnnotation.from_msfragger(pep, mods)
            if proforma == "PEPTIDE-[-0.9840]":
                # The mass is read as a float, so the trailing zero of the ProForma text is
                # gone from the annotation. The written string is unaffected now that the
                # writer formats with %.4f as Philosopher does, so the pair round trips
                # byte-identically where it used to come back as "C-term(-0.984)".
                assert back == pt.parse("PEPTIDE-[-0.984]")
            else:
                assert back == pt.parse(proforma)
            assert back.to_msfragger() == (pep, mods)

    @pytest.mark.parametrize(
        "annot",
        [
            pt.ProFormaAnnotation("PEPTIDE", isotope_mods={"13C": 6}),
            pt.ProFormaAnnotation("PEPTIDE", labile_mods={"Hex": 1}),
            pt.ProFormaAnnotation("PEPTIDE", unknown_mods={"+15.9949": 1}),
            pt.parse("(PEP)[+15.9949]TIDE"),
            pt.parse("PEPTIDE/2"),
            pt.parse("PEPTIDE/+2"),
            pt.ProFormaAnnotation("PEPTIDE", nterm_mods={"+42.0106": 2}),
            pt.ProFormaAnnotation("PEPTIDE", cterm_mods={"-0.9840": 2}),
            pt.ProFormaAnnotation("PEPTIDE", internal_mods={1: {"+15.9949": 2}}),
            pt.parse("PEPM[Oxidation]IDE"),
            pt.parse("PEPM[147]IDE"),
            pt.parse("[Acetyl]-PEPTIDE"),
            pt.parse("PEPTIDE-[Amidated]"),
        ],
    )
    def test_to_msfragger_refusals(self, annot: pt.ProFormaAnnotation):
        with pytest.raises(PeptacularError):
            annot.to_msfragger()

    def test_from_msfragger_refuses_unparseable_entries(self):
        bad_entries = [
            "foo",
            "4(15.9949)",
            "4T",
            "4T()",
            "4T(Oxidation)",
            "N-term()",
            "N-term(Oxidation)",
            "C-term()",
            "C-term(Oxidation)",
            "4T(15.9949), ",
            ", 4T(15.9949)",
        ]
        for bad in bad_entries:
            with pytest.raises(PeptacularError, match="Invalid MSFragger modification entry"):
                pt.ProFormaAnnotation.from_msfragger("PEPTIDE", bad)

    def test_from_msfragger_refuses_position_out_of_range(self):
        with pytest.raises(PeptacularError, match="out of range"):
            pt.ProFormaAnnotation.from_msfragger("PEPTIDE", "0P(15.9949)")
        with pytest.raises(PeptacularError, match="out of range"):
            pt.ProFormaAnnotation.from_msfragger("PEPTIDE", "8E(15.9949)")
        with pytest.raises(PeptacularError, match="out of range"):
            pt.ProFormaAnnotation.from_msfragger("PEPTIDE", "100P(15.9949)")

    def test_from_msfragger_refuses_residue_mismatch_naming_both_residues(self):
        with pytest.raises(PeptacularError) as exc_info:
            pt.ProFormaAnnotation.from_msfragger("PEPTIDE", "4Q(15.9949)")
        msg = str(exc_info.value)
        assert "Q" in msg and "T" in msg

        with pytest.raises(PeptacularError) as exc_info:
            pt.ProFormaAnnotation.from_msfragger("PEPTIDE", "1E(15.9949)")
        msg = str(exc_info.value)
        assert "E" in msg and "P" in msg

    def test_wrappers_batch_and_parallel(self):
        seqs = ["PEPTIDE", "PEPM[+15.9949]IDE"]
        expected_batch = [("PEPTIDE", ""), ("PEPMIDE", "4M(15.9949)")]

        assert pt.to_msfragger(seqs) == expected_batch
        assert pt.to_msfragger(seqs, n_workers=1, chunksize=1, method="sequential") == expected_batch

        assert pt.from_msfragger(expected_batch) == ["PEPTIDE", "PEPM[+15.9949]IDE"]
        assert pt.from_msfragger(expected_batch, n_workers=1, chunksize=1, method="sequential") == ["PEPTIDE", "PEPM[+15.9949]IDE"]

    def test_wrappers_invalid_input(self):
        with pytest.raises(PeptacularError):
            pt.from_msfragger("not_a_tuple")  # type: ignore[arg-type]
        with pytest.raises(PeptacularError):
            pt.from_msfragger(("only_one",))  # type: ignore[arg-type]
        with pytest.raises(PeptacularError):
            pt.from_msfragger([("only_one",)])  # type: ignore[list-item]
