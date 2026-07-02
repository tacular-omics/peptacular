import unittest
import warnings

import peptacular as pt
from peptacular.regex_utils import get_regex_match_indices

PROTEIN = "MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTE"


class TestRegexMatchIndices(unittest.TestCase):
    """Regression tests for cleavage-index computation on non-zero-length matches."""

    def test_multi_residue_motif_uses_end_index(self):
        """A multi-residue motif must cleave *after* the motif, not one past its start.

        Previously non-zero-length matches yielded ``start + 1``, so the two-residue
        motif ``WK`` in ``AWKB`` (match span [1, 3)) incorrectly returned ``2``
        (mid-motif) instead of ``3``.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(list(get_regex_match_indices("AWKB", "WK")), [3])

    def test_single_residue_motif_unchanged(self):
        """Single-residue raw patterns still cleave after the residue."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(list(get_regex_match_indices("AWKB", "K")), [3])

    def test_zero_width_lookaround_uses_start(self):
        """Zero-width lookaround patterns are unaffected and use the match position."""
        self.assertEqual(list(get_regex_match_indices("PEPTIDE", "(?=P)")), [0, 2])


class TestEmptyCleaveOn(unittest.TestCase):
    """Regression: an empty/None cleave_on must behave as non-specific cleavage."""

    def test_empty_cleave_on_is_non_specific(self):
        from peptacular.digestion.core import generate_regex, get_cleavage_sites

        annotation = pt.ProFormaAnnotation.parse("PEPTIDEKR")
        expected = list(range(len(annotation.stripped_sequence) + 1))
        for cleave_on in ("", None):
            sites = list(get_cleavage_sites(annotation, generate_regex(cleave_on=cleave_on)))
            self.assertEqual(sites, expected)


class TestDigest(unittest.TestCase):
    def test_trypsin_cleavage_sites(self):
        """Test getting cleavage sites for trypsin."""
        annotation = pt.ProFormaAnnotation.parse(PROTEIN)
        cleavage_sites = list(annotation.cleavage_sites(pt.PROTEASE_LOOKUP["trypsin"].pattern))
        expected_sites = [23, 31, 35, 42, 45, 47, 50, 56, 60, 62, 70, 73, 79, 88, 96]
        self.assertEqual(cleavage_sites, expected_sites)

    def test_digest_protein_missed_cleavages_2(self):
        """Test protein digestion with 2 missed cleavages."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        spans = annotation.digest(
            enzyme=pt.Proteases.TRYPSIN,
            missed_cleavages=2,
            min_len=0,
            max_len=100,
            semi=False,
        )
        peptides = {annotation[span].serialize() for span in spans}
        expected = {
            "TIDER",
            "TIDERTIDEK",
            "TIDERTIDEKTIDE",
            "TIDEK",
            "TIDEKTIDE",
            "TIDE",
        }
        self.assertEqual(peptides, expected)

    def test_digest_protein_missed_cleavages_1(self):
        """Test protein digestion with 1 missed cleavage."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        spans = annotation.digest(
            enzyme=pt.Proteases.TRYPSIN,
            missed_cleavages=1,
            min_len=0,
            max_len=100,
            semi=False,
        )
        peptides = {annotation[span].serialize() for span in spans}
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_with_k_at_beginning(self):
        """Test protein digestion with K at the beginning."""
        annotation = pt.ProFormaAnnotation.parse("KTIDERTIDEKTIDE")
        spans = annotation.digest(
            enzyme=pt.Proteases.TRYPSIN,
            missed_cleavages=1,
            min_len=0,
            max_len=100,
            semi=False,
        )
        peptides = {annotation[span].serialize() for span in spans}
        expected = {"K", "KTIDER", "TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_with_k_at_end(self):
        """Test protein digestion with K at the end."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDEK")
        spans = annotation.digest(
            enzyme=pt.Proteases.TRYPSIN,
            missed_cleavages=1,
            min_len=0,
            max_len=100,
            semi=False,
        )
        peptides = {annotation[span].serialize() for span in spans}
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDEK", "TIDEK"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_with_double_k(self):
        """Test protein digestion with consecutive K residues."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKKTIDE")
        spans = annotation.digest(
            enzyme=pt.Proteases.TRYPSIN,
            missed_cleavages=1,
            min_len=0,
            max_len=100,
            semi=False,
        )
        peptides = {annotation[span].serialize() for span in spans}
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKK", "K", "KTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_no_missed_cleavages(self):
        """Test protein digestion with no missed cleavages."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        spans = annotation.digest(
            enzyme=pt.Proteases.TRYPSIN,
            missed_cleavages=0,
            min_len=0,
            max_len=100,
            semi=False,
        )
        peptides = {annotation[span].serialize() for span in spans}
        expected = {"TIDER", "TIDEK", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_high_missed_cleavages(self):
        """Test protein digestion with very high missed cleavages."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        spans = annotation.digest(
            enzyme=pt.Proteases.TRYPSIN,
            missed_cleavages=10,
            min_len=0,
            max_len=100,
            semi=False,
        )
        peptides = {annotation[span].serialize() for span in spans}
        expected = {
            "TIDER",
            "TIDERTIDEK",
            "TIDERTIDEKTIDE",
            "TIDEK",
            "TIDEKTIDE",
            "TIDE",
        }
        self.assertEqual(peptides, expected)

    def test_semi_enzymatic_sequences_min_3_max_5(self):
        """Test semi-enzymatic sequences with min_len=3, max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.semi_spans(min_len=3, max_len=5)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"TIDE", "PTIDE", "PEPT", "PEPTI", "IDE", "PEP"}
        self.assertEqual(sequences, expected)

    def test_semi_enzymatic_sequences_min_2_max_5(self):
        """Test semi-enzymatic sequences with min_len=2, max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.semi_spans(min_len=2, max_len=5)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"TIDE", "PTIDE", "PEPT", "PEPTI", "IDE", "PEP", "PE", "DE"}
        self.assertEqual(sequences, expected)

    def test_semi_enzymatic_sequences_min_2_max_100(self):
        """Test semi-enzymatic sequences with min_len=2, max_len=100."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.semi_spans(min_len=2, max_len=100)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {
            "TIDE",
            "PTIDE",
            "PEPT",
            "PEPTI",
            "IDE",
            "PEP",
            "PE",
            "DE",
            "PEPTID",
            "EPTIDE",
        }
        self.assertEqual(sequences, expected)

    def test_semi_enzymatic_sequences_no_limits(self):
        """Test semi-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.semi_spans(min_len=None, max_len=None)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {
            "P",
            "E",
            "TIDE",
            "PTIDE",
            "PEPT",
            "PEPTI",
            "IDE",
            "PEP",
            "PE",
            "DE",
            "PEPTID",
            "EPTIDE",
        }
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_no_limits(self):
        """Test left semi-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.left_semi_spans(None, None)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"P", "PE", "PEP", "PEPT", "PEPTI", "PEPTID"}
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_min_3(self):
        """Test left semi-enzymatic sequences with min_len=3."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.left_semi_spans(3, None)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"PEP", "PEPT", "PEPTI", "PEPTID"}
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_max_5(self):
        """Test left semi-enzymatic sequences with max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.left_semi_spans(None, 5)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"P", "PE", "PEP", "PEPT", "PEPTI"}
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_min_3_max_4(self):
        """Test left semi-enzymatic sequences with min_len=3, max_len=4."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.left_semi_spans(3, 4)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"PEP", "PEPT"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_no_limits(self):
        """Test right semi-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.right_semi_spans(None, None)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"EPTIDE", "PTIDE", "TIDE", "IDE", "DE", "E"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_min_3(self):
        """Test right semi-enzymatic sequences with min_len=3."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.right_semi_spans(3, None)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"EPTIDE", "PTIDE", "TIDE", "IDE"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_max_5(self):
        """Test right semi-enzymatic sequences with max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.right_semi_spans(None, 5)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"PTIDE", "TIDE", "IDE", "DE", "E"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_min_3_max_4(self):
        """Test right semi-enzymatic sequences with min_len=3, max_len=4."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = annotation.right_semi_spans(3, 4)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"TIDE", "IDE"}
        self.assertEqual(sequences, expected)

    def test_non_enzymatic_sequences_no_limits(self):
        """Test non-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        spans = annotation.nonspecific_spans()
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"P", "E", "P", "T", "PE", "EP", "PT", "PEP", "EPT", "PEPT"}
        self.assertEqual(sequences, expected)

    def test_non_enzymatic_sequences_min_1_max_2(self):
        """Test non-enzymatic sequences with min_len=1, max_len=2."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        spans = annotation.nonspecific_spans(min_len=1, max_len=2)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"P", "E", "P", "T", "PE", "EP", "PT"}
        self.assertEqual(sequences, expected)

    def test_non_enzymatic_sequences_min_2_max_4(self):
        """Test non-enzymatic sequences with min_len=2, max_len=4."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        spans = annotation.nonspecific_spans(min_len=2, max_len=4)
        sequences = {annotation[span].serialize() for span in spans}
        expected = {"PT", "EP", "EPT", "PE", "PEP", "PEPT"}
        self.assertEqual(sequences, expected)

    def test_modified_sequences_with_missed_cleavages(self):
        """Test digestion with modified sequences and missed cleavages."""
        modified_seq = "[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)

        spans = annotation.digest(enzyme="(?<=[KR])", missed_cleavages=1)
        peptides = {annotation[span].serialize() for span in spans}
        expected = {
            "[Acetyl]-PEPTIDER[Phospho]",
            "[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K",
            "TIDEM[Oxidation]K",
        }
        self.assertEqual(peptides, expected)

    def test_modified_sequences_no_missed_cleavages(self):
        """Test digestion with modified sequences and no missed cleavages."""
        modified_seq = "[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)

        spans = annotation.digest(enzyme="(?<=[KR])", missed_cleavages=0)
        peptides = {annotation[span].serialize() for span in spans}
        expected = {"[Acetyl]-PEPTIDER[Phospho]", "TIDEM[Oxidation]K"}
        self.assertEqual(peptides, expected)

    def test_empty_sequence_regex_digest(self):
        """Test regex digestion with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(empty_annotation.digest(enzyme="(?<=[KR])"))
        self.assertEqual(result, [])

    def test_empty_sequence_semi_enzymatic(self):
        """Test semi-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(empty_annotation.semi_spans())
        self.assertEqual(result, [])

    def test_empty_sequence_left_semi_enzymatic(self):
        """Test left semi-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(empty_annotation.left_semi_spans())
        self.assertEqual(result, [])

    def test_empty_sequence_right_semi_enzymatic(self):
        """Test right semi-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(empty_annotation.right_semi_spans())
        self.assertEqual(result, [])

    def test_empty_sequence_non_enzymatic(self):
        """Test non-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(empty_annotation.nonspecific_spans())
        self.assertEqual(result, [])

    def test_single_letter_sequence_regex_digest(self):
        """Test regex digestion with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        spans = single_annotation.digest(enzyme="(?<=[KR])")
        result = list(single_annotation[span].serialize() for span in spans)
        self.assertEqual(result, ["K"])

    def test_single_letter_sequence_semi_enzymatic(self):
        """Test semi-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(single_annotation.semi_spans())
        self.assertEqual(result, [])

    def test_single_letter_sequence_left_semi_enzymatic(self):
        """Test left semi-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(single_annotation.left_semi_spans())
        self.assertEqual(result, [])

    def test_single_letter_sequence_right_semi_enzymatic(self):
        """Test right semi-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(single_annotation.right_semi_spans())
        self.assertEqual(result, [])

    def test_single_letter_sequence_non_enzymatic(self):
        """Test non-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(single_annotation.nonspecific_spans())
        # A single residue is itself a valid non-specific peptide (regression: used to be empty)
        self.assertEqual([single_annotation[sp].serialize() for sp in result], ["K"])

    def test_no_cleavage_sites_sequence(self):
        """Test digestion with sequence that has no cleavage sites."""
        peptide_annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        spans = peptide_annotation.digest(enzyme=pt.Proteases.TRYPSIN)
        result = list(peptide_annotation[span].serialize() for span in spans)
        self.assertEqual(result, ["PEPTIDE"])

    def test_sequential_digest_basic(self):
        """Test basic sequential digestion with trypsin and asp-n."""
        trypsin = pt.EnzymeConfig(
            enzyme_regex="(?<=[KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        asp_n = pt.EnzymeConfig(
            enzyme_regex="(?<=\\w)(?=D)",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        annotation = pt.ProFormaAnnotation.parse("PDEREKPKP")
        spans = annotation.sequential_digest(enzyme_configs=[trypsin, asp_n])
        peptides = {annotation[span].serialize() for span in spans}

        expected = {"P", "DER", "EK", "PK", "P"}
        self.assertEqual(peptides, expected)

    def test_sequential_digest_partial_digestion(self):
        """Test sequential digestion with partial digestion enabled."""
        partial_trypsin = pt.EnzymeConfig(
            enzyme_regex="(?<=[KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=False,
        )

        partial_asp_n = pt.EnzymeConfig(
            enzyme_regex="(?<=\\w)(?=D)",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=False,
        )

        annotation = pt.ProFormaAnnotation.parse("PEPRDK")
        spans = annotation.sequential_digest(enzyme_configs=[partial_trypsin, partial_asp_n])
        peptides = {annotation[span].serialize() for span in spans}

        expected = {"PEPR", "DK", "PEPRDK"}
        self.assertEqual(peptides, expected)

    def test_sequential_digest_modified_sequence(self):
        """Test sequential digestion with modified sequence."""
        trypsin = pt.EnzymeConfig(
            enzyme_regex="(?<=[KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        asp_n = pt.EnzymeConfig(
            enzyme_regex="\\w(?=D)",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        annotation = pt.ProFormaAnnotation.parse("PEP[Phospho]RDK")
        spans = annotation.sequential_digest(enzyme_configs=[trypsin, asp_n])
        peptides = {annotation[span].serialize() for span in spans}
        expected = {"PEP[Phospho]R", "DK"}
        self.assertEqual(peptides, expected)

    def test_digest_return_type_span(self):
        """Test digest with span return type."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")

        spans = list(
            annotation.digest(
                enzyme=pt.Proteases.TRYPSIN,
                missed_cleavages=0,
            )
        )

        # Should return Span tuples
        self.assertTrue(len(spans) > 0)
        for item in spans:
            self.assertIsInstance(item, tuple)

    def test_modified_sequence_left_semi_enzymatic(self):
        """Test left semi-enzymatic with modifications preserved."""
        modified_seq = "[1]-P[2]EPTIDE"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        spans = annotation.left_semi_spans()
        sequences = [annotation[span].serialize() for span in spans]
        expected = [
            "[1]-P[2]EPTID",
            "[1]-P[2]EPTI",
            "[1]-P[2]EPT",
            "[1]-P[2]EP",
            "[1]-P[2]E",
            "[1]-P[2]",
        ]
        self.assertEqual(sequences, expected)

    def test_isotopic_labeling_left_semi_enzymatic(self):
        """Test left semi-enzymatic with isotopic labeling."""
        labeled_seq = "<13C>TIDE"
        annotation = pt.ProFormaAnnotation.parse(labeled_seq)
        spans = annotation.left_semi_spans()
        sequences = [annotation[span].serialize() for span in spans]
        expected = ["<13C>TID", "<13C>TI", "<13C>T"]
        self.assertEqual(sequences, expected)

    def test_modified_sequence_right_semi_enzymatic(self):
        """Test right semi-enzymatic with modifications preserved."""
        modified_seq = "PEPTIDE[1]-[2]"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        spans = annotation.right_semi_spans()
        sequences = [annotation[span].serialize() for span in spans]
        expected = [
            "EPTIDE[1]-[2]",
            "PTIDE[1]-[2]",
            "TIDE[1]-[2]",
            "IDE[1]-[2]",
            "DE[1]-[2]",
            "E[1]-[2]",
        ]
        self.assertEqual(sequences, expected)

    def test_isotopic_labeling_right_semi_enzymatic(self):
        """Test right semi-enzymatic with isotopic labeling."""
        labeled_seq = "<13C>TIDE"
        annotation = pt.ProFormaAnnotation.parse(labeled_seq)
        spans = annotation.right_semi_spans()
        sequences = [annotation[span].serialize() for span in spans]
        expected = ["<13C>IDE", "<13C>DE", "<13C>E"]
        self.assertEqual(sequences, expected)

    def test_modified_non_enzymatic_sequences(self):
        """Test non-enzymatic sequences with modifications preserved."""
        modified_seq = "[Acetyl]-P[1.0]EP[1.0]-[Amide]"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        spans = annotation.nonspecific_spans()
        sequences = [annotation[span].serialize() for span in spans]
        expected = [
            "[Acetyl]-P[1.0]",
            "[Acetyl]-P[1.0]E",
            "[Acetyl]-P[1.0]EP[1.0]-[Amide]",
            "E",
            "EP[1.0]-[Amide]",
            "P[1.0]-[Amide]",
        ]
        self.assertEqual(sequences, expected)

    def test_isotopic_labeling_non_enzymatic(self):
        """Test non-enzymatic sequences with isotopic labeling."""
        labeled_seq = "<13C>PEP"
        annotation = pt.ProFormaAnnotation.parse(labeled_seq)
        spans = annotation.nonspecific_spans()
        sequences = [annotation[span].serialize() for span in spans]
        expected = ["<13C>P", "<13C>PE", "<13C>PEP", "<13C>E", "<13C>EP", "<13C>P"]
        self.assertEqual(sequences, expected)

    def test_cleavage_sites_enzyme_key_trypsin_p(self):
        """Test cleavage sites using trypsin/P enzyme key."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        sites = list(annotation.cleavage_sites(pt.Proteases.TRYPSIN))
        expected = [5, 10]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_no_match(self):
        """Test cleavage sites with enzyme that doesn't match sequence."""
        annotation = pt.ProFormaAnnotation.parse("TIDEPTIDEPTIDE")
        sites = list(annotation.cleavage_sites(pt.Proteases.TRYPSIN))
        expected = []
        self.assertEqual(sites, expected)

    def test_cleavage_sites_n_terminal(self):
        """Test cleavage sites with N-terminal cleavage (lys-n)."""
        annotation = pt.ProFormaAnnotation.parse("KPEPTIDEK")
        sites = list(annotation.cleavage_sites(pt.Proteases.LYS_N))
        expected = [0, 8]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_c_terminal(self):
        """Test cleavage sites with C-terminal cleavage (lys-c)."""
        annotation = pt.ProFormaAnnotation.parse("KPEPTIDEK")
        sites = list(annotation.cleavage_sites(pt.Proteases.LYS_C))
        expected = [1, 9]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_modified_sequence(self):
        """Test cleavage sites with modified sequence."""
        modified_seq = "[Acetyl]-TIDERT[1.0]IDEKTIDE-[Amide]"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        sites = list(annotation.cleavage_sites(pt.Proteases.TRYPSIN))
        expected = [5, 10]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_non_specific(self):
        """Test cleavage sites with non-specific cleavage."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sites = list(annotation.cleavage_sites(pt.Proteases.UNSPECIFIC))
        expected = [0, 1, 2, 3, 4, 5, 6, 7]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_overlapping_pattern(self):
        """Test cleavage sites with overlapping pattern."""
        annotation = pt.ProFormaAnnotation.parse("PPPP")
        sites = list(annotation.cleavage_sites("(?=PP)"))
        expected = [0, 1, 2]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_lookahead_pattern(self):
        """Test cleavage sites with lookahead pattern."""
        annotation = pt.ProFormaAnnotation.parse("PEPCTIDE")
        sites = list(annotation.cleavage_sites("(?=C)"))
        expected = [3]
        self.assertEqual(sites, expected)

    def test_digest_with_string_span_return_type(self):
        """Test digest with string-span return type."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        spans = list(annotation.digest(enzyme=pt.Proteases.TRYPSIN))
        expected = [(0, 5, 0), (5, 10, 0), (10, 14, 0)]
        self.assertEqual(spans, expected)

    def test_digest_complex_modified_sequence(self):
        """Test digest with complex modified sequence and length filters."""
        seq = "<13C>T[1][2]IDERTIDEKTIDE"
        annotation = pt.ProFormaAnnotation.parse(seq)
        spans = annotation.digest(enzyme="(?<=[KR])", missed_cleavages=2, min_len=6)
        result = [annotation[span].serialize() for span in spans]
        expected = [
            "<13C>T[1][2]IDERTIDEK",
            "<13C>T[1][2]IDERTIDEKTIDE",
            "<13C>TIDEKTIDE",
        ]
        self.assertEqual(result, expected)

    def test_digest_semi_enzymatic_with_modified_sequence(self):
        """Test digest with semi-enzymatic and modified sequence."""
        seq = "TIDERTIDEK[1]TIDE-[2]"
        annotation = pt.ProFormaAnnotation.parse(seq)
        spans = annotation.digest(
            enzyme="(?<=[KR])",
            missed_cleavages=1,
            min_len=9,
            semi=True,
        )
        result = [annotation[span].serialize() for span in spans]
        expected = ["TIDERTIDE", "TIDERTIDEK[1]", "IDERTIDEK[1]", "TIDEK[1]TIDE-[2]"]
        self.assertEqual(result, expected)

    def test_digest_non_specific_enzyme(self):
        """Test digest with non-specific enzyme."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        spans = annotation.digest(enzyme=pt.Proteases.UNSPECIFIC)
        result = [annotation[span].serialize() for span in spans]
        expected = ["P", "PE", "PEP", "PEPT", "E", "EP", "EPT", "P", "PT", "T"]
        self.assertEqual(result, expected)

    def test_sequential_digest_example_xxxkxxxdxxx(self):
        """Test sequential digest with example from docstring."""
        trypsin = pt.EnzymeConfig(
            enzyme_regex="(?<=[KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )
        asp_n = pt.EnzymeConfig(
            enzyme_regex="(?<=[D])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        annotation = pt.ProFormaAnnotation.parse("XXXKXXXDXXX")
        spans = annotation.sequential_digest(enzyme_configs=[trypsin, asp_n])
        result = [annotation[span].serialize() for span in spans]
        expected = ["XXXK", "XXXD", "XXX"]
        self.assertEqual(result, expected)

    def test_sequential_digest_partial_digestion_example(self):
        """Test sequential digest with partial digestion example."""
        partial_digest = [
            pt.EnzymeConfig(
                enzyme_regex="(?<=[KR])",
                missed_cleavages=0,
                semi_enzymatic=False,
                complete_digestion=False,
            ),
            pt.EnzymeConfig(
                enzyme_regex="(?<=[D])",
                missed_cleavages=0,
                semi_enzymatic=False,
                complete_digestion=False,
            ),
        ]

        annotation = pt.ProFormaAnnotation.parse("XXXKXXXDXXX")
        spans = annotation.sequential_digest(enzyme_configs=partial_digest)
        result = [annotation[span].serialize() for span in spans]
        expected = ["XXXK", "XXXKXXXD", "XXXKXXXDXXX", "XXX", "XXXD", "XXXDXXX", "XXX"]
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
