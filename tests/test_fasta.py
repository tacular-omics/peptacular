"""Tests for FASTA parsing (peptacular.fasta)."""

import io

import pytest

from peptacular.fasta import FastaSequence, parse_fasta, parse_fasta_text


class TestParseFastaText:
    def test_single_record(self):
        recs = parse_fasta_text(">sp|P1|NAME description\nPEPTIDE\n")
        assert recs == [FastaSequence(header="sp|P1|NAME description", sequence="PEPTIDE")]

    def test_multiple_records(self):
        recs = parse_fasta_text(">A\nPEPTIDE\n>B\nMKLV\n")
        assert [r.header for r in recs] == ["A", "B"]
        assert [r.sequence for r in recs] == ["PEPTIDE", "MKLV"]

    def test_wrapped_sequence_lines_are_joined(self):
        recs = parse_fasta_text(">h\nPEPT\nIDEK\nMKLV\n")
        assert recs[0].sequence == "PEPTIDEKMKLV"

    def test_lowercase_is_uppercased(self):
        recs = parse_fasta_text(">h\npeptide\n")
        assert recs[0].sequence == "PEPTIDE"

    def test_blank_lines_ignored(self):
        recs = parse_fasta_text(">h\n\nPEPTIDE\n\n")
        assert recs == [FastaSequence("h", "PEPTIDE")]

    def test_no_trailing_newline(self):
        recs = parse_fasta_text(">h\nPEPTIDE")
        assert recs[0].sequence == "PEPTIDE"

    def test_crlf_line_endings(self):
        recs = parse_fasta_text(">h\r\nPEPT\r\nIDE\r\n")
        assert recs == [FastaSequence("h", "PEPTIDE")]

    def test_header_with_no_sequence_is_dropped(self):
        # A header immediately followed by another header yields only the record
        # that actually has sequence data.
        recs = parse_fasta_text(">empty\n>real\nPEPTIDE\n")
        assert recs == [FastaSequence("real", "PEPTIDE")]

    def test_empty_text_raises(self):
        with pytest.raises(ValueError, match="Empty input text"):
            parse_fasta_text("   \n  \n")

    def test_empty_header_raises(self):
        with pytest.raises(ValueError, match="Empty header"):
            parse_fasta_text(">\nPEPTIDE\n")

    def test_sequence_before_header_raises(self):
        with pytest.raises(ValueError, match="Sequence data before header"):
            parse_fasta_text("PEPTIDE\n>h\nMKLV\n")

    def test_only_headers_no_sequences_raises(self):
        with pytest.raises(ValueError, match="No valid FASTA sequences"):
            parse_fasta_text(">a\n>b\n")


class TestParseFastaInputs:
    def test_parse_from_text_string(self):
        recs = parse_fasta(">h\nPEPTIDE\n")
        assert recs[0] == FastaSequence("h", "PEPTIDE")

    def test_parse_from_file_path_str(self, tmp_path):
        p = tmp_path / "seqs.fasta"
        p.write_text(">h1\nPEPTIDE\n>h2\nMKLV\n")
        recs = parse_fasta(str(p))
        assert [r.header for r in recs] == ["h1", "h2"]

    def test_parse_from_pathlib_path(self, tmp_path):
        p = tmp_path / "seqs.fasta"
        p.write_text(">h\nPEPTIDE\n")
        recs = parse_fasta(p)
        assert recs[0].sequence == "PEPTIDE"

    def test_parse_from_stream(self):
        stream = io.StringIO(">h\nPEPTIDE\n")
        recs = parse_fasta(stream)
        assert recs[0].sequence == "PEPTIDE"

    def test_parse_from_bytes_stream(self):
        stream = io.BytesIO(b">h\nPEPTIDE\n")
        recs = parse_fasta(stream)
        assert recs[0].sequence == "PEPTIDE"

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, ValueError)):
            parse_fasta(tmp_path / "does_not_exist.fasta")

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError):
            parse_fasta(12345)  # type: ignore[arg-type]
