"""Streaming FASTA behavior, compressed files, and resource ownership."""

import gzip
import io
from contextlib import closing

import pytest

import peptacular as pt

TEXT = ">first description\npept\nide\n>second\nMKR\n"


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig", "utf-16", "latin-1"])
def test_files_and_encoding(tmp_path, compressed, encoding):
    text = TEXT.replace("description", "café")
    path = tmp_path / ("data.fasta.gz" if compressed else "data.fasta")
    payload = text.encode(encoding)
    path.write_bytes(gzip.compress(payload) if compressed else payload)
    assert list(pt.iter_fasta(path)) == pt.parse_fasta_text(text)
    assert pt.parse_fasta(str(path)) == pt.parse_fasta_text(text)


def test_explicit_encoding_handles_non_ascii_late_in_file(tmp_path):
    path = tmp_path / "data.fasta"
    path.write_bytes((">first\n" + "A\n" * 5000 + ">café\nPEPTIDE").encode("latin-1"))
    assert list(pt.iter_fasta(path, encoding="latin-1"))[-1].header == "café"


def test_no_bulk_read_or_eager_consumption():
    class LineOnly(io.StringIO):
        def read(self, *args):
            raise AssertionError("The parser must iterate lines")

    stream = LineOnly(TEXT + ">\nINVALID\n")
    with closing(pt.iter_fasta(stream)) as records:
        assert next(records) == pt.FastaSequence("first description", "PEPTIDE")
        assert not stream.closed
        assert next(records) == pt.FastaSequence("second", "MKR")
        with pytest.raises(ValueError, match="line 6"):
            next(records)
    assert not stream.closed


@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("finish", ["exhaust", "close", "error"])
def test_caller_stream_remains_open(binary, finish):
    text = TEXT if finish != "error" else "PEPTIDE\n"
    stream = io.BytesIO(text.encode()) if binary else io.StringIO(text)
    records = pt.iter_fasta(stream)
    if finish == "exhaust":
        list(records)
    elif finish == "close":
        next(records)
        records.close()
    else:
        with pytest.raises(ValueError):
            next(records)
    assert not stream.closed


def test_file_iterator_closes_its_file(tmp_path):
    path = tmp_path / "data.fasta"
    path.write_text(TEXT)
    with closing(pt.iter_fasta(path)) as records:
        next(records)
    # Windows disallows deletion while the file is still open.
    path.unlink()


def test_text_empty_records_and_newlines_remain_compatible():
    text = ">empty\r>first\rpep\rtide\r>last\rMKR"
    assert list(pt.iter_fasta(text)) == pt.parse_fasta_text(text)
    assert [r.sequence for r in pt.iter_fasta(text)] == ["PEPTIDE", "MKR"]


def test_legacy_read_only_adapter_still_works():
    class Reader:
        def read(self):
            return TEXT

    assert pt.parse_fasta(Reader()) == pt.parse_fasta_text(TEXT)
