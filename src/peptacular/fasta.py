"""FASTA parsing with bounded memory and explicit stream ownership."""

import codecs
import gzip
import io
import pathlib
from collections.abc import Iterable, Iterator
from typing import BinaryIO, NamedTuple, Protocol, cast, runtime_checkable

FASTA_INPUT_TYPE = str | pathlib.Path | io.IOBase

__all__ = ["FASTA_INPUT_TYPE", "FastaSequence", "ReadableProtocol", "iter_fasta", "parse_fasta", "parse_fasta_text"]


@runtime_checkable
class ReadableProtocol(Protocol):
    def read(self, size: int = -1) -> str | bytes: ...


class FastaSequence(NamedTuple):
    header: str
    sequence: str


def _open_binary(path: pathlib.Path):
    return gzip.open(path, "rb") if path.suffix.lower() == ".gz" else path.open("rb")


def _detect_file_encoding(file_path: str | pathlib.Path) -> str:
    with _open_binary(pathlib.Path(file_path)) as stream:
        sample = stream.read(8192)
    if sample.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        return "utf-16"
    try:
        codecs.getincrementaldecoder("utf-8-sig")().decode(sample, final=False)
        return "utf-8-sig"
    except UnicodeDecodeError:
        return "latin-1"


def _iter_fasta_lines(lines: Iterable[str]) -> Iterator[FastaSequence]:
    header: str | None = None
    seq_lines: list[str] = []
    found_text = False
    emitted = False
    for line_number, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue
        found_text = True
        if line.startswith(">"):
            if header is not None and seq_lines:
                emitted = True
                yield FastaSequence(header, "".join(seq_lines).upper())
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Empty header found at line {line_number}")
            seq_lines = []
        else:
            if header is None:
                raise ValueError(f"Sequence data before header at line {line_number}")
            seq_lines.append(line)
    if header is not None and seq_lines:
        emitted = True
        yield FastaSequence(header, "".join(seq_lines).upper())
    if not found_text:
        raise ValueError("Empty input text")
    if not emitted:
        raise ValueError("No valid FASTA sequences found")


def iter_fasta(input_data: FASTA_INPUT_TYPE, *, encoding: str | None = None) -> Iterator[FastaSequence]:
    """Yield FASTA records from text, a path, or an open text/binary stream.

    Memory usage is proportional to the largest record. Paths ending in ``.gz``
    are decompressed automatically. Path encodings are sampled when ``encoding``
    is omitted. Binary streams default to UTF-8 with optional BOM removal.
    Text streams are already decoded, so ``encoding`` does not apply to them.

    Files opened here close on exhaustion, error, or iterator ``close()``.
    Caller-owned streams remain open. A binary wrapper may read ahead, so its
    stream position after early termination is not a record boundary guarantee.
    Use ``contextlib.closing`` when stopping an iterator over a path early.
    Empty records are skipped for compatibility. Errors in later records are
    raised when iteration reaches them, after earlier records may be yielded.

    :param input_data: FASTA text, a filesystem path, or an open stream.
    :param encoding: Optional file or binary-stream encoding override.
    :return: Iterator of header and uppercase sequence records.
    :raises ValueError: Empty input or invalid record structure.
    """
    if isinstance(input_data, str) and (not input_data.strip() or input_data.lstrip().startswith(">") or "\n" in input_data or "\r" in input_data):
        with io.StringIO(input_data, newline=None) as stream:
            yield from _iter_fasta_lines(stream)
    elif isinstance(input_data, (str, pathlib.Path)):
        path = pathlib.Path(input_data)
        selected_encoding = encoding or _detect_file_encoding(path)
        with _open_binary(path) as binary, io.TextIOWrapper(binary, encoding=selected_encoding) as stream:
            yield from _iter_fasta_lines(stream)
    elif isinstance(input_data, io.IOBase):
        if isinstance(input_data, (io.BufferedIOBase, io.RawIOBase)):
            wrapper = io.TextIOWrapper(cast(BinaryIO, input_data), encoding=encoding or "utf-8-sig")
            try:
                yield from _iter_fasta_lines(wrapper)
            finally:
                wrapper.detach()
        else:
            yield from _iter_fasta_lines(cast(Iterable[str], input_data))
    else:
        raise TypeError(f"Unsupported input type: {type(input_data)}")


def parse_fasta(input_data: FASTA_INPUT_TYPE, *, encoding: str | None = None) -> list[FastaSequence]:
    """Collect :func:`iter_fasta` into a list, preserving the existing return type."""
    # Preserve support for legacy read-only adapters. Streaming adapters should
    # expose a standard text or binary IOBase instead of only read().
    if not isinstance(input_data, (str, pathlib.Path, io.IOBase)) and isinstance(input_data, ReadableProtocol):
        content = input_data.read()
        text = content.decode(encoding or "utf-8-sig") if isinstance(content, bytes) else content
        return parse_fasta_text(text)
    return list(iter_fasta(input_data, encoding=encoding))


def parse_fasta_text(text: str) -> list[FastaSequence]:
    """Parse FASTA text with the same validation as :func:`iter_fasta`."""
    with io.StringIO(text, newline=None) as stream:
        return list(_iter_fasta_lines(stream))
