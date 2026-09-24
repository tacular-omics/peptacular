"""Every protease tacular knows resolves in peptacular, and unknown names say how to use a regex."""

import re

import pytest
from tacular import PROTEASE_LOOKUP, Protease

import peptacular as pt
from peptacular.digestion.core import resolve_enzyme

_ENTRIES = list(PROTEASE_LOOKUP)


@pytest.mark.parametrize("info", _ENTRIES, ids=[str(e.id.value) for e in _ENTRIES])
def test_every_protease_resolves_by_id_name_member_and_value(info):
    expected = info.pattern.pattern
    for key in (info.id, info.id.value, info.name):
        assert resolve_enzyme(key).pattern == expected, key


@pytest.mark.parametrize("member", list(Protease), ids=lambda m: m.value)
def test_every_protease_member_digests(member):
    seq = "MKRPEPTIDEKDRSTYWFLAGN"
    peptides = pt.digest(seq, member)
    assert peptides
    assert all(seq[span.start : span.end] == p for p, span in peptides)
    if member is not Protease.UNSPECIFIC:
        assert "".join(p for p, _ in peptides) == seq


def test_unknown_enzyme_hint_suggests_re_compile():
    with pytest.raises(pt.UnknownEnzymeError, match=re.escape("re.compile('[KR]')")) as info:
        pt.digest("PEPTIDEK", "[KR]")
    assert isinstance(info.value, KeyError)
    assert "trypsin" in str(info.value)


def test_compiled_pattern_passes_through():
    pattern = re.compile("(?<=[KR])")
    assert resolve_enzyme(pattern) is pattern
    assert [p for p, _ in pt.digest("AKBRC", pattern)] == ["AK", "BR", "C"]


def test_non_string_enzyme_is_type_error():
    with pytest.raises(TypeError):
        resolve_enzyme(42)  # ty: ignore[invalid-argument-type]
