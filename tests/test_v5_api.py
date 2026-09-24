"""Behaviour introduced by the 5.0 API cleanup (see docs/migration.rst)."""

import dataclasses
import re

import pytest

import peptacular as pt
from peptacular.digestion.core import resolve_enzyme


class TestDigestNaming:
    def test_annotation_span_methods_renamed(self):
        annot = pt.parse("PEPTIDEKAAR")
        for old in ("digest", "simple_digest", "sequential_digest"):
            assert not hasattr(annot, old), old
        assert [annot[s].serialize() for s in annot.digest_spans("trypsin")] == ["PEPTIDEK", "AAR"]
        assert [annot[s].serialize() for s in annot.simple_digest_spans(cleave_on="K")] == ["PEPTIDEK", "AAR"]

    def test_functional_digest_agrees_with_span_method(self):
        annot = pt.parse("PEPTIDEKAAR")
        assert pt.digest(annot, "trypsin") == [(annot[s].serialize(), s) for s in annot.digest_spans("trypsin")]

    def test_sequential_digest_spans(self):
        annot = pt.parse("AAKDBBRDCC")
        configs = [pt.EnzymeConfig(enzyme="trypsin"), pt.EnzymeConfig(enzyme=re.compile("(?=D)"))]
        assert [annot[s].serialize() for s in annot.sequential_digest_spans(configs)] == ["AAK", "DBBR", "DCC"]

    def test_unspecific_protease_cuts_everywhere(self):
        assert pt.cleavage_sites("PEP", "unspecific") == [0, 1, 2, 3]


class TestEnzymeConfig:
    def test_field_renamed(self):
        names = [f.name for f in dataclasses.fields(pt.EnzymeConfig)]
        assert names[0] == "enzyme"
        assert "enzyme_regex" not in names

    def test_frozen_and_slotted(self):
        cfg = pt.EnzymeConfig(enzyme="trypsin")
        with pytest.raises(dataclasses.FrozenInstanceError):
            cfg.enzyme = "lys_c"  # ty: ignore[invalid-assignment]
        assert not hasattr(cfg, "__dict__")


class TestResolveEnzyme:
    def test_pattern_passthrough(self):
        pattern = re.compile("(?<=K)")
        assert resolve_enzyme(pattern) is pattern

    def test_name_and_member(self):
        assert resolve_enzyme("trypsin") is resolve_enzyme(pt.Proteases.TRYPSIN)

    def test_unknown_error_message_names_known_proteases(self):
        with pytest.raises(pt.UnknownEnzymeError) as info:
            resolve_enzyme("trypsn")
        assert "trypsin" in str(info.value)
        assert not str(info.value).startswith("'")  # not KeyError's repr() formatting

    def test_unknown_enzyme_diagnostic_code(self):
        diag = pt.diagnose("PEPTIDE", "digest", enzyme="trypsn")
        assert diag is not None
        assert diag.code == "unknown_enzyme"


class TestKeywordOnlyParallelArgs:
    def test_n_workers_is_keyword_only(self):
        with pytest.raises(TypeError):
            pt.mass(["PEPTIDE"], "p", 0, True, None, None, False, 2)  # ty: ignore[too-many-positional-arguments]
        assert pt.mass(["PEPTIDE"], n_workers=1) == [pt.mass("PEPTIDE")]


def test_isotopic_distribution_takes_sequence_keyword():
    assert pt.isotopic_distribution(sequence="PEPTIDE", max_isotopes=3) == pt.isotopic_distribution("PEPTIDE", max_isotopes=3)
    with pytest.raises(TypeError):
        pt.isotopic_distribution(annotations="PEPTIDE")  # ty: ignore[unknown-argument]


class TestPeptacularErrors:
    def test_no_bare_value_or_index_errors_in_library_code(self):
        # 5.0: input errors raise PeptacularError (a ValueError) or a subclass.
        # The MCP layer is excluded: pydantic validators raise ValueError by contract.
        import ast
        import pathlib

        root = pathlib.Path(pt.__file__).parent
        offenders = []
        for path in root.rglob("*.py"):
            if "mcp" in path.relative_to(root).parts:
                continue
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.Raise) and isinstance(node.exc, ast.Call) and isinstance(node.exc.func, ast.Name):
                    if node.exc.func.id in {"ValueError", "IndexError", "KeyError"}:
                        offenders.append(f"{path.relative_to(root)}:{node.lineno}")
        assert offenders == []

    @pytest.mark.parametrize(
        "call",
        [
            lambda: pt.Terminal.from_str("middle"),
            lambda: pt.calc_property("PEPTIDE", scale="no_such_scale"),
            lambda: pt.estimate_isotopic_distribution(-1.0),
            lambda: pt.batch("no_such_operation", ["PEPTIDE"]),
            lambda: pt.parse("PEP[").serialize(),
        ],
    )
    def test_main_paths_raise_peptacular_error(self, call):
        with pytest.raises(pt.PeptacularError):
            call()

    def test_internal_index_errors_are_invalid_position_errors(self):
        assert issubclass(pt.InvalidPositionError, IndexError)
        assert issubclass(pt.InvalidPositionError, pt.PeptacularError)


def test_annotation_simple_cleavage_sites_uses_generated_pattern():
    annot = pt.parse("PEPTIDEKAAR")
    assert list(annot.simple_cleavage_sites("KR")) == [8, 11]
    assert list(annot.simple_cleavage_sites("KR")) == pt.simple_cleavage_sites("PEPTIDEKAAR", "KR")
