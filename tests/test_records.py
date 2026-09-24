"""digest_records / fragment_records: plain-dict rows for pandas or polars."""

import dataclasses
import importlib.util

import pytest

import peptacular as pt
from peptacular.sequence.records import _format_counts


@dataclasses.dataclass
class Entry:
    """Stand-in for a fastatacular/pefftacular entry: duck-typed .sequence and .accession."""

    sequence: str
    accession: str | None


SCALARS = (str, int, float, bool, type(None))


class TestDigestRecords:
    def test_key_set_and_order(self):
        rows = pt.digest_records("MKVLATSAGERTIDEK", "trypsin", missed_cleavages=1)
        assert rows
        for row in rows:
            assert tuple(row) == pt.DIGEST_RECORD_FIELDS
            assert all(isinstance(v, SCALARS) for v in row.values())

    def test_values_equal_digest(self):
        protein = "MKVLATSAGERTIDEK"
        rows = pt.digest_records(protein, "trypsin", missed_cleavages=1)
        pairs = pt.digest(protein, "trypsin", missed_cleavages=1)
        assert len(rows) == len(pairs)
        for row, (peptide, span) in zip(rows, pairs, strict=True):
            assert row["peptide"] == peptide
            assert (row["start"], row["end"], row["missed_cleavages"]) == tuple(span)
            assert row["stripped_sequence"] == protein[row["start"] : row["end"]]
            assert row["semi"] is False
            assert row["accession"] is None

    def test_modified_protein_keeps_proforma_and_strips(self):
        protein = "[Acetyl]-MKVLAT[Phospho]SAGERK"
        rows = pt.digest_records(protein, "trypsin", missed_cleavages=1)
        pairs = pt.digest(protein, "trypsin", missed_cleavages=1)
        assert [r["peptide"] for r in rows] == [p for p, _ in pairs]
        stripped = pt.parse(protein).stripped_sequence
        for row in rows:
            assert row["stripped_sequence"] == stripped[row["start"] : row["end"]]
        assert rows[0]["peptide"] == "[Acetyl]-MK"
        assert any(r["peptide"] == "VLAT[Phospho]SAGER" for r in rows)

    def test_half_open_span(self):
        rows = pt.digest_records("TIDERTIDEKTIDE", "trypsin")
        assert [(r["start"], r["end"]) for r in rows] == [(0, 5), (5, 10), (10, 14)]

    def test_semi_flag(self):
        protein = "MKVLATSAGERTIDEK"
        rows = pt.digest_records(protein, "trypsin", semi=True)
        sites = {0, len(protein), *pt.cleavage_sites(protein, "trypsin")}
        assert any(r["semi"] for r in rows) and any(not r["semi"] for r in rows)
        for row in rows:
            assert row["semi"] == (row["start"] not in sites or row["end"] not in sites)
        full = {(r["start"], r["end"]) for r in pt.digest_records(protein, "trypsin")}
        assert {(r["start"], r["end"]) for r in rows if not r["semi"]} == full

    def test_accession_is_duck_typed(self):
        rows = pt.digest_records(Entry("PEPTIDEKTIDE", "P12345"), "trypsin")
        assert [r["accession"] for r in rows] == ["P12345", "P12345"]
        assert pt.digest_records(Entry("PEPTIDEK", None), "trypsin")[0]["accession"] is None

    def test_non_str_accession_is_stringified(self):
        rows = pt.digest_records(Entry("PEPTIDEK", 42), "trypsin")  # type: ignore[arg-type]
        assert rows[0]["accession"] == "42"

    def test_batch_is_flat_in_input_order(self):
        rows = pt.digest_records([Entry("PEPKTIDE", "A"), Entry("SAMPLERK", "B")], "trypsin")
        assert [(r["accession"], r["peptide"]) for r in rows] == [("A", "PEPK"), ("A", "TIDE"), ("B", "SAMPLER"), ("B", "K")]

    def test_empty_digest(self):
        assert pt.digest_records("PEPTIDE", "trypsin", min_len=50) == []
        assert pt.digest_records([], "trypsin") == []

    def test_annotation_input(self):
        rows = pt.digest_records(pt.parse("PEPKTIDE"), "trypsin")
        assert [r["peptide"] for r in rows] == ["PEPK", "TIDE"]

    def test_unknown_enzyme(self):
        with pytest.raises(pt.UnknownEnzymeError):
            pt.digest_records("PEPTIDE", "not-an-enzyme")


class TestFragmentRecords:
    def test_key_set_and_scalar_values(self):
        rows = pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b", "y", "by", "p"), charges=(1,)))
        for row in rows:
            assert tuple(row) == pt.FRAGMENT_RECORD_FIELDS
            for key, value in row.items():
                if key == "position" and isinstance(value, tuple):
                    assert len(value) == 2 and all(isinstance(v, int) for v in value)
                else:
                    assert isinstance(value, SCALARS), key

    def test_multi_charge_values_equal_fragments(self):
        frags = pt.fragment("PEM[Oxidation]TIDEK/3", ion_types=("b", "y"), charges=(1, 2, 3))
        rows = pt.fragment_records(frags)
        assert len(rows) == len(frags) == 2 * 3 * 8
        assert {r["charge_state"] for r in rows} == {1, 2, 3}
        for row, frag in zip(rows, frags, strict=True):
            assert row["ion_type"] == frag.ion_type.value
            assert row["position"] == frag.position
            assert row["charge_state"] == frag.charge_state
            assert row["mz"] == frag.mz
            assert row["mass"] == frag.mass
            assert row["neutral_mass"] == frag.neutral_mass
            assert row["monoisotopic"] is frag.monoisotopic
            assert row["sequence"] == frag.sequence
            assert row["parent_sequence"] == frag.parent_sequence
            assert row["mzpaf"] == frag.to_mzpaf()
            assert row["losses"] == "" and row["isotopes"] == ""

    def test_losses_and_isotopes(self):
        frags = pt.fragment("PEPTIDE", ion_types=("y",), charges=(1,), neutral_deltas=["H2O"], isotopes=(0, 2, {"15N": 1}))
        rows = pt.fragment_records(frags)
        assert {r["isotopes"] for r in rows} == {"", "13C^2", "15N"}
        assert {r["losses"] for r in rows} == {"", "H-2O-1"}
        one_c13 = pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("y",), charges=(1,), isotopes=(1,)))
        assert {r["isotopes"] for r in one_c13} == {"13C"}

    def test_mass_delta_has_no_mzpaf(self):
        rows = pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b",), charges=(1,), deltas=({-17.0: 2},)))
        assert {r["losses"] for r in rows} == {"-17.0^2"}
        assert {r["mzpaf"] for r in rows} == {None}

    def test_without_parent_sequence(self):
        frag = pt.parse("PEPTIDE").frag("b", 1, position=3, _include_sequence=False)
        (row,) = pt.fragment_records([frag])
        assert row["sequence"] is None and row["parent_sequence"] is None
        assert row["mzpaf"] == "b3"

    def test_empty_and_type_error(self):
        assert pt.fragment_records([]) == []
        with pytest.raises(TypeError, match="Fragment"):
            pt.fragment_records(["b1"])  # type: ignore[list-item]

    def test_format_counts(self):
        assert _format_counts([("H-2O-1", 1), (1.5, 1), (-17.0, 3)]) == "H-2O-1,+1.5,-17.0^3"


class TestNoDataFrameDependency:
    def test_pandas_and_polars_are_not_imported(self):
        # Records are plain dicts; building them must not pull in a DataFrame library.
        pt.digest_records("PEPTIDEK", "trypsin")
        pt.fragment_records(pt.fragment("PEPTIDE"))
        import peptacular.sequence.records as records

        assert "pandas" not in vars(records) and "polars" not in vars(records)

    def test_not_declared_as_dependencies(self):
        import tomllib
        from pathlib import Path

        pyproject = tomllib.loads((Path(__file__).parent.parent / "pyproject.toml").read_text())
        declared = " ".join(pyproject["project"].get("dependencies", []))
        for extra in pyproject["project"].get("optional-dependencies", {}).values():
            declared += " " + " ".join(extra)
        assert "pandas" not in declared and "polars" not in declared

    @pytest.mark.skipif(importlib.util.find_spec("pandas") is None, reason="pandas not installed")
    def test_pandas_round_trip(self):  # pragma: no cover - depends on an optional local install
        import pandas as pd

        frame = pd.DataFrame(pt.digest_records("MKVLATSAGERTIDEK", "trypsin"))
        assert list(frame.columns) == list(pt.DIGEST_RECORD_FIELDS)

    @pytest.mark.skipif(importlib.util.find_spec("polars") is None, reason="polars not installed")
    def test_polars_round_trip(self):  # pragma: no cover - depends on an optional local install
        import polars as pl

        frame = pl.DataFrame(pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b", "y"), charges=(1, 2))))
        assert frame.columns == list(pt.FRAGMENT_RECORD_FIELDS)
