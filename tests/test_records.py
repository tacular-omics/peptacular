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


class PeffEntry:
    """Stand-in for a pefftacular entry: no .accession, the id is in .db_unique_id."""

    def __init__(self, prefix: str, db_unique_id: str, sequence: str) -> None:
        self.prefix = prefix
        self.db_unique_id = db_unique_id
        self.sequence = sequence


SCALARS = (str, int, float, bool, type(None))


class TestDigestRecords:
    def test_key_set_and_order(self):
        rows = pt.digest_records("MKVLATSAGERTIDEK", "trypsin", missed_cleavages=1)
        assert rows
        for row in rows:
            assert tuple(row) == pt.DIGEST_RECORD_KEYS
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

    def test_peff_entry_uses_db_unique_id(self):
        rows = pt.digest_records(PeffEntry("sp", "P12345", "PEPTIDEKTIDE"), "trypsin")
        assert [r["accession"] for r in rows] == ["P12345", "P12345"]

    def test_accession_wins_over_db_unique_id(self):
        entry = Entry("PEPTIDEK", "ACC")
        entry.db_unique_id = "UID"  # type: ignore[attr-defined]
        assert pt.digest_records(entry, "trypsin")[0]["accession"] == "ACC"

    def test_generator_input(self):
        entries = (e for e in [Entry("PEPKTIDE", "A"), PeffEntry("sp", "B", "SAMPLERK")])
        rows = pt.digest_records(entries, "trypsin")
        assert [(r["accession"], r["peptide"]) for r in rows] == [("A", "PEPK"), ("A", "TIDE"), ("B", "SAMPLER"), ("B", "K")]
        assert pt.digest_records(iter(["PEPKTIDE"]), "trypsin")[1]["peptide"] == "TIDE"

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
            assert tuple(row) == pt.FRAGMENT_RECORD_KEYS
            assert all(isinstance(v, SCALARS) for v in row.values())
            assert row["position"] is None or isinstance(row["position"], int)
            assert row["end_position"] is None or isinstance(row["end_position"], int)

    def test_internal_ion_start_and_end(self):
        frags = pt.fragment("PEPTIDE", ion_types=("by", "b", "p"), charges=(1,))
        rows = pt.fragment_records(frags)
        for row, frag in zip(rows, frags, strict=True):
            if isinstance(frag.position, tuple):
                assert (row["position"], row["end_position"]) == frag.position
            else:
                assert row["position"] == frag.position and row["end_position"] is None
        assert any(r["end_position"] is not None for r in rows)
        assert [r for r in rows if r["ion_type"] == "p"][0]["position"] is None

    def test_multi_charge_values_equal_fragments(self):
        frags = pt.fragment("PEM[Oxidation]TIDEK/3", ion_types=("b", "y"), charges=(1, 2, 3))
        rows = pt.fragment_records(frags)
        assert len(rows) == len(frags) == 2 * 3 * 8
        assert {r["charge_state"] for r in rows} == {1, 2, 3}
        for row, frag in zip(rows, frags, strict=True):
            assert row["ion_type"] == frag.ion_type.value
            assert row["position"] == frag.position and row["end_position"] is None
            assert row["charge_state"] == frag.charge_state
            assert row["mz"] == frag.mz
            assert row["mass"] == frag.mass
            assert row["neutral_mass"] == frag.neutral_mass
            assert row["monoisotopic"] is frag.monoisotopic
            assert frag.sequence is not None and frag.parent_sequence is not None
            assert row["sequence"] == pt.parse(frag.sequence).serialize(exclude_charge=True)
            assert row["parent_sequence"] == "PEM[Oxidation]TIDEK"
            assert "/" not in row["sequence"]
            assert row["mzpaf"] == frag.to_mzpaf()
            assert row["deltas"] == "" and row["isotopes"] == ""

    def test_losses_and_isotopes(self):
        frags = pt.fragment("PEPTIDE", ion_types=("y",), charges=(1,), neutral_deltas=["H2O"], isotopes=(0, 2, {"15N": 1}))
        rows = pt.fragment_records(frags)
        assert {r["isotopes"] for r in rows} == {"", "13C^2", "15N"}
        assert {r["deltas"] for r in rows} == {"", "H-2O-1"}
        one_c13 = pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("y",), charges=(1,), isotopes=(1,)))
        assert {r["isotopes"] for r in one_c13} == {"13C"}

    def test_mass_delta_mzpaf_is_signed_mass(self):
        rows = pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b",), charges=(1,), deltas=({-17.0: 2},)))
        assert {r["deltas"] for r in rows} == {"-17.0^2"}
        assert rows[2]["mzpaf"] == "b3{PEP}-34.0"

    def test_gain_is_negative_count(self):
        rows = pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b",), charges=(1,), deltas=({"H2O": -1},)))
        assert rows[2]["deltas"] == "H-2O-1^-1"
        assert rows[2]["mzpaf"] == "b3{PEP}+H2O"

    def test_mzpaf_none_when_unwritable(self, monkeypatch):
        def boom(self, include_sequence=True):
            raise pt.PeptacularError("no mzPAF")

        monkeypatch.setattr(pt.Fragment, "to_mzpaf", boom)
        (row,) = pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b",), charges=(1,))[:1])
        assert row["mzpaf"] is None

    def test_without_parent_sequence(self):
        frag = pt.parse("PEPTIDE").frag("b", 1, position=3, _include_sequence=False)
        (row,) = pt.fragment_records([frag])
        assert row["sequence"] is None and row["parent_sequence"] is None
        assert row["mzpaf"] == "b3"

    def test_empty_and_bad_input(self):
        assert pt.fragment_records([]) == []
        with pytest.raises(pt.PeptacularError, match="Fragment"):
            pt.fragment_records(["b1"])  # type: ignore[list-item]
        with pytest.raises(pt.PeptacularError, match="Fragment"):
            pt.fragment_records([1])  # type: ignore[list-item]
        with pytest.raises(pt.PeptacularError, match="inside a list"):
            pt.fragment_records([[1]])  # type: ignore[list-item]

    def test_wrong_top_level_type_is_type_error(self):
        # Same rule as pt.digest: a wrong top-level type is a TypeError.
        for bad in (None, 5):
            with pytest.raises(TypeError):
                pt.fragment_records(bad)  # type: ignore[arg-type]
            with pytest.raises(TypeError):
                pt.digest_records(bad, "trypsin")  # type: ignore[arg-type]
        with pytest.raises(TypeError):
            pt.digest_records([5], "trypsin")  # type: ignore[list-item]

    def test_mixed_sign_formula_delta_has_no_mzpaf(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"CH-2": 1})
        (row,) = pt.fragment_records([frag])
        assert row["deltas"] == "CH-2" and row["mzpaf"] is None

    def test_plain_formula_delta_is_a_gain(self):
        frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"C2H2O": 1})
        base = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3)
        (row,) = pt.fragment_records([frag])
        assert row["deltas"] == "C2H2O" and row["mzpaf"] == "b3{PEP}+C2H2O"
        assert frag.mass > base.mass

    def test_batch_nested_list_is_flattened(self):
        nested = pt.fragment(["PEPTIDE", "PEK/2"], ion_types=("b",), charges=(1,))
        rows = pt.fragment_records(nested)
        assert len(rows) == 7 + 3
        assert [r["parent_sequence"] for r in rows] == ["PEPTIDE"] * 7 + ["PEK"] * 3
        assert rows[:7] == pt.fragment_records(nested[0])

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
    def test_pandas_dtypes(self):
        import pandas as pd

        digest = pd.DataFrame(pt.digest_records([Entry("MKVLATSAGERTIDEK", "P1"), "PEPTIDEK"], "trypsin"))
        assert list(digest.columns) == list(pt.DIGEST_RECORD_KEYS)
        assert digest["start"].dtype == "int64" and digest["semi"].dtype == "bool"
        frags = pd.DataFrame(pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b", "by", "p"), charges=(1, 2))))
        assert list(frags.columns) == list(pt.FRAGMENT_RECORD_KEYS)
        assert frags["mz"].dtype == "float64" and frags["charge_state"].dtype == "int64"
        # position holds ints and None (precursor), so pandas makes it a float column, never object.
        assert frags["position"].dtype == "float64" and frags["end_position"].dtype == "float64"
        internal = frags[frags["ion_type"] == "by"]
        assert (internal["end_position"] >= internal["position"]).all()

    @pytest.mark.skipif(importlib.util.find_spec("polars") is None, reason="polars not installed")
    def test_polars_dtypes(self):
        import polars as pl

        digest = pl.DataFrame(pt.digest_records([Entry("MKVLATSAGERTIDEK", "P1"), "PEPTIDEK"], "trypsin"))
        assert digest.columns == list(pt.DIGEST_RECORD_KEYS)
        assert digest.schema["accession"] == pl.String and digest.schema["semi"] == pl.Boolean
        frags = pl.DataFrame(pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("b", "by", "p"), charges=(1, 2))))
        assert frags.columns == list(pt.FRAGMENT_RECORD_KEYS)
        assert frags.schema["position"] == pl.Int64 and frags.schema["end_position"] == pl.Int64
        assert frags.schema["mz"] == pl.Float64 and frags.schema["deltas"] == pl.String
        assert frags.filter(pl.col("ion_type") == "by")["end_position"].null_count() == 0
