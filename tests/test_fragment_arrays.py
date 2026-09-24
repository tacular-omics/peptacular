"""fragment_arrays() must return exactly the ions of fragment(), as numpy columns."""

import sys

import pytest

import peptacular as pt
from peptacular.interop import MissingOptionalDependencyError

np = pytest.importorskip("numpy")

SEQUENCES = [
    "PEPTIDE",
    "PEM[Oxidation]TIDEK/2",
    "[Acetyl]-PEPS[Phospho]TIDE-[Amidated]",
    "ACDEFGHIKLMNPQRSTVWY",
    "PEPTIDE/[Na:z+1,H:z+1]",
    "K",
    "PEPTIDE/-2",
    "PEPT[+79.966]IDE",
    "PEP[Formula:C2H2O]TIDE",
    "<13C>PEPTIDE",
    "<[Carbamidomethyl]@C>PEPCTIDEC",
    "{Glycan:Hex}PEPT[Phospho]IDE",
    "PEPTIDE/3",
    "PEPS[Phospho|+79.966]TIDE",
]


def _expected_rows(fragments_per_peptide):
    rows = []
    for index, fragments in enumerate(fragments_per_peptide):
        for fragment, record in zip(fragments, pt.fragment_records(fragments), strict=True):
            delta_mass = sum(
                (key.get_mass(monoisotopic=fragment.monoisotopic) if hasattr(key, "get_mass") else key) * count for key, count in fragment.deltas.items()
            )
            raw = fragment._isotopes
            isotope = raw if isinstance(raw, int) else int((raw or {}).get("13C", 0))
            rows.append(
                {
                    "peptide_index": index,
                    "ion_type": record["ion_type"],
                    "position": record["position"] or 0,
                    "end_position": record["end_position"] or 0,
                    "charge_state": record["charge_state"],
                    "mz": record["mz"],
                    "mass": record["mass"],
                    "isotope": isotope,
                    "isotope_label": record["isotopes"],
                    "delta_label": record["deltas"],
                    "delta_mass": delta_mass,
                }
            )
    return rows


def _assert_matches(columns, fragments_per_peptide):
    assert tuple(columns) == pt.FRAGMENT_ARRAY_KEYS
    expected = _expected_rows(fragments_per_peptide)
    lengths = {len(column) for column in columns.values()}
    assert lengths == {len(expected)}
    for key in pt.FRAGMENT_ARRAY_KEYS:
        # Exact equality, floats included: the arrays must be the same numbers as fragment().
        assert columns[key].tolist() == [row[key] for row in expected], key


CASES = [
    {"ion_types": ["b", "y"], "charges": [1, 2]},
    {"ion_types": ["a", "b", "c", "x", "y", "z"], "charges": [1, 2, 3]},
    {"ion_types": ["b", "y"], "charges": None},
    {"ion_types": ["b", "y", "p", "by", "i"], "charges": [1]},
    {"ion_types": ["b", "y"], "charges": [1, 2], "isotopes": [0, 1, 2]},
    {"ion_types": ["b", "y"], "charges": [1], "isotopes": [0, {"15N": 1}]},
    {"ion_types": ["b", "y"], "charges": [1], "deltas": [None, -18.010565, "H2O", -17.026549]},
    {"ion_types": ["b", "y"], "charges": [1], "deltas": [None, -18.010565], "max_ndeltas": 2},
    {"ion_types": ["b", "y"], "charges": [1, 2], "neutral_deltas": ["H2O", "NH3"], "max_ndeltas": 2},
    {"ion_types": ["b", "y"], "charges": [1, "Na:z+1", -1]},
    {"ion_types": ["b", "y"], "charges": [1, 2], "monoisotopic": False},
    {"ion_types": ["c-H", "z+H", "z."], "charges": [1, 2]},
]


@pytest.mark.parametrize("kwargs", CASES)
def test_matches_fragment(kwargs):
    columns = pt.fragment_arrays(SEQUENCES, **kwargs)
    _assert_matches(columns, pt.fragment(SEQUENCES, **kwargs))


def test_matches_fragment_with_composition():
    sequences = [s for s in SEQUENCES if "+79" not in s]
    kwargs = {"ion_types": ["b", "y"], "charges": [1, 2], "calculate_with_composition": True}
    _assert_matches(pt.fragment_arrays(sequences, **kwargs), pt.fragment(sequences, **kwargs))


def test_single_sequence_and_annotation_method():
    columns = pt.fragment_arrays("PEM[Oxidation]TIDE/2", ion_types=["b", "y"])
    _assert_matches(columns, [pt.fragment("PEM[Oxidation]TIDE/2", ion_types=["b", "y"])])
    assert set(columns["peptide_index"].tolist()) == {0}

    annotation = pt.parse("PEM[Oxidation]TIDE/2")
    kwargs = {"ion_types": ["b", "y", "by"], "charges": [1, 2], "min_length": 2, "max_length": 4}
    _assert_matches(annotation.fragment_arrays(**kwargs), [annotation.fragment(**kwargs)])


def test_dtypes_and_order():
    columns = pt.fragment_arrays(["PEPTIDE", "PEPTIDEK"], ion_types=["b"], charges=[1])
    for key in ("peptide_index", "position", "end_position", "charge_state", "isotope"):
        assert columns[key].dtype == np.int64, key
    for key in ("mz", "mass", "delta_mass"):
        assert columns[key].dtype == np.float64, key
    for key in ("ion_type", "isotope_label", "delta_label"):
        assert columns[key].dtype == object, key
    assert columns["peptide_index"].tolist() == [0] * 7 + [1] * 8
    assert columns["position"].tolist() == [*range(1, 8), *range(1, 9)]


def test_empty_inputs():
    columns = pt.fragment_arrays([], ion_types=["b", "y"], charges=[1])
    assert tuple(columns) == pt.FRAGMENT_ARRAY_KEYS
    assert all(len(column) == 0 for column in columns.values())


def test_parallel_matches_sequential():
    sequences = SEQUENCES * 5
    kwargs = {"ion_types": ["b", "y"], "charges": [1, 2]}
    sequential = pt.fragment_arrays(sequences, **kwargs)
    threaded = pt.fragment_arrays(sequences, method="thread", n_workers=3, chunksize=4, **kwargs)
    for key in pt.FRAGMENT_ARRAY_KEYS:
        assert threaded[key].tolist() == sequential[key].tolist(), key


def test_invalid_input_raises_like_fragment():
    with pytest.raises(pt.PeptacularError):
        pt.fragment("PEPTIDE", ion_types=["b"], charges=[1], deltas=[-10000.0])
    with pytest.raises(pt.PeptacularError):
        pt.fragment_arrays("PEPTIDE", ion_types=["b"], charges=[1], deltas=[-10000.0])


def test_polars_accepts_columns():
    pl = pytest.importorskip("polars")
    columns = pt.fragment_arrays(SEQUENCES, ion_types=["b", "y"], charges=[1, 2])
    frame = pl.DataFrame(columns)
    assert frame.columns == list(pt.FRAGMENT_ARRAY_KEYS)
    assert frame.height == len(columns["mz"])
    assert frame["mz"].to_list() == columns["mz"].tolist()
    assert frame["ion_type"].dtype == pl.String


def test_pyarrow_accepts_columns():
    pa = pytest.importorskip("pyarrow")
    columns = pt.fragment_arrays(SEQUENCES, ion_types=["b", "y"], charges=[1, 2])
    table = pa.table(columns)
    assert table.column_names == list(pt.FRAGMENT_ARRAY_KEYS)
    assert table.num_rows == len(columns["mz"])
    assert table["delta_label"].type == pa.string()


def test_missing_numpy_says_how_to_install(monkeypatch):
    monkeypatch.setitem(sys.modules, "numpy", None)
    with pytest.raises(MissingOptionalDependencyError, match=r'pip install "peptacular\[numpy\]"'):
        pt.fragment_arrays("PEPTIDE")
    with pytest.raises(MissingOptionalDependencyError):
        pt.parse("PEPTIDE").fragment_arrays()
