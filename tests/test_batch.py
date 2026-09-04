"""Batch ordering, bounded consumption, diagnostics, and process portability."""

import multiprocessing as mp
from contextlib import closing

import pytest

import peptacular as pt
from peptacular.sequence.parallel import parallel_apply_internal


@pytest.mark.parametrize("method", ["sequential", "thread", "process"])
def test_collect_results_are_ordered_across_batches(method):
    inputs = ["PEPTIDE", "PEP[NoSuchModification]TIDE", "PEPTIDES", "PEP[", "MKR"]
    results = pt.batch("mass", iter(inputs), errors="collect", batch_size=2, method=method, n_workers=2)
    assert [r.index for r in results] == list(range(5))
    assert [r.input for r in results] == inputs
    assert [r.ok for r in results] == [True, False, True, False, True]
    assert results[0].value == pytest.approx(pt.mass(inputs[0]))
    assert results[1].error.code == "unresolved_modification"
    assert results[3].error.code == "invalid_notation"
    assert results[4].value == pytest.approx(pt.mass("MKR"))


@pytest.mark.parametrize("method", ["spawn", "forkserver"])
def test_process_context(method):
    if method not in mp.get_all_start_methods():
        pytest.skip("Start method unavailable on this platform")
    before = mp.get_start_method()
    results = pt.batch("mass", ["PEPTIDE", "MKR"], method="process", start_method=method, n_workers=2)
    assert all(r.ok for r in results)
    assert mp.get_start_method() == before


def test_lazy_bounded_consumption_and_close():
    seen = []

    def inputs():
        for i in range(100):
            seen.append(i)
            yield "PEPTIDE"

    with closing(pt.iter_batch("mass", inputs(), batch_size=3, method="thread", n_workers=2)) as results:
        assert next(results).index == 0
        assert seen == [0, 1, 2]
    assert seen == [0, 1, 2]


def test_default_errors_raise_and_iterator_errors_are_not_collected():
    with pytest.raises(ValueError):
        pt.batch("mass", ["PEPTIDE", "PEP[NoSuchModification]TIDE"])

    def inputs():
        yield "PEPTIDE"
        raise OSError("input source failed")

    results = pt.iter_batch("mass", inputs(), batch_size=1, errors="collect")
    assert next(results).ok
    with pytest.raises(OSError, match="input source failed"):
        next(results)


@pytest.mark.parametrize("operation,kwargs", [("wrong", {}), ("mass", {"typo": 1}), ("digest", {})])
def test_configuration_errors_raise_even_for_empty_input(operation, kwargs):
    with pytest.raises((TypeError, ValueError)):
        pt.batch(operation, [], errors="collect", **kwargs)


@pytest.mark.parametrize("key", ["n_workers", "chunksize", "batch_size"])
@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_bad_execution_settings(key, value):
    with pytest.raises(ValueError, match=key):
        pt.batch("mass", [], **{key: value})


@pytest.mark.parametrize("key", ["n_workers", "chunksize"])
@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_existing_parallel_api_validates_settings(key, value):
    with pytest.raises(ValueError, match=key):
        parallel_apply_internal(abs, [], **{key: value})


def test_small_auto_batches_allow_local_functions():
    def double(value):
        return value * 2

    assert parallel_apply_internal(double, [1, 2]) == [2, 4]


@pytest.mark.parametrize(
    "sequence,operation,kwargs,code,stage",
    [
        ("PEP[", "mass", {}, "invalid_notation", "parse"),
        ("PEP[NoSuchModification]TIDE", "mass", {}, "unresolved_modification", "calculate"),
        ("PEP[+42]TIDE", "comp", {}, "unavailable_composition", "calculate"),
        ("PEPTIDE", "mass", {"isotopes": 1000}, "invalid_adjustment", "calculate"),
        ("PEPTIDE", "fast_fragment", {"ion_types": ["by"]}, "unsupported_operation", "calculate"),
    ],
)
def test_diagnostics(sequence, operation, kwargs, code, stage):
    diagnostic = pt.diagnose(sequence, operation, **kwargs)
    assert diagnostic.code == code
    assert diagnostic.stage == stage
    assert diagnostic.message
    assert diagnostic.exception_type


def test_diagnose_success_and_annotation_not_mutated():
    annotation = pt.parse("PEPTIDE/2")
    assert pt.diagnose(annotation, "mass", charge=1) is None
    assert annotation.serialize() == "PEPTIDE/2"
    assert pt.diagnose("PEP[+42]TIDE", "mass") is None


def test_batch_parse_and_other_operations():
    assert pt.batch("parse", ["PEPTIDE"])[0].value.serialize() == "PEPTIDE"
    assert pt.batch("comp", ["PEPTIDE"])[0].value == pt.comp("PEPTIDE")
    assert pt.batch("mz", ["PEPTIDE"], charge=2)[0].value == pt.mz("PEPTIDE", charge=2)
    assert pt.batch("fragment", ["PEPTIDE"], ion_types=["b"], charges=[1])[0].value
    assert pt.batch("digest", ["PEPTIDEK"], enzyme=pt.Proteases.TRYPSIN)[0].value


def test_collect_does_not_catch_unexpected_type_errors():
    with pytest.raises(TypeError):
        pt.batch("fragment", ["PEPTIDE"], errors="collect", ion_types=object())


class CrashingAnnotation(pt.ProFormaAnnotation):
    def copy(self):
        import os

        if mp.current_process().name == "MainProcess":
            raise AssertionError("This input must run only in a worker")
        os._exit(17)


def test_worker_failure_propagates_in_collect_mode():
    from concurrent.futures.process import BrokenProcessPool

    with pytest.raises(BrokenProcessPool):
        pt.batch("mass", [CrashingAnnotation("PEPTIDE")], errors="collect", method="process", n_workers=1)


def test_parse_defaults_preserve_unresolved_modifications():
    sequence = "PEP[UnresolvedModification]TIDE"
    assert pt.diagnose(sequence, "parse") is None
    assert pt.batch("parse", [sequence])[0].value.serialize() == sequence
    assert pt.diagnose(sequence, "parse", validate=True).code == "invalid_annotation"


def test_invalid_residue_reports_validation_stage():
    error = pt.diagnose("PEPTIDE1")
    assert error is not None
    assert error.code in {"invalid_notation", "invalid_annotation"}


def test_invalid_input_and_execution_modes():
    assert pt.batch("mass", [123], errors="collect")[0].error.code == "invalid_input"
    with pytest.raises(TypeError):
        pt.batch("mass", [123])
    with pytest.raises(TypeError):
        pt.batch("mass", "PEPTIDE")
    with pytest.raises(ValueError):
        pt.batch("mass", [], errors="ignore")
    with pytest.raises(ValueError):
        pt.batch("mass", [], method="invalid")
    with pytest.raises(ValueError):
        pt.batch("mass", [], method="thread", start_method="spawn")
