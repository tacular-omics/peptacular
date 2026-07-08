"""Coverage for peptacular.sequence.parallel: start-method helpers, GIL detection,
and the parallel_apply_internal dispatcher.
"""

import multiprocessing as mp

from peptacular.sequence import parallel as parallel_mod
from peptacular.sequence.parallel import (
    _is_gil_disabled,
    get_available_start_methods,
    get_start_method,
    parallel_apply_internal,
    set_start_method,
)


def _double(x: int) -> int:
    return x * 2


class TestSetStartMethod:
    def test_none_is_noop(self):
        set_start_method(None)

    def test_valid_method(self):
        set_start_method(mp.get_start_method())

    def test_runtime_error_is_caught_and_warned(self, monkeypatch, caplog):
        def _raise(*args, **kwargs):
            raise RuntimeError("context already set")

        monkeypatch.setattr(parallel_mod.mp, "set_start_method", _raise)
        with caplog.at_level("WARNING"):
            set_start_method("fork")
        assert "Could not set start method" in caplog.text


class TestGetStartMethod:
    def test_returns_string(self):
        assert get_start_method() in get_available_start_methods()


class TestGetAvailableStartMethods:
    def test_returns_list(self):
        methods = get_available_start_methods()
        assert isinstance(methods, list)
        assert len(methods) > 0


class TestIsGilDisabled:
    def test_pre_313_is_false(self, monkeypatch):
        monkeypatch.setattr(parallel_mod.sys, "version_info", (3, 12, 0, "final", 0))
        assert _is_gil_disabled() is False

    def test_313_with_gil_enabled_attr_true(self, monkeypatch):
        monkeypatch.setattr(parallel_mod.sys, "version_info", (3, 13, 0, "final", 0))
        monkeypatch.setattr(parallel_mod.sys, "_is_gil_enabled", lambda: True, raising=False)
        assert _is_gil_disabled() is False

    def test_313_with_gil_enabled_attr_false(self, monkeypatch):
        monkeypatch.setattr(parallel_mod.sys, "version_info", (3, 13, 0, "final", 0))
        monkeypatch.setattr(parallel_mod.sys, "_is_gil_enabled", lambda: False, raising=False)
        assert _is_gil_disabled() is True

    def test_313_without_gil_enabled_attr(self, monkeypatch):
        monkeypatch.setattr(parallel_mod.sys, "version_info", (3, 13, 0, "final", 0))
        monkeypatch.delattr(parallel_mod.sys, "_is_gil_enabled", raising=False)
        assert _is_gil_disabled() is False

    def test_313_attribute_error_is_caught(self, monkeypatch):
        def _raise():
            raise AttributeError("boom")

        monkeypatch.setattr(parallel_mod.sys, "version_info", (3, 13, 0, "final", 0))
        monkeypatch.setattr(parallel_mod.sys, "_is_gil_enabled", _raise, raising=False)
        assert _is_gil_disabled() is False


class TestParallelApplyInternal:
    def test_empty_input(self):
        assert parallel_apply_internal(_double, []) == []

    def test_sequential(self):
        assert parallel_apply_internal(_double, [1, 2, 3], method="sequential") == [2, 4, 6]

    def test_process_pool(self):
        result = parallel_apply_internal(_double, [1, 2, 3, 4], method="process", n_workers=2)
        assert result == [2, 4, 6, 8]

    def test_thread_pool(self):
        result = parallel_apply_internal(_double, [1, 2, 3, 4], method="thread", n_workers=2)
        assert result == [2, 4, 6, 8]

    def test_verbose_logs_worker_and_chunksize(self, caplog):
        with caplog.at_level("DEBUG", logger=parallel_mod.logger.name):
            parallel_apply_internal(_double, [1, 2, 3, 4], method="thread", n_workers=2, verbose=True)
        assert "n_workers" in caplog.text
        assert "chunksize" in caplog.text

    def test_auto_method_none(self):
        result = parallel_apply_internal(_double, [1, 2, 3], method=None)
        assert result == [2, 4, 6]

    def test_explicit_chunksize(self):
        result = parallel_apply_internal(_double, [1, 2, 3, 4], method="sequential", chunksize=2)
        assert result == [2, 4, 6, 8]

    def test_func_kwargs_forwarded(self):
        def _add(x: int, offset: int) -> int:
            return x + offset

        result = parallel_apply_internal(_add, [1, 2, 3], method="sequential", offset=10)
        assert result == [11, 12, 13]
