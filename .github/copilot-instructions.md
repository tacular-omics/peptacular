# Copilot instructions for peptacular

Read [`CLAUDE.md`](../CLAUDE.md) at the repo root. It is the canonical guide: commands,
architecture, public API, conventions and gotchas. Library usage is in
[`llms-full.txt`](../llms-full.txt).

Key rules:

1. Target ProForma 2.1. A new operation gets a `ProFormaAnnotation` method and a
   functional wrapper in `src/peptacular/sequence/` that accepts
   `str | ProFormaAnnotation | list` and the `n_workers` / `chunksize` / `method` kwargs.
2. Sphinx-style docstrings (`:param:` / `:return:`), not Google style. Type
   annotations on all public functions. Python >= 3.12, ruff line length 160.
3. Mutators take `inplace: bool = True` and return the annotation so calls chain.
4. Tests go in `tests/`. Use `pytest.approx` and `tmp_path`, and do not mock
   internals. Branch coverage must stay at or above 79% (`just test-cov`).
5. Before a commit: `uv run ruff check src tests`, `uv run ruff format --check src tests`,
   `uv run ty check src`, `uv run pytest tests`.
6. Never bump the version, tag, or publish. Only the tacular-omics overseer releases.
