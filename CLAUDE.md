# Peptacular — Claude Code Guide

## Project Overview

Peptacular is a ProForma 2.1 compliant Python library (v3.0.0) for peptide sequence parsing and annotation. It is in active JOSS submission. The main public API is `peptacular` (imported as `pt`), built on top of the `tacular` package.

Key entry points:
- `pt.parse(seq)` → returns a `ProFormaAnnotation` object (OOP API)
- `pt.mass(seq_or_list)` → functional API, supports batch/parallel processing
- `ProFormaAnnotation` in `src/peptacular/annotation/annotation.py` — the central class (~2000+ lines)

## Commands

```bash
just install        # uv sync --all-extras (dev + all extras)
just test           # pytest tests/
just test-cov       # pytest with branch coverage report
just lint           # ruff check src/
just format         # ruff isort + F401 + format
just ty             # ty type check src/
just check          # lint + test + ty
just docs           # sphinx-build docs/
just examples       # run all examples/*.py scripts
just paper          # build JOSS paper PDF via Docker (openjournals/inara)
just clean          # remove __pycache__, .pytest_cache, build artifacts
```

## Architecture

```
src/peptacular/
  annotation/       # ProFormaAnnotation class + modification helpers
  sequence/         # Sequence manipulation, parallel processing
  digestion/        # Enzyme digestion logic
  proforma_components/  # ProForma 2.1 data model primitives
  property/         # Mass/mz calculation, isotope properties
  fasta.py          # FASTA parsing
  chem.py           # Chemical formula utilities
  isotope.py        # Isotopic distribution
  constants.py      # Amino acid/element constants
```

## Docstring Style

Use Sphinx style consistently (matches existing docs and RTD setup):

```python
def foo(seq: str, inplace: bool = False) -> "ProFormaAnnotation":
    """One-line summary.

    :param seq: Description of seq.
    :type seq: str
    :param inplace: If True, modifies in place; if False, returns a new copy.
    :type inplace: bool
    :return: The modified annotation.
    :rtype: ProFormaAnnotation
    :raises ValueError: If seq is invalid.
    """
```

## Key Design Patterns

- **Factory / method chaining**: `set_*` methods accept `inplace: bool = False`. When `False` they return a new `ProFormaAnnotation` copy, enabling chains like `pt.parse("PEM[Oxidation]TIDE").set_charge(2).serialize()`.
- **Batch parallelism**: Functional API functions (e.g. `pt.mass([...])`) auto-parallelize for list inputs via `parallel.py`.
- **ProForma 2.1**: All serialization/deserialization targets the ProForma 2.1 spec.

## Tooling

- **Python**: ≥3.12
- **Package manager**: `uv`
- **Linter/formatter**: `ruff` (line length 160, targets py312)
- **Type checker**: `ty`
- **Test runner**: `pytest` with `pytest-cov`
- **Docs**: Sphinx + RTD theme + `sphinx-autodoc-typehints`
- **Build**: `hatchling`; version sourced from `src/peptacular/__init__.py`

## Testing

- Tests live in `tests/`
- Run `just test-cov` to check coverage; target is ≥83%
- Use `tmp_path` fixture (not `tempfile`) for file-based tests
- Do not mock the database or internal state — integration tests should exercise real code paths

## Linting Notes

- `__init__.py` files ignore `F403`/`F401` (star imports are intentional for the public API)
- `data.py` and `constants.py` ignore `E741` (ambiguous variable names are unavoidable for amino acid tables)
