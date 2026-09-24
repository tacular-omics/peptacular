# peptacular — Claude Code Guide

## Project overview

peptacular is a ProForma 2.1 compliant Python library for peptide sequences: parse and
serialize ProForma, edit modifications, and calculate mass, m/z, elemental composition,
fragment ions, isotopic envelopes, enzymatic digests and physicochemical properties. It
is imported as `pt` and is under JOSS review (`paper/`).

Place in the tacular-omics graph: tier 1. It depends on `tacular>=1.2,<2` (all
modification, amino-acid, element, ion-type, protease and neutral-loss data) and is used
by `paftacular` (optional extra), `spxtacular`, `peff_digest` and `pepbit`. A breaking
change here must be noted for those. `__init__` imports every public name explicitly and
`pt.__all__` is the public API (a test checks `dir(pt)` against it). Of tacular only the
enums `IonType`, `NeutralDelta` and `Proteases` are re-exported; import lookups such as
`UNIMOD_LOOKUP` or `PROTEASE_LOOKUP` from `tacular` directly.

Key entry points:
- `pt.parse(seq)` returns a `ProFormaAnnotation` (OOP API).
- `pt.mass(seq_or_list)`, `pt.fragment(...)`, `pt.digest(...)` and friends are the
  functional API (`sequence/`): they take a string, an annotation or a list of either,
  and parallelize lists.
- `ProFormaAnnotation` lives in `src/peptacular/annotation/annotation.py` (~4900 lines).

## Commands

```bash
just install        # uv sync --all-extras (dev + all extras)
just test           # uv run pytest tests/          (~2860 tests, ~15 s)
just test-cov       # pytest with branch coverage + scripts/check_branch_coverage.py (min 79%)
just lint           # ruff check src/ tests/ (same as CI)
just format-check   # ruff format --check src/ tests/ (same as CI)
just format         # ruff isort + F401 fix + format on src/tests/profile  -- WRITES FILES
just ty             # ty check src/
just check          # lint + format-check + test + ty
just docs           # sphinx-build -W docs -> docs/_build/html
just docs-test      # sphinx doctest build (the rst testcode blocks; 16 tests)
just examples       # run every examples/*.py (CI runs this)
just paper          # build the JOSS paper PDF via Docker (openjournals/inara)
just check-version  # scripts/release_version.py check
just clean          # remove caches and build artifacts
```

Also: `install-prod`, `sync`, `codecov-tests`, `docs-clean`, `docs-open`, `upgrade`
(pyupgrade --py312-plus), `all` (clean + install + test), and the release recipes below.

Docstring examples (`>>>`) in `src/` run via `tests/test_doctest.py` as part of pytest.

CLI: the `peptacular-mcp` console script (extra `mcp`) is a stdio MCP server with 12
tools. `peptacular-mcp --check` prints the installed tools, versions and request limits
as JSON. Setup is in `docs/mcp.rst`.

CI (`.github/workflows/ci.yml`): `ruff check src tests`, `ruff format --check src tests`,
`ty check src`, `release_version.py check`, pytest on 3.12-3.14 / macOS / Windows,
lowest-direct resolution, the built wheel per extra (`scripts/test_wheel.py`), then
`just test-cov`, `just docs`, `just docs-test`, `just examples`.

## Architecture

```
src/peptacular/
  __init__.py            star-imports everything below plus `from tacular import *`
  annotation/            the OOP API
    annotation.py        ProFormaAnnotation: storage, get/set/append/extend/pop/remove/clear
                         per mod type, mass/mz/comp, fragment, digest, slicing, conversions
    parser.py            ProForma string -> ProFormaAnnotation (syntax only)
    serializer.py        ProFormaAnnotation -> ProForma string
    mod.py               Mod (value + count) and Mods (one mod type's collection), Interval
    cached_comps.py      lru_cached isotope / delta / charge-carrier resolution
    frag.py              Fragment dataclass (.mz, .to_mzpaf(), .composition)
    slicing.py, manipulation.py, combinatorics.py, ambiguity.py, randomizer.py,
    mod_builder.py, positions.py, utils.py
  sequence/              the functional API: one wrapper per operation that accepts
                         str | annotation | list and dispatches through parallel.py
    parallel.py          AUTO_PARALLEL_MIN_ITEMS = 1000; below that lists run sequentially
    basic.py, mass_funcs.py, fragmentation.py, digestion.py, isotope.py, properties.py,
    mod_builder.py, combinatoric.py, subseqs.py, transformations.py, converters.py
  proforma_components/   ProForma 2.1 data model (comps.py), string parsers (lru_cached)
                         and serializers for mod tags, formulas, glycans, charge carriers
  digestion/             EnzymeConfig, DigestProtocol, cleavage-site logic
  property/              amino-acid property scales (data.py) and calculation (core.py,
                         prop.py: AnnotationProperties behind `annot.prop`)
  spans.py               Span(start, end, missed_cleavages) and span builders
  isotope.py             BRAIN isotope envelopes, averagine estimates, IsotopicData
  chem.py                chem_mass / chem_comp / chem_formula / parse_formula
  fasta.py               parse_fasta, iter_fasta (streaming, .gz), FastaSequence
  batch.py               batch / iter_batch / diagnose with per-item error collection
  diagnostics.py         Diagnostic, UnknownModificationError, CompositionError, ...
  proforma_json.py       versioned lossless JSON; schema in schemas/proforma-json-v1.schema.json
  constants.py           ModType, ParallelMethod, PROTON/ELECTRON/NEUTRON masses
  _regex_utils.py (private), utils.py
  interop/               optional pyteomics / psm_utils / alphabase converters (lazy imports)
  mcp/                   optional MCP server (cli.py, server.py, operations.py, contracts.py)
```

Data flow: `parse()` checks syntax only and stores each modification as an interned
string with a count in a `Mods` collection (`internal_mods` is `{index: Mods}`). Nothing is
resolved against tacular until you ask for mass, m/z or composition. Then
`proforma_components.parsers` turns each mod string into a component (`TagName`,
`TagAccession`, `TagMass`, `ChargedFormula`, `GlycanTag`, ...) that looks up tacular. These
parsers are lru_cached, so repeated mods are cheap. Plain precursor/neutral masses take a
fast path in `ProFormaAnnotation._mass_and_charge`. Everything else goes through `frag()`,
which builds a `Fragment`.

## Public API (everything is on `pt`)

- **Parse / serialize**: `parse`, `parse_chimeric`, `serialize`, `serialize_chimeric`,
  `validate`, `ProFormaAnnotation`, `join`, `split`.
- **Mass**: `mass`, `mz`, `comp`, `chem_mass`, `chem_comp`, `chem_formula`,
  `parse_formula`, `add_composition`, `merge_compositions`.
- **Fragments**: `fragment`, `frag`, `fast_fragment`, `Fragment`.
- **Digestion**: `digest`, `simple_digest`, `semi_digest`, `left_semi_digest`,
  `right_semi_digest`, `nonspecific_digest`, `cleavage_sites`, `simple_cleavage_sites`,
  `EnzymeConfig`, `DigestProtocol`. Spans: `Span`, `build_spans`,
  `build_enzymatic_spans`, `build_semi_spans`, `build_left_semi_spans`,
  `build_right_semi_spans`, `build_non_enzymatic_spans`, `calculate_span_coverage`,
  `span_to_sequence`.
- **Isotopes**: `isotopic_distribution`, `brain_isotopic_distribution`,
  `estimate_isotopic_distribution`, `merge_isotopic_distributions`, `averagine_comp`,
  `estimate_averagine_comp`, `IsotopicData`.
- **Modifications**: `modify` (static + variable combinatorics), `get_mods`, `set_mods`,
  `append_mods`, `extend_mods`, `pop_mods`, `remove_mods`, `filter_mods`, `strip_mods`,
  `condense_static_mods`, `condense_to_peptidoform`, `is_modified`, `get_mod_type`,
  `Mod`, `Mods`, `Interval`, `ModType`.
- **Sequence utilities**: `sequence_length`, `count_residues`, `percent_residues`,
  `is_ambiguous`, `reverse`, `shift`, `shuffle`, `sort`, `permutations`,
  `combinations`, `combinations_with_replacement`, `product`, `coverage`,
  `percent_coverage`, `modification_coverage`, `find_subsequence_indices`,
  `is_subsequence`, `annotate_ambiguity`, `generate_random`.
- **Properties**: `hydrophobicity`, `pi`, `charge_at_ph`, `aromaticity`,
  `secondary_structure`, `calc_property`, `calc_window_property` and ~20 more named
  scales, plus the scale enums (`HydrophobicityScale`, `HPLCScale`, ...) and the
  `annot.prop` object (`AnnotationProperties`).
- **Converters**: `convert_ip2_sequence`, `convert_diann_sequence`,
  `convert_casanovo_sequence`, `to_ms2_pip`, `from_ms2_pip`.
- **FASTA / batch**: `parse_fasta`, `parse_fasta_text`, `iter_fasta`, `FastaSequence`,
  `batch`, `iter_batch`, `diagnose`, `BatchResult`, `Diagnostic`.
- **JSON**: `to_proforma_json`, `from_proforma_json`, `to_proforma_dict`,
  `from_proforma_dict`, `get_proforma_json_schema`.
- **ProForma components** (`proforma_components`): `ChargedFormula`, `FormulaElement`,
  `GlycanTag`, `GlycanComponent`, `ModificationTags`, `TagName`, `TagAccession`,
  `TagMass`, `TagCustom`, `TagInfo`, `GlobalChargeCarrier`, `IsotopeReplacement`,
  `FixedModification`, `PositionRule`, `Peptidoform`, `PeptidoformIon`,
  `CompoundPeptidoformIon`, `CV`.
- **Errors**: `UnknownModificationError`, `CompositionError`,
  `InvalidAdjustmentError`, `UnsupportedOperationError` (all `ValueError`s).
- **Parallel**: `set_start_method`, `get_start_method`, `get_available_start_methods`,
  `ParallelMethod`.
- **Constants**: `PROTON_MASS`, `ELECTRON_MASS`, `NEUTRON_MASS`, `C13_NEUTRON_MASS`,
  `PEPTIDE_AVERAGINE_NEUTRON_MASS`, `AVERAGINE_RATIOS`, `PROFORMA_JSON_SCHEMA_ID`.
- **Optional** `peptacular.interop` (not star-imported): `to/from_pyteomics`,
  `to/from_psm_utils`, `to/from_alphabase_row`, `to/from_alphabase_dataframe`,
  `LossPolicy`.

## Conventions

- **Docstrings: Sphinx style** (`:param:` / `:type:` / `:return:` / `:rtype:` /
  `:raises:`). They match the RTD setup and most of the code. About 70 older functions
  still use Google `Args:` blocks. Convert them when you touch them; do not add new ones.

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
- **Method chaining**: every `set_*` / `append_*` / `extend_*` / `remove_*` method takes
  `inplace: bool = True` and returns the annotation (itself, or a copy when
  `inplace=False`), so `pt.parse("PEM[Oxidation]TIDE").set_charge(2).serialize()` works.
- **Functional API parity**: a new operation gets an annotation method and a
  `sequence/` wrapper that accepts `str | ProFormaAnnotation | Sequence[...]` and takes
  `n_workers`, `chunksize`, `method`.
- Type annotations on every public function. Python >= 3.12, ruff line length 160.
- `__init__.py` files ignore `F403`/`F401`, because star imports build the public API.
  `data.py` and `constants.py` ignore `E741`, because amino-acid tables need
  single-letter names.
- Tests: in `tests/`, use `pytest.approx` for floats, use the `tmp_path` fixture (not
  `tempfile`), and do not mock internal state: exercise real code paths. Branch coverage
  has a minimum of 79% and a target of 83% (`just test-cov`). The combined figure that
  coverage.py prints is a different metric.

## Gotchas

- **`parse()` checks syntax, not meaning.** `pt.parse("PEP[Foo]TIDE")` succeeds;
  `.mass()` then raises `UnknownModificationError`. Pass `validate=True` to `parse` (or
  call `pt.validate`) to resolve mods up front.
- **`mass()` includes the annotation's charge.** `pt.parse("PEPTIDE/2").mass()` is
  `M + 2 * proton` (801.37), not the neutral mass. Use `neutral_mass()` for the
  uncharged value, and `mz()` for m/z.
- **Delta-mass mods have no composition.** `pt.comp("PEP[+15.995]TIDE")` raises
  `CompositionError`, while `mass()` works. Isotope envelopes need a composition, so use
  `estimate_isotopic_distribution(mass)` for these.
- **`fast_fragment` returns a dict** `{(IonType, charge): [mz, ...]}`, not `Fragment`
  objects. Its values agree with `fragment()` to about 1e-8 Da, not bit for bit.
- **Two digest styles.** The functional `pt.digest(seq, enzyme=...)` (and every `pt.*digest`)
  returns `[(sequence, Span), ...]`. `ProFormaAnnotation` methods ending in `_spans`
  (`digest_spans`, `simple_digest_spans`, `sequential_digest_spans`, `semi_spans`, ...)
  yield `Span`s; slice the annotation with them (`annot[span]`).
- **`enzyme` is a protease name or a compiled pattern.** A `str` is only looked up in
  `PROTEASE_LOOKUP` (`"trypsin"`, `Proteases.TRYPSIN`); an unknown string raises
  `UnknownEnzymeError`. For a custom rule pass `re.compile("(?<=[KR])")`.
- **Lists under 1000 items run sequentially** unless you pass `n_workers` or `method`
  (`AUTO_PARALLEL_MIN_ITEMS`). Process pools use the platform default start method,
  which is `fork` on Linux before Python 3.14. With `fork`, a process that already runs
  threads gets a DeprecationWarning (the test run shows it). `set_start_method("spawn")`
  avoids the warning.
- **`pt.parse` of one string cannot hold chimeric or cross-linked input**
  (`PEPTIDE+ELVIS` raises). Use `parse_chimeric`.
- **Every public module has `__all__`** (a test enforces it). A new public name must be
  added to the module's `__all__` and, if it belongs on `pt`, to `peptacular/__init__.py`.
- **`ProFormaAnnotation` is unhashable** because it is mutable. Use `annot.serialize()`
  as a dict key or set member.
- **Errors**: input errors raise `PeptacularError` (a `ValueError`) or a subclass
  (`ProFormaFormatError`, `UnknownModificationError`, `UnknownEnzymeError`,
  `InvalidPositionError`, ...). Do not add a bare `raise ValueError`; a test fails on it.
  `n_workers`, `chunksize` and `method` are keyword-only on every function.
- `secondary_structure()` returns fractions (0-1) keyed by `SecondaryStructureType`, not
  percentages.
- `just format` rewrites files in place. Run it only on your own branch.

## Releasing

Only the tacular-omics overseer bumps versions or publishes. See `just --list`
(`set-version`, `sync-version`, `check-version`).

- Version source: `__version__` in `src/peptacular/__init__.py`
  (`[tool.hatch.version]`). `scripts/release_version.py` (synced from the workspace
  `templates/scripts/`, do not edit here) copies it to `CITATION.cff`, dates the
  `CHANGELOG.md` `## [Unreleased]` section, and rejects `version`/`grants` in
  `.zenodo.json`.
- Publishing is the `publish.yml` workflow on a GitHub release (PyPI trusted publishing).
- The JOSS paper is in `paper/` (`just paper`).

## Workspace note

This repo is also developed inside the tacular-omics uv workspace
(`~/Repos/tacular-omics/packages/peptacular`). There, `uv run` uses the shared `.venv`
and the root `uv.lock` (with the local `tacular` checkout), not this repo's own lock. See
the workspace CLAUDE.md.
