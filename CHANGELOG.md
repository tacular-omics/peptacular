
# Changelog

All notable changes to this project will be documented in this file.

## [3.3.0] (2026-09-04)

### Added
- Optional local MCP integration with 12 stateless tools for theoretical calculations, annotation transformations, and reference lookup. Small inline batches return bounded structured results with diagnostics.
- Optional Pyteomics, psm_utils, and AlphaBase interoperability adapters with separate installation extras.
- Versioned ProForma JSON serialization for annotations and structured components, with a bundled JSON Schema.
- Streaming `iter_fasta()` with gzip paths, encoding overrides, line-numbered structure errors, and caller-owned stream preservation.
- `batch()` and bounded `iter_batch()` APIs with ordered input indexes, opt-in error collection, reusable executors, and local multiprocessing contexts.
- `diagnose()` and structured diagnostics for parsing, validation, and calculation failures, with specific ValueError subclasses for common scientific errors.

### Fixed
- Materialized batch digestion spans so lazy failures are collected and process workers can serialize their results.
- Included intrinsic charge when applying electron-mass corrections to annotation isotope distributions.
- Hardened JSON field and number validation, expanded the bundled schema, and rejected duplicate JSON keys.
- Guarded optional conversions against annotation data loss, invalid element counts, and truncated AlphaBase charges.
- Precursor (`p`) and neutral (`n`) fragments generated through `fragment()` now leave their position unset, allowing the default lazy `.composition` and `.sequence` properties to represent the complete parent sequence instead of raising an invalid-position `ValueError`.
- Curly-brace glycan mass components now reject non-finite and malformed numeric values such as `{+nan}`, `{+inf}`, and `{1e309}` instead of propagating `NaN` or infinity through mass calculations.
- Programmatically constructed integer glycan masses are normalized to floats and serialize consistently instead of being mistaken for charged formulas.
- Corrected John R. Yates III's given name, surname, and suffix fields in the JATS paper metadata.
- Unified global isotope and elemental adjustment handling across mass and composition paths, including ion terminal atoms and impossible isotope or loss counts.
- Fast fragmentation now accounts for annotation isotope labels, intrinsic charges, labile precursor modifications, and average masses through regular-path fallbacks. Unsupported ion series and invalid charges now raise instead of returning misleading values.
- Derived ion mass offsets from elemental formulas so charge states share consistent mass precision.
- Forwarded `keep=False` correctly when filtering modifications on a copy.
- Protected cached delta and charge mappings against caller mutation, bounded helper caches, preserved negative charge-carrier atom counts, and rejected non-finite delta values and invalid counts.

### Changed
- Expanded package CI to cover pull requests and release-relevant configuration files, removed redundant test-suite executions, added package-build validation, and made the publishing workflow run lint, type, test, and build checks before uploading to PyPI.
- Small automatic batches run sequentially, worker counts are capped to available work, and worker/chunk settings are validated before execution.
- Added Python 3.13/3.14 and macOS/Windows CI coverage, documentation and example checks, an installed-wheel smoke check, and an explicit 79% branch coverage baseline with an 83% improvement target.
- Migrated development dependencies to standard dependency groups and resolved the new optional dependency extras.

### Deferred
- Cross-link calculations remain on the separate feature branch.
- Select random modifications from the reference databases.
- Review terminal modification handling for w, v, and d ions.
- Review string interning and modification cache reuse.

## [3.2.0]
### Added
- Glycan composition components given as a molecular formula or monoisotopic mass in curly braces, intermixed with named monosaccharides (ProForma 2.1 §10.2): `Glycan:{C8H13N1O5}1Hex2` (formula), `Glycan:{C8H13[15N1]O5}1Hex2` (isotope-labelled formula), `Glycan:{C8H13N1O5Na1:z+1}1Hex2` (charged formula, level 3), and `Glycan:{+203.079}1Hex2` (bare mass). Formula/charged-formula components contribute an elemental composition (and charge); a bare-mass component contributes mass and routes through the delta-mass path (so `mass()` works and `comp()` raises the same "cannot calculate composition with delta mass" error as any other bare-mass modification). `GlycanComponent` now accepts `Monosaccharide | ChargedFormula | float` — the previous `NotImplementedError` on a formula component is removed — and serializes formula components as `{...}` (no `Formula:` prefix) and masses as `{+m}`.

### Fixed
- `Fragment.composition` on the lazy path (`calculate_composition=False`) now applies the fragment's ion-type offset instead of returning the sub-sequence's *precursor* composition. A b-ion's lazily computed composition was heavier than its own `.mass` by a full water (e.g. `EVTKLE` `b4`: composition implied 476.27 Da vs the fragment's 458.26 Da); y-ions coincidentally matched because a y-ion's neutral composition equals its C-terminal sub-sequence's precursor. The lazy path now agrees element-for-element with the eager (`calculate_composition=True`) path across b/y/a/c/x/z ions, charges, mods, custom deltas and global isotope labels.
- Glycan compositions with whitespace between monosaccharide/count tokens (e.g. `Glycan:Hex5 HexNAc4`) are now parsed instead of raising `Could not parse glycan composition`. Per ProForma 2.1 §10.2 ("Spaces MAY be used around the names or numbers"), `_parse_glycan_composition` now skips whitespace between tokens and between a name and its count.

## [3.1.2]
### Changed
- Bumped the `tacular` dependency floor to `>=1.1.0`, which fixes several upstream data-consistency bugs: isotope-labelled modification compositions (e.g. `UNIMOD:536`, `Label:13C(2)15N(1)`) that dropped their isotope atoms while keeping the correct mass, all 9 internal-fragment-ion mass offsets (previously shifted so the default "by" internal fragment was `-CO` instead of `0`), and two neutral-loss formulas (Formic acid, Formamide) that were parsed with a dropped repeated-element count.

### Added
- Clearer, more actionable parse/validation error messages aimed at both humans and AI agents: they now name the offending value, state what was expected, and (for parse errors) point at the exact position. Unknown-modification errors include a hint listing valid ways to specify a modification (name, CV accession, formula, glycan, or delta mass).

### Fixed
- Repeated fragment charge carriers (e.g. `frag(..., charge=["Na:z+1", "Na:z+1"])`) no longer collapse to a single carrier: `adjust_mass_mz`/`adjust_comp` built the fragment's adduct tuple from `_mods.keys()`, dropping the per-carrier occurrence count, so the returned `Fragment` reported one Na (with a `neutral_mass` off by a full sodium) and serialized to `[M+Na]` where the mass reflected two. `Fragment.to_mzpaf` now also folds a carrier's collection-level count into the mzPAF repeat prefix, so a doubly-repeated carrier renders `[M+2Na]` consistently whether written as two list entries or as `Na:z+1^2`
- `Fragment.sequence` now serializes the fragment's external charge instead of `charge_state`, which includes charge intrinsic to an internal formula modification (e.g. `[Formula:CH2:z+1]`); a fragment with one external proton and an internal `+1` formula charge previously emitted `/2` (double-counting the internal charge) contradicting its own mass/mz and the `comp` property
- `SequenceElement.get_composition` now routes its modification merge through the negative-count-safe `add_composition` helper instead of `Counter` `+=`; a residue carrying an atom-removing mod (e.g. `N[Formula:O-1]`) previously kept the removed atom, giving a composition/mass too heavy by that atom (the parallel `comp()` path was already fixed)
- `GlobalChargeCarrier.to_mz_paf`/`Fragment.to_mzpaf` now include a repeated adduct's occurrence count (e.g. `[M+2Na]` for two sodium atoms, per mzPAF spec section 4.7's own example) instead of silently dropping it, and multiple distinct adducts are now alphabetized (`[M+H+Na]` rather than `[M+Na+H]`) as the spec recommends
- `Fragment.to_mzpaf` now raises a clear `ValueError` for a numeric neutral-loss/gain delta instead of emitting a bare mass (e.g. `y1{E}+15.99490`) that isn't valid mzPAF syntax; mzPAF neutral losses (spec section 4.5) must be a chemical formula or a named reference group, and there is no representation for an arbitrary unnamed mass delta
- `ChargedFormula.from_mz_paf` no longer requires a `Formula:` prefix (real mzPAF chemical formulas never carry one) and its negative-loss branch no longer relies on a "prepend a 0" trick that never actually parsed; `ChargedFormula.to_mz_paf`'s own output (e.g. `-H2O`) previously could not be round-tripped back through `from_mz_paf` at all
- `set_internal_mods`/`set_nterm_mods`/`set_cterm_mods`/etc. now strip whitespace from a modification string before storing/interning it; a padded value like `"  Oxidation  "` previously serialized to non-canonical ProForma (`PEP[  Oxidation  ]TIDE`) and, since it was interned as-is, prevented peptides using differently-padded-but-equivalent mod strings from sharing the same interned string and downstream parser cache entry
- Internal-fragment mzPAF labels (`_INTERNAL_MASS_DIFFS`) now match tacular>=1.1.0's corrected ion offsets; 5 of the 9 non-default internal ion types (`ax`, `bx`, `az`, `bz`, `cy`) previously carried a stale label whose implied mass no longer matched the actual computed fragment mass
- `GlycanComponent.get_composition` now multiplies by occurrence like `get_mass` does; a glycan monosaccharide count > 1 (e.g. `Glycan:Hex3`, or any real N-glycan) previously produced a composition whose mass disagreed with `mass()` by the count of dropped units
- Cached charge-carrier/delta/isotope-loss accessors (`ChargeCarrierInfo.composition`/`to_fragment_mapping`/`to_explicit_fragment_mapping`, `DeltaInfo.to_fragment_mapping`, `get_losses`) no longer return a shared mutable container from a process-wide `@cache` singleton; mutating a returned dict/Counter could previously corrupt every future caller with the same inputs
- `Fragment.neutral_mass` now undoes the per-charge electron-mass correction baked into `.mass`, so it is charge-invariant again (previously off by `charge * electron_mass`, with the sign flipping between positive and negative charge states)
- mzPAF serialization of negatively charged fragments no longer includes a minus sign (e.g. `y3{IDE}^1` instead of `y3{IDE}^-1`), per mzPAF spec section 4.8 ("the charge state component ... MUST NOT include the minus sign")
- `modify` now enumerates all positional isomers by default (previously the public API forced `unique_peptidoforms=True`, collapsing e.g. the 6 single-phospho placements on S/T/Y down to 2); the `unique_peptidoforms` flag is now exposed on the functional API
- `modify` now treats a bare string modification value as a single modification instead of shredding it into per-character mods (e.g. `{'S': 'Phospho'}` no longer becomes `S[P]`, `S[h]`, …)
- `condense_mods_to_intervals(inplace=False)` no longer mutates the original annotation; `copy()`/`update()` now deep-copy interval objects instead of sharing them
- Non-specific digestion now includes the full-length peptide, and a single-residue sequence now yields itself instead of an empty list (previously `build_non_enzymatic_spans` capped at length−1)
- Semi-enzymatic digestion no longer drops valid in-range peptides when `max_len` is smaller than a missed-cleavage parent span (the `max_len` filter was applied to parent spans before deriving semi sub-peptides)
- Charge carriers must now be a bare charged formula (e.g. `/[Na:z+1]`, `/[C2H6:z+2]`) per ProForma 2.1 §11.5; a `Formula:`/`Glycan:`-prefixed charge carrier (e.g. `/[Formula:C2H6:z+2]`) is now rejected with a clear message instead of parsing and then crashing with `Failed to parse element 'Fo'` on mass calculation. Corrected the `ChargedFormula`/`GlobalChargeCarrier` docstrings that advertised the invalid prefixed form.
- Unterminated modification brackets (e.g. `PEP[Oxidation`) are now rejected instead of being silently completed to `PEP[Oxidation]`
- Empty modifications (e.g. `PEP[]TIDE`, `[]-PEPTIDE`) are now rejected at parse time instead of producing an object that raised only later on mass/composition access
- A dangling charge separator (e.g. `PEPTIDE/` or `PEPTIDE/x`) is now rejected with a clear message instead of being silently ignored
- `mass()` and `mz()` now apply global isotope labels (`<13C>`, `<15N>`); previously the fast mass path ignored them and returned the unlabeled mass while `comp()` applied them
- `comp()` no longer drops atom-removing modifications (e.g. Amidated's `O:-1`); a `Counter` accumulation was silently discarding negative element counts, inflating the composition mass
- The `D`/`T` (deuterium/tritium) shorthands inside a chemical formula (e.g. `[Formula:CD3]`, `[Formula:T2]`) now keep their isotope (H-2 / H-3) instead of collapsing to protium; previously the formula parser mapped `D`/`T` to a bare `H` and dropped the isotope number, so a deuterated formula computed the natural-hydrogen mass (off by ~1.006 Da per D). The `<D>`/`<T>` global-isotope forms were already correct; this fixes only the `Formula:` path. Per ProForma 2.1 §11.3.1, `D`/`T` are the sanctioned shorthands for `2H`/`3H`
- `set_charge` now correctly accepts a `Mod`-wrapped charge carrier (an advertised input type): it serializes the wrapped carrier (e.g. `H:z+1^2`) instead of the raw dataclass repr, so the result matches passing the bare `GlobalChargeCarrier` and no longer produces an uncomputable charge like `PEPTIDE/[Mod(value=GlobalChargeCarrier(...))]`
- `set_charge(0)` now clears the charge to `None` (a charge of 0 is a neutral peptidoform per ProForma 2.1 §11.5) instead of storing a literal `0`, so a 0-charge peptide compares equal to an uncharged one; a `bool` (an `int` subclass) is now rejected instead of serializing as `PEPTIDE/True`
- `GlobalChargeCarrier.to_mz_paf` now renders the adduct sign correctly for a negatively charged (deprotonated) carrier — e.g. `charged_proton(-2)` serializes to `M-2H` instead of the malformed `M+-2H`; the +/- direction is the combination of the charged formula's own sign and the sign of the occurrence count
- `condense_ambiguity_to_xnotation` no longer adds a spurious proton (~1.007 Da) to each condensed region's mass
- `estimate_isotopic_distribution` (averagine) now anchors the monoisotopic peak to the requested mass instead of returning the averagine composition's own mass, which drifted by 10-34 Da; only the envelope shape comes from averagine
- `chem_formula` now handles a single formula string correctly (previously a `str` was treated as a batch of characters and raised `ValueError`)
- `isotopic_distribution` now applies the per-charge electron-mass correction for all formulas; `charge_state` was previously ignored for integer formulas, so every charge state returned identical masses
- `get_regex_match_indices` now uses the match end index for non-zero-length (multi-residue) enzyme patterns, fixing incorrect cleavage sites for motifs longer than one residue
- `SequenceRegion.from_string` now parses ambiguous regions (`(?...)`) correctly instead of raising (previously a residue was skipped for every element)
- `parse_modification` now routes cross-link definitions (e.g. `XLMOD:02001#XL1`, `#BRANCH`) to the cross-linker parser via the label rather than mis-classifying them as ambiguous modifications; modifications containing `|` are no longer misrouted to the cross-linker parser
- `isotopic_distribution` no longer aborts convolution early on a single low-abundance pairing (`break` → `continue`), which could drop significant isotope peaks
- Isoelectric point (`pi`) is no longer clamped to the `[4.05, 12.0]` range, so strongly acidic/basic peptides get correct values; the duplicated pI implementation in the functional API now delegates to the annotation property
- `get_cleavage_sites` now treats an empty enzyme pattern (from an empty/`None` `cleave_on`) as non-specific cleavage, consistent with the explicit `"()"` sentinel
- `TagMass` now preserves the `C:` custom-mass CV prefix through a parse -> serialize round trip (previously it was dropped, indistinguishable from a bare mass)
- Parallel processing now honours the documented auto-detection: with no `method` specified it selects threads on free-threaded (no-GIL) Python and processes otherwise, instead of always using processes
- Standardized the internal isotope convolution abundance threshold to a single value (`1e-14`) across all entry points (previously a mix of `10e-15` and `1e-15`)
- `comp()` now scales a repeated modification's composition by its occurrence count, matching `mass()`; applying the same modification 2-3x at one position previously left `comp()` unchanged, producing a composition off by thousands of Da for common exotic-mod peptides
- The functional `fragment()` API now matches `.fragment()`'s smart charge-state default instead of hardcoding `charges=(1,)`; it previously silently dropped higher-charge fragments and returned the wrong sign of charge for negative-mode precursors (e.g. `/-3` returned `charge_state=1` instead of `[-1, -2]`)
- `find_subsequence_indices`/`coverage`/`modification_coverage` now find overlapping subsequence matches (e.g. `'II'` in `'IIIII'`); a plain `re.finditer` call was skipping past each match and missing overlaps
- `count_residues`/`percent_residues` no longer mutate a caller-supplied `ProFormaAnnotation`; they combined a non-copying `get_annotation_input` with an in-place `condense_static_mods` call
- `convert_ip2_sequence` no longer crashes or produces unparseable output for modifications adjacent to other brackets; leading (N-terminal), trailing (C-terminal), and internal-residue bracket runs are now each handled with the correct dash placement instead of one blanket substitution
- Randomly generated ambiguity intervals (`ProFormaAnnotation.random`) can now reach the sequence's final residue; an off-by-one in the interval-generation range meant no randomly generated interval ever included the last residue
- `generate_partitions` no longer silently overlaps windows when `aa_overlap=0` is requested and the sequence doesn't divide evenly by `num_windows`; the "windows don't fit" fallback now honors the requested step for every window boundary, clamping any window (not just the last) whose start would otherwise run past the sequence end. Previously only the last window was clamped, so an interior window could collapse to an empty span and silently report a property value of `0.0`
- Fragment composition adjustment (`adjust_comp`) now merges ion-type and charge-adduct compositions element-by-element instead of via `Counter.__iadd__`, which silently discarded any element whose resulting count was <= 0; a charge adduct or ion type removing more atoms of an element than the base composition has is now correctly rejected with a `ValueError` instead of silently producing an incomplete composition
- `_validate_coverage_lengths` (used by `annotate_ambiguity`) now correctly rejects mismatched coverage-vector lengths; a chained comparison (`a != b != c`) meant forward/reverse coverage vectors that matched each other in length but disagreed with the sequence length were never caught
- `annotate_ambiguity`'s newly created ambiguity interval (for a mass shift with no pre-existing matching interval) now receives the mass shift as a modification; it previously hardcoded `mods=None`, silently dropping the shift
- `set_charge` and the `charge`/`charge_state`/`charge_adducts` accessors no longer collapse repeated identical charge carriers; keying adducts into a `{carrier: count}` dict with a hardcoded count of 1 meant two identical adducts (e.g. `set_charge(["Na:z+1", "Na:z+1"])`, or a `Mods`/`Mod` carrying an occurrence count > 1) were silently deduplicated to one, so `charge_state`, `mz()`, `mass()`, and serialization were all computed for a singly-charged ion. Repeated carriers are now tallied into their true occurrence count and serialize correctly (`PEPTIDE/[Na:z+1,Na:z+1]`)
- `Fragment.charge_adducts` likewise no longer collapses repeated identical adducts, so a fragment carrying two of the same carrier now reports the correct neutral mass and composition instead of subtracting only one carrier
- `set_charge(Mod(carrier, 0))` (a `Mod`-wrapped carrier with an occurrence count of 0) now clears the charge to `None` (a neutral peptidoform) instead of storing an empty list and serializing the malformed `PEPTIDE/[]` with a dangling slash
- `Mods.get_composition` now merges element compositions negative-count safe, so an atom-removing modification (e.g. `Formula:C2H-2`, Amidated's `O:-1`) no longer has its removed atoms silently dropped; the negative-safe merge is now a single shared helper (`add_composition`/`merge_compositions`) used by every additive composition merge instead of three divergent re-implementations
- ProForma header names (`(>name)`) containing an unbalanced parenthesis (e.g. `(>a(b)PEPTIDE`) parse again; the switch to balanced-parenthesis scanning (added so names like `(>my (special) peptide)` keep their inner parens) had inadvertently made a stray `(` raise `ValueError`. Balanced names still work and a name with no closing `)` at all is still an error
- `set_internal_mods_at_index(validate=True)` no longer re-scans the entire annotation on every single-index set, which made residue-by-residue construction O(n^2); the global ambiguous-label invariant is now only re-validated when the newly-set modifications actually introduce a concrete `#label` (unlabelled modifications cannot violate it)
- `Mods.mods` is cached again (it was briefly a plain `@property`), so the stored modification strings are parsed into `Mod` objects once per instance instead of on every mass/composition/charge access

## [3.1.1]
- Renamed `fragment_masses` to `fast_fragment` and updated related references
- Added mzPAF label serialization to Fragment class
- Updated logo URL and features in README
- Added CITATION.cff for software citation information
- Added Contributor Covenant Code of Conduct
- Added CONTRIBUTING.md with development setup, code style, testing, and PR guidelines
- Added `[project.urls]` to pyproject.toml (Homepage, Documentation, Repository, Issues, Changelog)

## [3.1.0]
- Fast fragment ion calculation
- Fixed slice and modification handling bugs

## [3.0.0]
- Major refactor / overhaul
- Proforma 2.1 compatible
- Proforma annotation methods now return proforma objects such that methods can be chained (factory pattern)
- Most functionality is now accessible through annotation objetcs
- Fasta reader
- Split now splits unambiguous segemnts of the annotation (intervals are not split)
- Extensive tests
- uv backend
- Decoy protein generation methods
- Auto multiproccessing/threading via functional api
- Improved copy perforamnce
- Mod objects are frozen and cached
- Proforma components dataclasses
- Removed custom errors
- Removed score / spectra (will make a seperate apckage termed spextacular)
- Added seperate tacular dependancy to handle element/aa/obo lookup/parsing

## [2.5.1]
- added ambiguity support to coverage funcs

## [2.5.0]
- Added features to annotate peptides depending on fragment ions
- moved convertor functions to peptacular.sequence.convertors
- formatting with black

## [2.4.0]
- added modification_coverage function to sequence_funcs.py
- bug fix for 'convert_ip2_sequence'

## [2.3.0]
- bug fixes for isotopes.py
- added support to isotopes.py for:
  - scaling the intensity of isotopic distribution by an intensity factor
  - custom neutron values
  - optionally return mass of isotopic distributions calculated with neutron offsets
  - merge_isotopic_distributions
- added C13_NEUTRON_MASS = 1.003350 & PEPTIDE_AVERAGINE_NEUTRON_MASS = 1.002856 to constants.py


## [2.2.1]
- removed labile mod check from pt.contains_sequence_ambiguity()

## [2.2.0]
### Added:
- condense_mods function to condense_to_mass_mods.py
- [potentially breaking] updated digestion and span functions to return generator objects
- added sequential_digest, digest_from_config and EnzymeConfig to digest.py
- regex strings which don't have the same start/end site will give a warning such as ([KR])
- fixed regex bug with nterm enzymes
- simplified supported enzyme regexes
- added simple fasta_parser since I kept recreating it in other projects, works with several input types
- some more tests
- added condense_to_mass_mods to mass_calc.py which condenses modifications to a single +/- mass value
- Fixed Docs
- Linting

## [2.0.0]
### Added:
- Full ProForma2.0 support
- proforma.py for handling ProForma strings (full support for ProForma2.0)
- Support for all types of internal fragment ions (ax, ay, bx, bx...)
- isotope.py for generating isotope distributions
- apply_static_mod and apply_variable_mods now support n/c term mods
- gno, resid, and xlmod support
- randomizer.py for generating random proforma sequences
- added mods module to handle loading obo files and finding mods
- added fragmenter to fragment.py

### Changed:
- Terminal modifications notation has been changed to use []- and -[] for N- and C-terminal modifications, respectively
- All internal modifications now use [] notation
- Element masses/isotopes are generated using physics.nist.gov db
- Move static/var mod builders to mod_builder.py
- Moved combinatorics funcs to combinatorics.py
- All public functions are accessible from peptacular base (suggest using import peptacular as pt)
- Most functions now support a ProFormaAnnotation object
- Improved digest and fragment performance
- Improved docs

## [1.3.0]
### Added:
- Permutation / Combination / Product functions in sequence.py
- Immonium Ion support to fragment.py

## [1.2.0] 
### Added:
- Added support for custom aa masses to mass.py and fragment.py

## [1.1.1] 
### Added:
- Added isotopes and loss to fragment.py

## [1.1.0] 
### Added:
- speed_test.py to examples and more examples!
- more term functions! pop_c_term & pop_n_term & _get_c_term_index & _get_n_term_index
- split_sequence function to sequence.py, which splits a sequence into a list of modified single residues
- pop_modifications function to sequence.py, since it's useful!

### Changed
- fragment.fragment now returns a list rather than a generator
- identify_cleavage_sites no longer returns cleavage sites at the beginning or end of the sequence (x)PEPTIDE(x)
- removed the x's appended onto sequence when calculating sites
- fixed bug with non-specific digestion
- made term functions more robust, using _get_c_term_index & _get_n_term_index when possible
- split term.py into term.modification and term.residue (it just got too large)

## [1.0.1] 

### Added:
- term.py - a module for handling terminal modifications
- readthedocs 'sphynx-style' documentation
- examples folder

### Changed:
- split up sequence module into sequence and digest
- split mass module into fragment and mass
- updated term modification notation to use square brackets
- renamed most functions... again
- updated documentation to be compatible with sphynx

## [0.3.0]

### Added:
- score.py - a module for scoring peptide matches

### Changes:
- made Fragment dataclass frozen

### Removed:
- removed indexing code from mass.py


## [0.2.0]

### Added:
- Fragmenter.py - an easier to use peptide fragmenter (returns dataclasses)
- Internal fragment ion support 

## Changes:
- fixed bug with C-Term PTMs

## [0.1.0]

### Added:
- spans.py to handle generating peptide spans from cleavage site info
- masses.py to handle generating peptide masses from sequences

### Removed:
- refseq module - was a dumb idea

## [0.0.6]

### Added:
- fragment ion
- peptide mass and mz functions

### Changes
- converted peptide.py to sequence.py (makes more sense since proteins are also sequences of amino acids)

## [0.0.5]

### Added:
- license
- changelog
- documentation

### Changes
- src-based layout
- pyproject.toml

### Removed
- removed fasta module -> independent FastaFrames package
- calculate_peptide_mass(
- moved streamlit components to root

## [0.0.3]

### Added
- added semi, and non-enzymatic digestion
