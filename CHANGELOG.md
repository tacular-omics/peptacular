
# Changelog

All notable changes to this project will be documented in this file.

### TODO (Next Release?)
- Take valid mod values from the respective dbs for randomizer
- W/V/D iosn should pop the terminal mods if present? and/or internal mods on first/last aa?
- ensure str values are properly handles with intern and that mod values are cached

## [3.1.2]
### Changed
- Bumped the `tacular` dependency floor to `>=1.1.0`, which fixes several upstream data-consistency bugs: isotope-labelled modification compositions (e.g. `UNIMOD:536`, `Label:13C(2)15N(1)`) that dropped their isotope atoms while keeping the correct mass, all 9 internal-fragment-ion mass offsets (previously shifted so the default "by" internal fragment was `-CO` instead of `0`), and two neutral-loss formulas (Formic acid, Formamide) that were parsed with a dropped repeated-element count.

### Added
- Clearer, more actionable parse/validation error messages aimed at both humans and AI agents: they now name the offending value, state what was expected, and (for parse errors) point at the exact position. Unknown-modification errors include a hint listing valid ways to specify a modification (name, CV accession, formula, glycan, or delta mass).

### Fixed
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