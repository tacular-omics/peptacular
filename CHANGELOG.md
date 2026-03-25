
# Changelog

All notable changes to this project will be documented in this file.

### TODO (Next Release?)
- Take valid mod values from the respective dbs for randomizer
- W/V/D iosn should pop the terminal mods if present? and/or internal mods on first/last aa?
- ensure str values are properly handles with intern and that mod values are cached

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