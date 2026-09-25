
# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [5.0.0] (2026-09-24)

**Breaking:** major API cleanup (removed/renamed names, keyword-only options, `PeptacularError` exceptions, tacular 2.0). See `docs/migration.rst` for an old -> new table.

### Removed
- 82 names from the top-level `peptacular` namespace: 52 that belong to tacular and 30 that were internal or moved (the full list is in `docs/migration.rst`). The tacular lookups, `*Info`/`*Lookup` classes and literal types (`AA_LOOKUP`, `ELEMENT_LOOKUP`, `UNIMOD_LOOKUP`, `PSIMOD_LOOKUP`, `PROTEASE_LOOKUP`, `FRAGMENT_ION_LOOKUP`, `ElementInfo`, `FragmentIonInfo`, `IonTypeProperty`, `parse_composition`, ...) must be imported from tacular. `IonType`, `NeutralDelta` and `Protease` are still re-exported. Also removed: `Any`, `SEQUENCE_TYPE`, `MODIFICATION_*_TYPE`, `GLOBAL_CHARGE_TYPE`, `ModLocation`, `MassPropertyMixin`, `OboEntity`, `OntologyLookup`.
- The FASTA module: `peptacular.fasta`, `parse_fasta`, `parse_fasta_text`, `iter_fasta`, `FastaSequence`, `FASTA_INPUT_TYPE` and `FastaFormatError`. Read FASTA with fastatacular (`read_fasta`, `FastaReader`); its entries can be passed straight to peptacular's sequence functions (see `HasSequence` below).
- `peptacular.regex_utils`, `get_regex_match_indices` and `get_regex_match_range` (now private `peptacular._regex_utils`).
- `CV_TO_NAME_PREFIX`, `CV_TO_ACCESSION_PREFIX`, `CV_TO_MASS_PREFIX` (now private).
- `ReadableProtocol`, `SupportsStr`, `handle_number_and_intern_mod` and `utils.get_mods`.
- The deprecated `FLIXIBILITY_SCALES` alias; use `FLEXIBILITY_SCALES`.
- The `peptacular.isotope.isotopic_distribution` alias; use `brain_isotopic_distribution`. `pt.isotopic_distribution` is unchanged.
- The `enzyme_regex=` keyword of `digest`, `cleavage_sites` and `EnzymeConfig` (used by `sequential_digest`).

### Changed
- Requires `tacular>=2.0,<3` (was `>=1.2,<2`). The re-exported protease enum follows tacular's rename: `pt.Proteases` -> `pt.Protease`.
- `PROTON_MASS`, `ELECTRON_MASS` and `NEUTRON_MASS` are re-exported from `tacular.constants` (CODATA 2018) instead of being defined here, and `C13_NEUTRON_MASS` is `tacular.constants.C13_C12_MASS_DIFF` (unchanged, 1.00335483507). The values move slightly: `PROTON_MASS` 1.00727646688 -> 1.007276466621 (-2.6e-10 Da), `ELECTRON_MASS` 0.00054857990946 -> 0.000548579909065 (-4.0e-13 Da), `NEUTRON_MASS` 1.00866491597 -> 1.00866491595 (-2.0e-11 Da). An m/z at charge z shifts by at most 2.6e-10 Da, far below any instrument's resolution.
- The Da/ppm switch is named `tolerance_unit` everywhere, as across tacular-omics, and typed with `tacular.types.ToleranceUnit` (no local copies). Calls into tacular use its renamed `tolerance_unit=` keyword.
- Every public module declares `__all__`; `peptacular.__all__` is explicit and tested.
- `parallelMethod`/`parallelMethodLiteral` are renamed `ParallelMethod`/`ParallelMethodLiteral`.
- `n_workers`, `chunksize` and `method` are keyword-only on every parallel function.
- `isotopic_distribution` and `estimate_isotopic_distribution` take `sequence=` instead of `annotations=`.
- Digestion takes `enzyme=`, a protease name or a compiled `re.Pattern`. A string is only looked up as a protease name; an unknown name raises the new `UnknownEnzymeError` (a `PeptacularError`) instead of being used as a regex. `digestion.core.resolve_enzyme` does the lookup. `EnzymeConfig` has `enzyme` instead of `enzyme_regex` and is frozen.
- `ProFormaAnnotation.digest`/`simple_digest`/`sequential_digest` are renamed `digest_spans`/`simple_digest_spans`/`sequential_digest_spans`, since they return spans while the functional `pt.digest` returns sequence/span pairs. The batch `"digest"` operation calls `digest_spans`.
- `ProFormaAnnotation` is unhashable (`__hash__ = None`): it is mutable, so a hash could change while it sits in a set or dict. Key on `annot.serialize()`.
- The library raises `PeptacularError` subclasses instead of bare `ValueError`, and `InvalidPositionError` (still an `IndexError` subclass) instead of a bare `IndexError`. The ProForma component parsers raise `ProFormaFormatError`. The MCP layer still raises `ValueError` for pydantic.
- Optional parameters are keyword-only across the public API. Only the input, required arguments and a short list of natural second arguments (`charge` for `mass`/`mz`/`comp`/`isotopic_distribution`, ion type(s) and charge(s) for `frag`/`fragment`/`fast_fragment`, `mods`, `size`, `pH`, ...) stay positional. Every parallel option and every `inplace`/`validate` flag is keyword-only.
- Typo renames in `peptacular.property.data`: `AMIGUOUS_AMINO_ACID_MAP` -> `AMBIGUOUS_AMINO_ACID_MAP`, `surface_accessiblility_janin` -> `surface_accessibility_janin`, `hphob_agros` -> `hphob_argos` (enum member `AGROS` -> `ARGOS`), `hphob_adoberin` -> `hphob_aboderin` (`ADOBERIN` -> `ABODERIN`).
- `get_mod_type` and `_resolve_mod_types` raise `TypeError` (not `ValueError`) for an argument of the wrong type. A string naming no mod type still raises `PeptacularError`.
- `enzyme=""` no longer means a nonspecific digest; it raises `UnknownEnzymeError`. Use `enzyme="unspecific"` or `pt.nonspecific_digest`.
- `calculate_composition=` is renamed `calculate_with_composition=` wherever it appears (`mass`, `mz`, `fragment`, ...).
- `brain_isotopic_distribution(chemical_formula, ..., charge_state=)` is now `brain_isotopic_distribution(formula, *, ..., charge=)`.
- `coverage`, `percent_coverage` and `modification_coverage` methods take `subsequences=` instead of `annotations=`.
- `annot[int]` raises `UnsupportedOperationError` with a hint (`annot[i:i+1]` or `annot.stripped_sequence[i]`) instead of a bare `TypeError`.
- Bare `KeyError`/`ValueError`/`TypeError` from user input at entry points are `PeptacularError` subclasses. `Mods` snapshots its mapping so its hash stays stable. `BatchResult` is compared by value and is not hashable when it holds an annotation, list or dict.
- `Fragment` is immutable (`__slots__`; assignment raises `dataclasses.FrozenInstanceError`). Build a changed copy with `fragment.replace(mass=...)`, which takes the constructor names and rejects unknown keys. Fragments compare and hash by value (the cached composition is not part of the value). Pickle and copy still work.
- MCP (`peptacular-mcp`): tool arguments and fragment rows use the library's names, with no aliases. `fragment_peptides` takes `ion_types` (was `ion_series`) and `isotopes` (was `isotope_offsets`), and `include` offers `"mzpaf"` (was `"label"`). `enumerate_modifications` takes `max_variable_mods` (was `max_variable_modifications`). Fragment rows use `FRAGMENT_RECORD_KEYS` names: `ion_type`, `position`, `charge_state`, `mass`, `neutral_mass`, `mzpaf` (were `ion_series`, `ordinal`, `charge`, `ion_mass_da`, `neutral_mass_da`, `label`); `get_reference(topic="ions")` rows say `ion_type`. `analyze_peptides` and `compare_peptides` measurements and row keys, and `isotope_envelopes` axis values, are `neutral_mass` and `mass` (were `neutral_mass_da` and `ion_mass_da`); `find_modifications` rows use `mass` and `mass_error` (were `mass_da` and `mass_error_da`; `mass_error_ppm` is unchanged). No MCP name carries a `_da` suffix. The response `contract_version` is `"2.0"`. See the MCP table in `docs/migration.rst`.
- `Fragment.losses` is renamed `Fragment.deltas`, to match the constructor argument: the property, `asdict()`, the string form and the MCP fragment rows (`losses` -> `deltas`) all use `deltas`. `Fragment` options after `charge_state` are keyword-only.
- A monoisotopic proton charge carrier weighs CODATA `PROTON_MASS`, as mzPAF (section 4.4.1), pyteomics and OpenMS use; 4.x `fragment()` and `mass()` used H - e, 1.43e-8 Da lighter. Monoisotopic `mass()`, `mz()` and `fragment()` values move by +1.43e-8 Da per charge in mass and +1.43e-8 Da in m/z (-1.43e-8 Da for deprotonated ions), and `fast_fragment` (which already used `PROTON_MASS`) now matches `fragment()` to within 1e-9 Da. A mass summed from a charged composition adds `HYDROGEN_BINDING_MASS` per proton, so the composition and mass paths agree. Average masses still use average H minus an electron and are unchanged.
- Faster: `fragment()` terminal series sum a per-residue mass vector instead of slicing the annotation per ion, and build the isotope/delta/loss products and ion-type lookups once per series instead of once per ion (49 -> 1,298 peptides/s for modified peptides, b/y at charges 1 and 2: indicative, direct run on one core); digest functions return substrings for plain sequences and skip parsing plain strings; `comp` counts residues once and scales each residue composition. Composition mode, formula deltas, isotope swaps and static/isotope/charged mods still take the slicing path; results are equal to within 1e-9 Da.
- `multiprocessing` and `concurrent.futures` are imported lazily, which speeds up `import peptacular`.
- mzPAF neutral-loss labels use the canonical mzPAF names for the known losses (`-NH3`, `-H2O`, `-H3PO4`, `-HPO3`, `-HCONH2`, `-HCOOH`, ...) and write any other formula in Hill order (`+NaS`, `+[13C2]H2`). 4.2.0 wrote `-H3CON` and `-H2CO2`, and tacular 2.0 would have given `-H3N`; every other label equals 4.2.0 on a 400-peptide x 10-loss differential (`tests/reference/mzpaf_neutral_losses_4_2_0.json`).
- A `neutral_deltas=` loss that is impossible for an ion (H3PO4, HPO3 or SO3 from an unmodified S/T/Y) skips that ion instead of aborting the whole `fragment()` call with `InvalidAdjustmentError`. An explicit `deltas=` entry still raises.
- More keyword-only options: every `ProFormaAnnotation` constructor option after `sequence`; `Interval(ambiguous=, mods=, validate=)`; `validate`/`inplace` on `Interval.set_mods`/`append_mod`/`extend_mods`; `Fragment.to_mzpaf(include_sequence=)` and `serialize(format=, include_sequence=)`; the options of `AnnotationProperties.calc_property`/`property_windows`/`property_partitions`; `get_mass(*, monoisotopic=)` on `Mod`, `Mods` and every ProForma component class (as in tacular 2.0), and the options of `ChargedFormula.from_string`/`serialize`/`from_composition`, `FormulaElement.from_string` and `ModificationTags.validate`. `tests/test_v5_signatures.py` now walks every public class.
- `EnzymeConfig.semi_enzymatic` is renamed `semi` (the name `digest()` uses) and its options are keyword-only.
- `Interval.append_mod` returns the interval, or the new copy when `inplace=False` (it returned None).
- `pt.parse` raises a `TypeError` naming the accepted inputs for bytes, None or a number, and accepts an object with a str `sequence` (a FASTA entry).
- `fragment()` and `fast_fragment()` accept a single ion type, charge, isotope or neutral delta without a list. A multi-letter ion type string such as `"by"` is one ion type; 4.x iterated it letter by letter.
- `PeptidoformIon.get_mass` and `CompoundPeptidoformIon.get_mass` raise `UnsupportedOperationError` with a hint instead of a bare `NotImplementedError`; their dead `get_composition` methods are removed.
- `ProFormaAnnotation`'s mass and fragment internals moved to private modules. A subclass that overrides private methods (`_frag*`, `_fragment*`, `_base_*`, `_mass_and_charge`) is no longer called by `mass()`, `frag()` or `fragment()`; override the public methods instead.

### Added
- `fragment_arrays()` and `ProFormaAnnotation.fragment_arrays()`: the ions of `fragment()` as a dict of numpy columns (`FRAGMENT_ARRAY_KEYS`: `peptide_index`, `ion_type`, `position`, `end_position`, `charge_state`, `mz`, `mass`, `isotope`, `isotope_label`, `delta_label`, `delta_mass`), one row per ion, in the same order and with the same float values. `pl.DataFrame`/`pa.table`/`pd.DataFrame` take the dict as is. Plain a/b/c/x/y/z series are computed with numpy prefix sums (about 10x faster than `fragment()` plus a table for 10k peptides); other ions go through `fragment()`. Needs the new `numpy` extra (`pip install "peptacular[numpy]"`); without numpy it raises `MissingOptionalDependencyError` naming that command.
- `ProFormaAnnotation.map_isotopes()`: the global isotope labels as `{element: isotope}` (`<13C>PEP` gives `{C: 13C}`). It replaces the internal `_map_isotopes()`, which is removed (no alias).
- `HYDROGEN_BINDING_MASS` (`PROTON_MASS - (HYDROGEN_MASS - ELECTRON_MASS)`, 1.43e-8 Da), the term a composition-based charged mass adds per proton so it lands on `PROTON_MASS` (see `docs/mass_calculation.rst`).
- `Fragment.replace(**changes)`, `Fragment.__eq__` and `__hash__` (by value). `replace` and the `Fragment` constructor accept the values `deltas`, `isotopes` and `charge_adducts` return.
- `digest_records` and `fragment_records` return digest and fragment results as lists of plain dicts (one per peptide or ion) for `pandas.DataFrame(...)` or `polars.DataFrame(...)`; column names are in `DIGEST_RECORD_KEYS` and `FRAGMENT_RECORD_KEYS`. Every value is a str, int, float, bool or None, one type per column: an internal ion's span is split into `position` and `end_position`. Digest rows carry the input's `accession` (or a PEFF entry's `db_unique_id`); `digest_records` takes one protein, a list or a generator, and `fragment_records` takes one peptide's fragments or the nested list a batch `fragment` returns. peptacular does not depend on pandas or polars. See `docs/records.rst`.
- Localization isomers: `localization_isomers` (also an annotation method) expands unknown-position mods (including `^n`), ranges and `#label` groups into concrete placements, one mod per residue (a `#label` group that lists a residue already carrying another mod raises `PeptacularError`), keeping each group's label and score on the placed mod; `max_isomers=` caps the expansion (default `DEFAULT_MAX_ISOMERS` = 10,000, `None` for no limit). `candidate_sites(peptide, mod, *, residues)` places a mod on residues you name. `site_determining_ions` returns the fragments unique to each isomer against all others, and `pairwise_site_determining_ions` returns, for each ordered pair `(i, j)`, the ions of `i` that `j` cannot explain (the Ascore/PhosphoRS comparison). Both match within `tolerance` in `tolerance_unit` (`"da"` or `"ppm"`, typed `tacular.types.ToleranceUnit`). Candidate sites come only from the ProForma string. See `docs/localization.rst`.
- `PeptacularKeyError` (a `PeptacularError` and a `KeyError`) and its subclass `UnknownElementError`; `UnknownEnzymeError` now subclasses `PeptacularKeyError`.
- `HasSequence`: a protocol for any object with a `.sequence` string. Sequence functions, `batch`, `iter_batch` and `diagnose` accept such objects (fastatacular and PEFF entries) directly, with no dependency on those packages.

### Fixed
- MCP `fragment_peptides`: a formula delta that some ions cannot lose (`"H3PO4"` on y1 of `PEPS[Phospho]TIDE`) no longer fails the whole call. Those ions are skipped and listed in one short `impossible_ions_skipped` diagnostic. `"-H3PO4"` is accepted as the same loss, `"+HPO3"` is a gain, and a formula that does not parse is rejected up front with an actionable message. The tool description states the sign rule. An impossible adjustment elsewhere reports a one-line message instead of the full composition.
- Loss and isotope ions of peptides with named modifications keep the mods' listed database masses. With `deltas=` (formula losses such as `"H2O"`), `neutral_deltas=` or `isotopes=`, `frag()`, `fragment()` and `mass()` summed the mods' elemental compositions instead, so `PEM[Oxidation]TIDEK` y6-H2O was 3.8e-7 Da off the plain y6 minus water. A loss or isotope peak is now exactly the plain ion plus its delta. **Negative-mode masses change:** negative charges took the same composition path even for plain ions, so every negative-charge mass, m/z and fragment moves by up to about 3e-6 Da monoisotopic and about 1e-3 Da average (`PEM[Oxidation]TIDEK/-2` mass goes from 975.4230099731 to 975.4230103537, now exactly the neutral mass minus two protons, as in positive mode). Composition is still used for `calculate_with_composition=True` and under global isotope labels (`<13C>`).
- `Fragment.to_mzpaf` wrote a formula gain (`deltas={"H2O": -1}`) as a loss (`-H2O`); it now writes `+H2O`. A numeric mass delta is written as a signed fixed-point mass rounded to 6 decimals (`b2-34.0` for `{-17.0: 2}`, `+0.00001` for `1e-5`) instead of raising.
- `Fragment.to_mzpaf` and `serialize(format="mzpaf")` write a negative charge signed (`y3{IDE}^-1`, `b3{PEP}^-2`) instead of as a bare magnitude (`^1`), which parsed back as a positive ion 2 Da heavier per charge. This reverses the 4.x change. The new `signed_charge=False` keyword writes the magnitude only (mzPAF 1.0.1 section 4.8), matching paftacular's `serialize(signed_charge=False)`. That form is only valid alongside negative-mode spectrum metadata (at z=-1 it has no charge suffix).
- The mzPAF label of an immonium ion dropped a terminal modification on its residue, so the label's mass was wrong: `[Acetyl]-PEPTIDE` at position 1 gave `IP`. The terminal mod is now written as the immonium modification (`IP[Acetyl]`, `IE[Amidated]`), as paftacular 2.0 does. A global fixed mod that applies to the residue is folded in the same way (`<[Oxidation]@P>` gives `IP[Oxidation]`, not `IP`, which was 15.995 Da off). A global isotope label is written as isotope shifts, one per labelled atom of the final ion, counted after the immonium offset, formula deltas and any atoms a charge carrier removes (`<13C>P` gives `IP+4i13C`, not `IP`, which was 4.013 Da off; `<D>P` gives `IP+7i2H`, and `IP+6i2H^-1` at charge -1; `<15N>K` with an NH3 loss gives `IK-NH3+i15N`). Mass-only deltas, `+i` shifts and adduct atoms are not labelled, so a heavy carrier is counted once, in its adduct (`<2H>P` with `D:z+1` gives `IP+7i2H[M+[2H]]`). An isotope label plus a mass-only mod on the residue raises, since the atoms cannot be counted. More than one modification to write (for example `[Acetyl]-P[Oxidation]`) raises `PeptacularError`. The mod tag is now the plain name, so `P[U:Oxidation]` gives `IP[Oxidation]` (was `IP[U:Oxidation]`, same mass).
- Full-length satellite ions ignored the terminal modification on the residue whose side chain is cleaved: v, w, wa, wb (and their residue-specific variants) dropped the N-terminal mod, and d, da, db dropped the C-terminal mod, so `[Acetyl]-PPA` v3 had the same mass as `PPA` v3. The terminal mod sits on the backbone and now stays, so every full-length ion type carries both terminal mods like a/b/c/x/y/z. Shorter ions are unchanged.
- A hydride charge carrier (`charge="H:z-1"`) counted as a proton: `is_protonated` was True, so the fragment's mzPAF label was `y3{IDE}^-1`, which reads back as the deprotonated ion 2.016 Da lighter. `ChargedFormula.is_protonated` (and so `GlobalChargeCarrier.is_protonated`) is now True only for plain hydrogen carrying +1 per atom. A hydride is kept as an adduct and written `y3{IDE}[M+H]^-1` (`[M+2H]^-2` for two), which paftacular 2.0 parses to the same m/z. The fragment mass was already right. Charge carriers that sum to zero (`["H:z-1", "Na:z+1"]`) now make `to_mzpaf` raise `PeptacularError` instead of writing a label that reads as +1, as paftacular does.
- d, da, db, w, wa and wb ions of a residue that carries a modification (explicit, or from a global fixed mod such as `<[Oxidation]@V>`) were built as if the mod stayed on the side-chain remnant: `PEPV[Oxidation]K` d4 gave `d4{PEPV[Oxidation]}`. These ions keep part of the cleaved side chain, so they are not defined there: `frag()` raises `PeptacularError`, as paftacular 2.0 does, and `fragment()` leaves them out. v ions lose the whole side chain and still drop the mod.
- `fragment()` raised `InvalidAdjustmentError` for a whole series when one ion could not exist, such as the one-residue d1 of `V-[Amidated]` or a1 of `G-[Amidated]` at charge -1. Any series now leaves that ion out; an explicit `frag()` of it, or a `deltas=` the ion cannot take, still raises. An ion the caller's `isotopes=` do not fit (`{"15N": 3}` on b1) is left out too, and `fragment()` raises only when no requested ion is left.
- `Fragment.to_mzpaf` and `serialize(format="mzpaf")` of an uncharged fragment (charge 0, or `frag()` without a charge) wrote a label with no charge suffix, which reads as +1. They raise `PeptacularError`; build the ion with `charge=1`. `fragment_records` gives `mzpaf=None` for it.
- Several charge carriers on a labelled ion gave a composition that depended on their order (`["H-1:z-1", "H:z+1", ...]` on `<D>PEK` could raise). Carriers are summed per element before they are applied.
- Negative charge on a deuterium-labelled peptide (`<D>PEK`, `<2H>PEK`) raised `InvalidAdjustmentError`: deprotonation removed an ¹H the labelled composition did not hold. A charge carrier that removes atoms now takes the isotope the ion holds (a deuteron here), so the composition stays valid and matches the mass. A carrier that names an isotope (`D-1:z-1`) never takes a light atom instead: on `PEK` it still raises.
- The fixed offsets of the ax and bx internal ions are written in Hill order like every other delta: `-H2` and `+CO-H2` instead of `-2H` and `+CO-2H`. This matches paftacular.

## [4.2.0] (2026-09-23)

### Fixed
- `fast_fragment` put static N-terminal mods (`<[Carbamidomethyl]@N-term>`) on the last residue and C-terminal ones on the second-to-last. It now matches `frag()`.
- Immonium and internal fragments counted neutral-loss sites on the whole parent sequence instead of the fragment, so they got losses the fragment cannot have (or raised `InvalidAdjustmentError`).
- `simple_digest`/`generate_regex` dropped a site at the sequence ends when `restrict_before` (C-terminal cleavage) or `restrict_after` (N-terminal cleavage) was set: `KAAAKAAA` with K not after P now cleaves at 1 and 5.
- `sequential_digest` reported only the first enzyme's missed cleavages. The count is now the uncut sites of every enzyme inside the span.
- `is_subsequence(..., order=False)` raised `KeyError` when the subsequence had a residue the sequence lacks; it returns `False`.
- `ProFormaAnnotation.__hash__` depended on modification order while `==` did not, so equal annotations could hash differently. Comparing an annotation with a non-annotation now returns `False` instead of raising `NotImplementedError`.
- Modification names ending in `(...)` lost that suffix: it was read as a localisation score. `K[U:Label:13C(6)]` raised `UnknownModificationError`, `N[HexNAc(2)]` weighed as one HexNAc, and `C[L-cystine (cross-link)#XL1]` did not resolve. Per ProForma 2.0 a score only follows a `#group` label, so `(n)` is now part of the name unless it follows one.
- `isotopes=n` failed with `InvalidAdjustmentError` on ions with no light atoms left to swap (y1 of `K[Formula:[13C6]C-6]`) and was silently ignored under a global `<13C>` label. The offset is now added as a mass delta when no composition is requested; asking for more heavy atoms than the ion has still raises.
- The C-terminal pKa values of Glu and Gln were swapped (E 2.17, Q 2.19). They now match the cited table (E 2.19, Q 2.17), which shifts `charge_at_ph` and `pi` slightly for peptides ending in E or Q.
- `parse_chimeric` read cross-linked peptidoforms (`A//B`) as separate chimeric ions, so `serialize_chimeric` wrote them back as `A+B`, a different meaning. Cross-links are not supported, so it now raises `UnsupportedOperationError`.
- Satellite ions (d, v, w and the residue-specific da/db/wa/wb variants) summed every residue of the fragment and then added the side-chain remnant, counting the cleaved residue twice. They now use the mzPAF 1.0.1 sum of the other residues (n-1 for d, c-1 for v and w) plus the remnant; the cleaved residue's modifications leave with its side chain. Requesting `d` or `w` now also yields the generic ion (it returned only da/db or wa/wb, so `fragment("SAMPLER", ion_types=["d"])` was empty), db now forms on the last residue like da and wb on the first like wa, a residue-specific match no longer stops later positions in `fragment()`, and `frag(position=...)` checks the fragment's own terminal residue. Masses also follow tacular's corrected d/v/w offsets (tacular 5ffebec).
- `Fragment.to_mzpaf` wrote the Biemann z ion as mzPAF `z`, which mzPAF 1.0.1 defines as the z-dot radical, so the label parsed back 1.008 Da heavy; z-dot was written as `z.`, which is not mzPAF. z-dot is now `z`, and Biemann z, z+H and c-H are written as `z-H`, `z+H` and `c-H`, matching paftacular's `to_mzpaf`.
- `C13_NEUTRON_MASS` was rounded to 1.003350; it is now the AME2020 13C-12C difference, 1.00335483507.
- `parse` and `parse_chimeric` raised a bare `ValueError` for invalid ProForma, and `parse` raised `ValueError` for chimeric or cross-linked input. They now raise `ProFormaFormatError` (a `ValueError` subclass, so existing `except ValueError` still works) and `UnsupportedOperationError` respectively.
- `parse` looped forever on a `?` that does not follow a modification (`?[Phospho]PEPTIDE`). It now raises `ProFormaFormatError`.
- Malformed modification, glycan, isotope, static-mod and adduct strings parse lazily, so they surfaced from `mass()` as a bare `ValueError` (or `KeyError` for an unknown isotope such as `<113C>`). They now raise `ProFormaFormatError`. A non-numeric adduct multiplier (`/[Na:z+1^x]`) was read as `^1`; it is now a parse error.
- `pt.shift("PEPTIDE", "a")` raised `TypeError: not all arguments converted during string formatting`. A non-integer `n` now raises `TypeError: n must be an int`.
- `fragment(ion_types=["q"])`, `frag(ion_type="q")` and `mass(..., ion_type="q")` leaked the enum's `'q' is not a valid IonType`. They raise `UnsupportedOperationError` listing the valid ion types.
- `pt.digest(seq, "notanenzyme")` silently used an unknown protease name as a regex and returned the sequence uncut. A string that is neither a known protease nor contains regex metacharacters now emits a `UserWarning`; behaviour is otherwise unchanged.
- `append_mods` wrote a list or tuple value as one bracketed string (`pt.append_mods("PEMTIDE", {1: ["Oxidation"]})` gave `PE[['Oxidation']]MTIDE`, invalid ProForma), and `extend_mods` iterated a bare string character by character (`{1: "Oxidation"}` gave `PE[O][x][i]...`). `append_mods` now adds each item of a list/tuple (a `(mod, count)` pair is still one mod with a count), and `extend_mods` treats a string, number or `Mod` as one modification. The same applies to the `ProFormaAnnotation.append_*`/`extend_*` methods.
- Bare `ValueError`s on main public paths are now typed (all still `ValueError` subclasses): `mass("")` raises `CompositionError`; out-of-range slices and `frag(position=...)` raise `InvalidPositionError`; invalid `parse_fasta_text`/`iter_fasta` input raises `FastaFormatError`.

### Added
- `ProFormaFormatError`, raised for strings that are not valid ProForma.
- `PeptacularError(ValueError)`, the common base of `ProFormaFormatError`, `UnknownModificationError`, `CompositionError`, `InvalidAdjustmentError`, `UnsupportedOperationError` and the two new classes below. `except ValueError` still catches all of them.
- `InvalidPositionError` (slice index or fragment position outside the sequence) and `FastaFormatError` (invalid FASTA text), both `PeptacularError` subclasses.
- Docstrings for the property functions (`pi`, `charge_at_ph`, `hydrophobicity`, ...), `set_mods`/`append_mods`/`extend_mods`/`remove_mods`, `isotopic_distribution`, `left_semi_digest`/`right_semi_digest` and `ModType`; `set_start_method` has a `-> None` return annotation. The error classes are documented in the streaming guide, README, `llms.txt` and `llms-full.txt`.

### Changed
- Requires `tacular>=1.2,<2` (was `>=1.1.0`): tacular 1.2 bundles PSI-MOD 1.039.0 and the corrected d/v/w satellite-ion offsets.
- An unknown isotope label (`<113C>PEPTIDE`) now raises `ProFormaFormatError` instead of `KeyError`. Code that caught `KeyError` for this case must catch `ProFormaFormatError` (or `ValueError`).
- `Diagnostic.code` from `pt.diagnose`/`pt.batch(errors="collect")` follows the new types: an empty sequence reports `unavailable_composition` and an unknown ion type `unsupported_operation` (both were `calculation_error`).
- Reference-value tests (`tests/reference/`) against pyteomics, Biopython, ExPASy ProtScale and the ProForma 2.0 spec examples, and Hypothesis property tests (`hypothesis` added to the dev group).

## [4.1.0] (2026-09-23)

### Fixed
- `pt.get_mods` is the documented functional `get_mods(sequence, mods)` again. The internal `peptacular.utils.get_mods` (a `ModType` list helper) had shadowed it through the star imports.
- `chem`, `constants`, `regex_utils`, `spans` and `utils` define `__all__`, so `from peptacular import *` and `pt.` no longer expose stdlib/typing helpers (`re`, `sys`, `Counter`, `Sequence`, `Final`, ...).
- `POLARITY_SCALES` was defined twice in `property/data.py`; the unused string-keyed copy is removed. The enum-keyed table (the one already in effect) is unchanged.
- The MCP server instructions and `peptacular://conventions` resource name `spxtacular` (not "Spectacular") as the package for observed spectra.
- Added docstrings to `ProFormaAnnotation`, `Interval`, `Fragment` and `get_mod_type`.
- `just lint` now checks `tests/` as CI does, and `just check` runs the new `just format-check`.

### Deprecated
- `FLIXIBILITY_SCALES` is renamed `FLEXIBILITY_SCALES`. The old name still works and emits a `DeprecationWarning`.

## [4.0.1] (2026-09-23)

### Changed
- Capped the tacular requirement below the next major version (`tacular>=1.1.0,<2`) so a breaking tacular release cannot reach installs before it is tested.

### Maintenance
* Publish from GitHub Actions with PyPI trusted publishing (`publish.yml`);
  release metadata is checked against the tag.
* Keep `__version__`, `CITATION.cff` and `.zenodo.json` in sync with
  `scripts/release_version.py` (`just set-version X.Y.Z`).
* CI tests Python 3.12-3.14 on Linux plus macOS and Windows, the lowest
  direct dependency versions, and the built wheel.

## [4.0.0] (2026-09-13)

### Added
- Added one-command version synchronization for package and citation metadata, with CI checks for release tags, changelog dates, and built distribution versions.

### Changed
- Accelerated ordinary precursor and neutral mass/m/z calculations by avoiding fragment allocation and charge-override copies, and using direct residue mass lookups.
- Replaced exact-mass isotope fine-structure convolution with the BRAIN Newton-Girard recurrence. Isotope distributions now contain one aggregated peak per populated nominal neutron offset with an exact probability-weighted center mass.
- Isotope envelopes now use adaptive sizing when `max_isotopes` is omitted. Weak leading peaks are retained through the last peak that meets the relative-abundance threshold.
- Peptide averagine now subtracts the fixed terminal composition before applying mass-scaled elemental ratios.

### Fixed
- Validate isotope counts before cache lookup so cached integer counts cannot cause boolean or floating-point counts to be accepted.
- Report a clear validation error when calculating mass, m/z, composition, or fragments from an empty sequence.

### Removed
- Removed the `distribution_resolution`, `use_neutron_count`, and `conv_min_abundance_threshold` isotope arguments. Every `IsotopicData` result now provides both center mass and nominal neutron offset directly.
- Removed `IsotopeLookup`. BRAIN calculations use bounded internal caches and no longer require coarse 50 Da mass bins.

### Migration from 3.x
- Remove `distribution_resolution`, `use_neutron_count`, and `conv_min_abundance_threshold` from isotope calculation calls.
- Use each result's `mass` for its probability-weighted center mass and `neutron_count` for its nominal isotope position. Results are aggregated nominal peaks, so callers expecting exact-mass fine structure must use a different calculation method.
- Replace `IsotopeLookup` with direct isotope calculation calls. Existing mass and m/z APIs remain available.

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
