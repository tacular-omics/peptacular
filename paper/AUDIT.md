# Manuscript audit, 13 September 2026 (verification refreshed 24 September 2026)

## Update for the 5.0 release, 24 September 2026

The paper was changed only where 5.0 made it wrong or stale. No new features were
added to the text.

- The FASTA sentence was replaced. 5.0 removed `peptacular.fasta`; FASTA is read by
  fastatacular, and any object with a `sequence` attribute (the `HasSequence`
  protocol, which covers fastatacular and pefftacular entries) passes directly to
  sequence functions. Checked by digesting a fastatacular `SequenceEntry` and a
  pefftacular entry from `tests/fixtures/minimal.peff`.
- The Tacular citation now uses the Zenodo concept DOI 10.5281/zenodo.18475556
  (from tacular's CITATION.cff; resolves, record title "tacular: Proteomics
  ontology and reference-data lookups in Python") instead of the v1.0.1 version DOI.
- The composition-path clause now says that composition is used for global isotope
  labels or on request (`calculate_with_composition=True`), that neutral losses are
  added to the listed modification masses, and that protons use the CODATA mass
  (CHANGELOG 5.0, `src/peptacular/constants.py`).
- The date is 24 September 2026. `paper.pdf` and `jats/paper.jats` were rebuilt.
  The orphaned `jats/fig1.png` (the paper has no figure) was removed.

The sections below "Verification performed" record the 13 September revision and
describe the 4.0.0 branch at that time. Their "streaming FASTA" and "20 entries"
statements are superseded by this update.

The manuscript has been revised against the current source, changes since the
v3.0.0 submission, the finalized ProForma 2.1 specification, current competing
software documentation, and the public JOSS review record. The resulting draft
describes this release branch, including the unreleased BRAIN and scalar mass
changes. It should accompany the corresponding release when those changes ship.

## Substantive corrections

| Area | Finding and revision | Evidence |
| --- | --- | --- |
| Title and scope | Removed the unqualified compliance claim. Distinguished notation representation from the calculations available for each annotation. Aligned the paper, CITATION.cff, and Zenodo metadata titles. | `src/peptacular/annotation/annotation.py`, `PROFORMA_COMPLIANCE.md`, `COMPLICANCE_NOTES.md` |
| ProForma standard | Replaced old level names with levels 1, 2, and 3. Added labile-location and mixed-glycan examples. Placement constraints are preserved without being enforced. | [Final specification and implementation listing](https://www.psidev.info/proforma) |
| Competing software | Removed outdated characterizations of Pyteomics and RustyMS. Pyteomics handles ProForma mass, composition, and fragmentation. mzcore and its related libraries provide extensive chemistry and fragmentation capabilities. | [Pyteomics sequence documentation](https://pyteomics.readthedocs.io/en/latest/sequences.html), [mzcore](https://github.com/rusteomics/mzcore) |
| OpenMS comparison | Acknowledged current OpenMS ProForma parsing and conversion. Avoided claiming that every C++ interface is exposed in pyOpenMS. | [OpenMS ProForma class](https://openms.de/current_doxygen/html/classOpenMS_1_1ProForma.html), [pyOpenMS peptide guide](https://pyopenms.readthedocs.io/en/latest/user_guide/peptides_proteins.html) |
| Object model | Replaced the factory-pattern description with the actual editable annotation behavior and `inplace=False` copy option. | `src/peptacular/annotation/annotation.py` |
| Batch execution | Explained the 1,000-input automatic threshold, process/thread selection, temporary functional pools, and executor reuse within `iter_batch()`. Removed the claim that functional pools are cached across calls. | `src/peptacular/sequence/parallel.py`, `src/peptacular/batch.py` |
| Mass calculations | Described direct residue lookups, the scalar precursor path, and consistent isotope, charge, and electron-mass handling. No hardware-independent speedup claim was added. | `src/peptacular/annotation/annotation.py`, `cached_comps.py`, `utils.py` |
| Isotope envelopes | Replaced the earlier convolution description with BRAIN aggregated nominal peaks and probability-weighted center masses. Explicitly stated that fine structure is not resolved. | `src/peptacular/isotope.py`, [BRAIN 2.0](https://doi.org/10.1007/s13361-013-0796-5) |
| Recent interfaces | Added concise coverage of JSON validation, streaming FASTA, batch diagnostics, optional adapters, and the local MCP interface. | `CHANGELOG.md`, `src/peptacular/batch.py`, `docs/json_serialization.rst`, `docs/mcp.rst` |
| Dependencies and validation | Matched Python requirements, core dependencies, optional extras, and CI configurations to the repository. Removed unsupported numeric coverage and compile-time safety claims. | `pyproject.toml`, `.github/workflows/python-package.yml` |
| External research use | Added an independently published use of Peptacular for theoretical ion matching in plasma proteomics. The evidence is specifically the Figure 4 caption. | [Malsagova et al., 2026](https://www.nature.com/articles/s41598-026-44729-5) |
| Downstream integration | Added Spxtacular's use of Peptacular fragments in matching and scoring. Identified it as part of the authors' ecosystem, so it is not presented as independent adoption. | [Spxtacular](https://github.com/tacular-omics/spxtacular), its dependency declaration and matching/scoring code, author confirmation |
| AI disclosure | Recorded the author-supplied models: Opus 5 and Fable 5.1 through Claude Code, and Sol and Astra through Codex. Described the assistance and verification performed. Added the author's confirmation of completed human review, editing, validation, and core design decisions. | Author reply, [JOSS AI policy](https://joss.readthedocs.io/en/latest/submitting.html#ai-usage-policy) |

The description deliberately remains a short software paper. Detailed numerical
benchmarks, migration instructions for removed isotope arguments, individual bug
fixes, and the complete API inventory belong in the changelog and documentation.

## Examples and bibliography

- Corrected `<C13>ARE` to `<13C>ARE`. The neutral mass is 388.2384 Da and the
  charge-2 m/z is 195.1265 at the displayed precision.
- Corrected `PEM[Oxidation]TIDE` to 849.343 Da and charge-2 m/z 425.679 at three
  decimal places. The old values were truncated incorrectly for those comments.
- Replaced the unresolvable RESID and XLMOD examples with available entries.
  Replaced the globally labelled `CARBON` example, which contains mass-ambiguous B.
- Used batch mass assignment in the pandas example and retained an explicitly
  modified serine example. Synced the duplicated README calculations.
- Updated the companion compliance checklist to match the paper's terminology,
  examples, and placement-control limitation.
- Added the ProForma 2.0 publication and final 2.1 specification, BRAIN reference,
  current competing-tool documentation, external research paper, and Spxtacular.
- Corrected the PSI-MOD link, which previously pointed to the PSI-MS vocabulary.
- Corrected pyOpenMS to its 2014 journal issue year and the textbook chapter's
  publication type and book title. Completed the Angel et al. page range using
  the [publisher record](https://pubs.rsc.org/en/content/articlelanding/2012/cs/c2cs15331a).
- Removed unused bibliography entries. All 20 remaining entries are cited and
  every manuscript citation resolves in the generated proof.

## Verification performed

Refreshed 24 September 2026 against peptacular `origin/main` `fe76a8a` (the
unreleased 5.0) in a clean worktree, with the paper edits above applied. Examples
ran with Python 3.13.7, pandas 3.0.6 and Tacular from its local `origin/main`
checkout `4a10a5f` (the upcoming 2.0; its version string still reads 1.2.0 because
the bump has not been applied).

| Check | Result |
| --- | --- |
| Three Python code blocks | All executed. `PEM[Oxidation]TIDE` gives 849.343 Da and charge-2 m/z 425.679; the chained edit serializes as `(>Peptacular)PEM[Oxidation]TIDE/2`; the batch gives [928.4026, 388.2384, 451.2454] and m/z [465.2086, 195.1265, 225.6227]; the DataFrame assignment works. The CODATA proton change does not show at the printed precision. |
| Table examples | 42 examples in 40 rows. All 39 marked Y or P passed `validate=True` and exact serialization round trips; `NEEYN+SEQUEN` was checked through `parse_chimeric()`. |
| Table mass behavior | 38 of the 39 Y/P examples return a mass. `BZJX` raises `PeptacularError`, as expected. The two `//` linked examples raise `UnsupportedOperationError`. The intrachain cross-link `EVTK[X:DSS#XL1]LEK[#XL1]SEFD` returns a mass (1461.7239 Da, the DSS mass counted once, equal to `EVTK[X:DSS]LEKSEFD`) and fragments without error, ignoring the link. |
| Citation integrity | 26 entries, all cited, no missing keys. The build reported no citation warnings. |
| Manuscript length | About 1,113 prose words, excluding headings, table, code blocks, bibliography and YAML; about 1,505 words including headings, table and code. Counted from Pandoc's plain-text output. |
| JOSS build | `openjournals/inara` (`sha256:a0414b8b72fd8923917ede614d340d98dc7aa3102aabc2e955e9c02a6100fd62`), `-o pdf,jats`, networking disabled: exit 0, no warnings, 7 pages (as before). |
| Visual proof | Pages 1 and 2 inspected; the changed sentences and the Tacular reference render correctly. |
| Author metadata | Both names and ORCIDs present in JATS; John R. Yates III has separate given-name, surname and suffix fields. |

The regenerated files are `paper/paper.pdf` and `paper/jats/paper.jats`. The PDF
remains a local JOSS draft with journal-generated placeholder publication metadata.

These checks establish the behavior of the displayed examples. They are not a
complete ProForma conformance certification, a full-suite test run, or an
independent validation of all scientific algorithms.

## Items to resolve before the final submission update

1. **Release and review version.** The [public review](https://github.com/openjournals/joss-reviews/issues/11277)
   identified v3.0.0 at audit time. 4.0.0 and 4.2.0 have since been released and
   `origin/main` is heading to 5.0.0. After the 5.0 release, update the JOSS review
   version and its software archive through the editorial workflow.
2. **Human review, resolved.** The author confirmed that human authors reviewed,
   edited, and validated all AI-assisted work and made the core design decisions.
   The manuscript now states this explicitly. Opus 5, Fable 5.1, Sol, and Astra
   are recorded as supplied by the author.
3. **Tacular archive citation, resolved.** The paper now cites Tacular's Zenodo
   concept DOI (10.5281/zenodo.18475556), which always resolves to the latest
   archived version.
4. **Author-controlled statements, resolved.** Funding is disclosed through the
   NIH grant list in Acknowledgements. The author confirmed that the funders
   provided financial support only and that the authors have no competing
   interests. Both statements are now included in the manuscript.

The [pre-review discussion](https://github.com/openjournals/joss-reviews/issues/10355)
requested a shorter manuscript, removal of the Mathematics section, and an
explicit References heading. Those structural changes are preserved. The revised
length is consistent with the editor's approximate 1,500-word target and the
[current paper guidance](https://joss.readthedocs.io/en/latest/paper.html).
