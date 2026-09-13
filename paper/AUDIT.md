# Manuscript audit, 13 September 2026

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

The source checked was HEAD `c81f282` plus the current working-tree changes on
`codex/release/4.0.0`. Example execution used Python 3.12.3, Tacular 1.1.0,
pandas 3.0.0, and the local Peptacular source.

| Check | Result |
| --- | --- |
| Three Python code blocks | Executed successfully. Checked the displayed rounded mass and m/z values, serialized string, and DataFrame assignment. |
| Table examples | 42 examples total. All 39 marked supported or preserved passed `validate=True` and serialization round trips. The chimeric example was checked through `parse_chimeric()`. |
| Table mass behavior | Calculated masses for the supported/preserved examples. The `BZJX` example raised the expected ambiguity error. Three linked-chain examples were excluded in accordance with the stated limitation. |
| Citation integrity | 20 entries, no missing citation keys and no unused entries. |
| Manuscript length | Approximately 1,058 prose words, excluding headings, table, code blocks, bibliography, and YAML. Approximately 1,444 words including headings, table, and code, still excluding bibliography and YAML. Counts were extracted from Pandoc's document structure. |
| JOSS build | PDF and JATS generated successfully with the installed `openjournals/inara:latest` image, with networking disabled for the render. |
| Visual proof | All six pages inspected. Corrected overlapping table headers and literal escape characters in ProForma examples. |
| Author metadata | Both names and ORCIDs checked in JATS. John R. Yates III has separate given-name, surname, and suffix fields. |

The image used for the build was
`sha256:2415076f0ef85d98dca68707ec3b1263f913487f55aa59e3448704ab7087789d`.
The regenerated files are `paper/paper.pdf` and `paper/jats/paper.jats`.
The PDF remains a local JOSS draft, including journal-generated placeholder
publication metadata. Its front-page author names and JATS metadata are correct.
The journal-generated footer citation should also be checked in the official proof.

These checks establish the behavior of the displayed examples. They are not a
complete ProForma conformance certification, a new full-suite test run, or an
independent validation of all scientific algorithms.

## Items to resolve before the final submission update

1. **Release and review version.** The [public review](https://github.com/openjournals/joss-reviews/issues/11277)
   identified v3.0.0 at audit time. Release preparation now aligns the package,
   citation metadata, and dated changelog to 4.0.0. After publication, update the
   JOSS review version and its software archive through the editorial workflow.
   The archive was still pending in the review record when checked.
2. **Human review, resolved.** The author confirmed that human authors reviewed,
   edited, and validated all AI-assisted work and made the core design decisions.
   The manuscript now states this explicitly. Opus 5, Fable 5.1, Sol, and Astra
   are recorded as supplied by the author.
3. **Tacular archive citation.** The retained citation points to historical
   v1.0.1, while the tested dependency is 1.1.0. A historical software citation is
   distinct from the tested environment, but the preferred archive should be
   confirmed. The Zenodo record could not be independently fetched during this
   audit, so its existing metadata was retained rather than guessed.
4. **Author-controlled statements, resolved.** Funding is disclosed through the
   NIH grant list in Acknowledgements. The author confirmed that the funders
   provided financial support only and that the authors have no competing
   interests. Both statements are now included in the manuscript.

The [pre-review discussion](https://github.com/openjournals/joss-reviews/issues/10355)
requested a shorter manuscript, removal of the Mathematics section, and an
explicit References heading. Those structural changes are preserved. The revised
length is consistent with the editor's approximate 1,500-word target and the
[current paper guidance](https://joss.readthedocs.io/en/latest/paper.html).
