---
title: 'Peptacular: A Python package for amino acid sequence analysis with ProForma 2.1'
tags:
  - Python
  - Proteomics
  - Mass Spectrometry
  - ProForma
  - Bioinformatics 
authors:
  - name: Patrick T. Garrett 
    orcid: 0000-0002-8434-9693 
    affiliation: 1
  - given-names: John R.
    surname: Yates
    suffix: III
    orcid: 0000-0001-5267-1672 
    corresponding: true
    affiliation: 1
affiliations:
  - name: The Scripps Research Institute, United States
    index: 1
date: 07 February 2026
bibliography: paper.bib
---
# Summary

Mass spectrometry identifies and characterizes proteins by measuring molecules and their fragments. Interpreting these measurements requires software that accounts for chemical modifications, charge states, and isotope composition [@angel-2012]. **Peptacular** is a Python library for representing and analyzing peptide and protein sequences using ProForma notation. It supports sequence editing, mass and mass-to-charge ratio (m/z) calculations, elemental compositions, isotope envelopes, enzymatic digestion, theoretical fragmentation, and physicochemical properties. Its calculation APIs operate on individual peptide chains, while separate interfaces represent chimeric assignments and structured notation. Parsing support and calculation limits are distinguished in Table 1.

# Statement of Need

ProForma standardizes how peptide and protein sequences describe modifications and ambiguity [@leduc-2022]. Version 2.1 extends this notation with additional chemical and structural annotations [@proforma-2026]. Proteomics developers still need to carry those annotations through sequence editing, digestion, and chemical calculations without silently discarding information. Peptacular targets researchers building analysis scripts, theoretical peptide libraries, and proteomics software that require consistent behavior across these operations.

The library provides scalar and batch interfaces for serialized sequences and parsed annotations. Batch results can be assigned to tabular data, including pandas DataFrames [@team-2025]. Streaming input and optional collection of calculation errors support workflows in which some annotations are unresolved or unsuitable for a requested operation.

# State of the Field

Existing packages address overlapping needs. **Pyteomics** [@goloborodko-2013] supports ProForma parsing, mass and composition calculations, and fragment generation alongside its historical modX format [@pyteomics-docs]. **Biopython** [@cock-2009] provides general sequence analysis, including protein properties. The **mzcore** and related RustyMS libraries provide ProForma support, chemical representations, and theoretical fragmentation, with Python bindings for selected components [@Schulte_mzcore]. **pyOpenMS** [@rost-2014] exposes a broader mass spectrometry toolkit. Current OpenMS documentation also describes ProForma parsing and conversion [@openms-proforma].

Peptacular's contribution is a common, editable ProForma annotation model shared by sequence transformations and chemical calculations in Python. This design allows researchers to inspect and modify annotations directly while using consistent charge, modification, and ambiguity handling across operations. A dedicated package keeps this interface focused on sequence analysis. Optional adapters connect it to Pyteomics, psm_utils, and AlphaBase for their supported representations, with checks for conversion losses. These integrations complement the specialized capabilities of the surrounding ecosystem.

# Software Design

Peptacular offers functional and object-oriented APIs. `ProFormaAnnotation` objects support inspection, serialization, and chained edits. Many editing methods modify the object by default and accept `inplace=False` to return a copy. Functional operations accept strings or annotations and support sequence batches.

Execution can be sequential, threaded, or process-based. Automatic functional calls use sequential processing below 1,000 inputs, avoiding worker startup costs for small batches. Larger batches use processes with the GIL enabled and threads with it disabled. Explicit settings override this selection. Functional calls create temporary pools, while `iter_batch()` reuses an executor across chunks within a call and returns ordered results with optional diagnostics.

Lazy modification parsing and bounded caches reduce repeated work. Direct residue mass lookups and a scalar path for ordinary precursor calculations avoid unnecessary fragment objects and annotation copies. Composition-based calculations handle isotope labels and elemental adjustments. Both calculation paths account for intrinsic modification charge and electron mass. BRAIN recurrences calculate aggregated isotope envelopes with a probability-weighted center mass per nominal isotope peak [@dittwald-2014]. Fine structure is not resolved. A separate averagine API estimates envelopes from mass alone.

Shared reference data are supplied by **Tacular** [@garrett-2026-tacular], including Unimod [@creasy-2004], PSI-MOD [@hupo-psi-mod], RESID [@resid], XLMOD [@hupo-psi-xlmod], and GNOme [@gnome]. Embedded data avoid runtime ontology downloads. Calculations require a resolvable mass or composition as appropriate. Unresolved annotations can still be represented, while diagnostics distinguish parsing, validation, and calculation failures.

Versioned JSON serialization preserves annotation structure and validates input against a closed set of supported types. Streaming FASTA input supports plain and gzip files. An optional local Model Context Protocol interface exposes the same sequence operations to agent clients. The core package requires Python 3.12 or later and Tacular, with additional dependencies installed only for optional integrations. Type annotations support static analysis. Continuous integration runs tests, linting, type checks, and package builds across Python 3.12-3.14 and Linux, macOS, and Windows configurations.

# Research impact statement

Peptacular has been used outside its development group. Malsagova and colleagues used it to match theoretical b- and y-ions, including water and ammonia losses, when presenting peptide identifications from plasma proteomics in a study of exercise intensity [@malsagova-2026, Figure 4]. HUPO-PSI also lists Peptacular among Python implementations of ProForma [@proforma-2026]. Within the authors' software ecosystem, Spxtacular uses Peptacular's theoretical fragments in spectrum matching and scoring [@spxtacular-2026]. Peptacular has additionally been used by its developers to prepare figures for a proteomics textbook chapter [@garrett-2025].

# Example Usage

## Object-oriented API

```python
import peptacular as pt

# Parse a sequence into a ProFormaAnnotation
peptide: pt.ProFormaAnnotation = pt.parse("PEM[Oxidation]TIDE")

# Calculate mass and m/z
mass: float = peptide.mass() # 849.343
mz: float = peptide.mz(charge=2) # 425.679

# Chained edits modify the annotation
print(peptide.set_charge(2).set_peptide_name("Peptacular").serialize())
# (>Peptacular)PEM[Oxidation]TIDE/2
```

## Functional API

```python
import peptacular as pt

peptides = ['[Acetyl]-PEPTIDES', '<13C>ARE', 'SICK/2']

# Calculate mass and m/z for all peptides
masses: list[float] = pt.mass(peptides) # [928.4026, 388.2384, 451.2454]
mzs: list[float] = pt.mz(peptides, charge=2) # [465.2086, 195.1265, 225.6227]
```

## Tabular workflow

```python
import peptacular as pt
import pandas as pd

df = pd.DataFrame(
    {
        "seq": ["PEM[Oxidation]TIDE", "ACDEFGHIK", "AS[Phospho]TPEK"],
    }
)

df["mass"] = pt.mass(df["seq"].tolist())
```

# Notation support and calculation limits

**Table 1: Representative ProForma 2.1 notation support**

| S | Feature                    | Example                                      | § [Support] |
| - | -------------------------- | -------------------------------------------- | ----------- |
| Y | Amino acids (+UO)          | `AAHCFKUOT`                                  | 6.1 [1]     |
| Y | Unimod names               | `PEM[Oxidation]AT`                           | 6.2.1 [1]   |
| Y | PSI-MOD names              | `PEM[monohydroxylated residue]AT`            | 6.2.1 [1]   |
| Y | Unimod numbers             | `PEM[UNIMOD:35]AT`                           | 6.2.2 [1]   |
| Y | PSI-MOD numbers            | `PEM[MOD:00425]AT`                           | 6.2.2 [1]   |
| Y | Delta masses               | `PEM[+15.995]AT`                             | 6.2.3 [1]   |
| Y | N-terminal modifications   | `[Carbamyl]-QPEPTIDE`                        | 6.3 [1]    |
| Y | C-terminal modifications   | `PEPTIDEG-[Methyl]`                          | 6.3 [1]     |
| Y | Labile modifications       | `{Glycan:Hex}EM[U:Oxidation]EV`              | 6.4 [1]     |
| Y | Multiple modifications     | `MPGNW[Oxidation][Carboxymethyl]PESQE`       | 6.5 [1]     |
| Y | Information tag            | `ELV[INFO:AnyString]IS`                      | 6.6 [1]     |
| Y | Ambiguous amino acids      | `BZJX`                                       | 7.1  [2]    |
| Y | Prefixed delta masses      | `PEM[U:+15.995]AT`                           | 7.2  [2]    |
| Y | Mass gap                   | `PEX[+147.035]AT`                            | 7.3  [2]    |
| Y | Formulas                   | `PEM[Formula:O]AT`, `PEM[Formula:[17O1]]AT`  | 7.4  [2]    |
| Y | Mass with interpretation   | `PEM[+15.995|Oxidation]AT`                   | 7.5  [2]    |
| Y | Unknown mod position       | `[Oxidation]?PEMAT`                          | 7.6.1  [2]  |
| Y | Set of positions           | `PEP[Oxidation#1]M[#1]AT`                    | 7.6.2  [2]  |
| Y | Range of positions         | `PRT(ESFRMS)[+19.0523]ISK`                   | 7.6.3  [2]  |
| Y | Position scores            | `PEP[Oxidation#1(0.95)]M[#1(0.05)]AT`        | 7.6.4  [2]  |
| Y | Range position scores      | `(PEP)[Oxidation#1(0.95)]M[#1(0.05)]AT`      | 7.6.5  [2]  |
| Y | Amino acid ambiguity       | `(?VCH)AT`                                   | 7.7  [2]    |
| Y | Modification prefixes      | `PEPM[U:Oxidation]AS[M:O-phospho-L-serine]`  | 7.8  [2]    |
| Y | Labile locations           | `{Phospho#g1}EMEVS[#g1]`                    | 7.9 [2]    |
| Y | RESID modifications        | `EM[R:L-methionine sulfone]EM[RESID:AA0251]` | 8.1 [T]     |
| Y | Names                      | `(>Heavy chain)EVQLVESG`                     | 8.2 [T]     |
| Y | XL-MOD modifications       | `EVTK[X:DSS]LEK[XLMOD:02001]SEFD`     | 9.1 [X]     |
| N | Cross-linkers (intrachain) | `EVTK[X:DSS#XL1]LEK[#XL1]SEFD`               | 9.2.1 [X]   |
| N | Cross-linkers (interchain) | `EVTK[X:DSS#XL1]L//EK[#XL1]SEFD`             | 9.2.2 [X]   |
| N | Branches                   | `ED[MOD:00093#BRANCH]//D[#BRANCH]ATR`        | 9.3 [X]     |
| Y | GNO modifications          | `NEEYN[GNO:G59626AS]K`                       | 10.1 [G]    |
| Y | Glycan compositions        | `NEEYN[Glycan:Hex5HexNAc4NeuAc1]K`           | 10.2 [G]    |
| Y | Mixed glycan components    | `N[Glycan:Hex{H2O}{+204.068}]K`             | 10.2 [G]   |
| Y | Charged formulas           | `SEQUEN[Formula:Zn1:z+2]CE`                  | 11.1 [3]    |
| P | Controlling placement      | `PTI(MERMERME)[+32|Position:E]PTIDE`         | 11.2 [3]    |
| Y | Global isotope             | `<13C>PEPTIDE`                                | 11.3.1 [3]  |
| Y | Fixed modifications        | `<[Oxidation]@M>ATPEMILTCMGCLK`              | 11.3.2 [3]  |
| Y | Chimeric spectra           | `NEEYN+SEQUEN`                               | 11.4 [3]    |
| Y | Charges                    | `SEQUEN/2`, `SEQUEN/[Na:z+1,H:z+1]`          | 11.5 [3]    |
| Y | Ion notation               | `SEQUEN-[b-type-ion]`                        | 11.6 [3]    |

**Table 1** summarizes notation handling against the finalized ProForma 2.1 specification [@proforma-2026]. In column S, Y indicates supported notation, P indicates a preserved annotation whose placement constraints are not applied, and N indicates unsupported linked-chain calculation semantics. Levels 1, 2, and 3 follow the specification, with top-down (T), cross-linking (X), and glycan (G) extensions. This table is not a claim that every represented sequence supports every calculation.

Mass-ambiguous residues B and Z cannot produce a unique mass. Delta-mass tags and bare-mass glycan components lack elemental compositions. Unknown localization and ambiguous intervals restrict fragmentation. Labile modifications contribute to precursor mass but are omitted from fragment ions. Chimeric assignments are handled component-wise through `parse_chimeric()`. The structured model and JSON format can represent linked notation, but the calculation APIs do not implement cross-links or branches. The documentation provides operation-specific limits and examples.

# AI usage disclosure

Opus 5 and Fable 5.1 through Claude Code, and Sol and Astra through OpenAI Codex, assisted with code generation, refactoring, tests, debugging, documentation, and manuscript revision. Verification included automated tests, static analysis, execution of manuscript examples, and checks against primary references. The human authors reviewed, edited, and validated all AI-assisted work and made the core design decisions.

# Availability

Peptacular is distributed through PyPI (<https://pypi.org/project/peptacular/>) and available as open-source software on GitHub (<https://github.com/tacular-omics/peptacular>). Documentation is accessible at <https://peptacular.readthedocs.io>. The software is released under the MIT license.

# Acknowledgements

This work was supported by the National Institutes of Health under grants R01 AG077046 (Analysis of protein interactions in neurodegenerative disease), R01 MH132570 (Brain-wide mapping of neuronal inhibition by novel inverse activity markers), R01 MH100175 (Proteogenetics of Autism Spectrum Disorders), R01 HL165168 (The CFTR Interactome), and U01 AG088679 (Understanding Gene-Environment Interactions in Brain Aging and Alzheimer's Disease (AD) and AD-Related Dementias (ADRD)).

The funders provided financial support only.

# Competing interests

The authors report no competing interests.

# References
