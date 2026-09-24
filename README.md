# Peptacular

<div align="center">
  <img src="https://raw.githubusercontent.com/tacular-omics/peptacular/main/peptacular_logo.png" alt="Peptacular Logo" width="400" style="margin: 20px;"/>

[![Python package](https://github.com/tacular-omics/peptacular/actions/workflows/ci.yml/badge.svg)](https://github.com/tacular-omics/peptacular/actions/workflows/ci.yml)
[![codecov](https://codecov.io/github/tacular-omics/peptacular/graph/badge.svg?token=1CTVZVFXF7)](https://codecov.io/github/tacular-omics/peptacular)
[![PyPI version](https://badge.fury.io/py/peptacular.svg)](https://badge.fury.io/py/peptacular)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15054278.svg)](https://doi.org/10.5281/zenodo.15054278)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-g.svg)](https://opensource.org/licenses/MIT)

</div>

Peptacular parses [ProForma 2.1](https://github.com/HUPO-PSI/ProForma) peptide
sequences and calculates their masses, fragments, and isotopic distributions.
It's for anyone working with peptide-level proteomics data in Python who wants
exact masses and fragment ions without hand-rolling ProForma parsing and mass
tables. It's built on [tacular](https://github.com/tacular-omics/tacular)'s
lookup data, and its fragments export directly as mzPAF strings readable by
[paftacular](https://github.com/tacular-omics/paftacular).

## Why peptacular?

- **Full ProForma 2.1 parsing** into a chainable, editable `ProFormaAnnotation`
  object — or use the functional API directly on strings.
- **Mass, m/z, composition, and predicted isotopic distributions**, with
  monoisotopic and average mass support.
- **Enzymatic digestion** with missed cleavages, semi-specific, and
  non-specific modes.
- **Fragment ion generation** for 20+ ion types, exportable straight to
  mzPAF strings for [paftacular](https://github.com/tacular-omics/paftacular).
- **Batch-friendly**: functional API calls on lists of sequences parallelize
  automatically, with per-item error collection. Pass FASTA/PEFF entries from
  [fastatacular](https://github.com/tacular-omics/fastatacular) straight in: any object with a
  `sequence` string works.
- **Type-annotated throughout**, plus optional Pyteomics, psm_utils, AlphaBase,
  and MCP integrations.

## Install

```bash
pip install peptacular
```

Optional integrations install as extras:

```bash
pip install "peptacular[pyteomics]"
pip install "peptacular[psm-utils]"
pip install "peptacular[alphabase]"
pip install "peptacular[mcp]"
pip install "peptacular[numpy]"   # fragment_arrays(): ions as numpy columns
```

See the [interoperability guide](https://peptacular.readthedocs.io/en/latest/interoperability.html)
for supported conversions.

## Quick example

```python
import peptacular as pt

# Parse a sequence into a ProFormaAnnotation
peptide = pt.parse("PEM[Oxidation]TIDE")

# Calculate mass and m/z
print(peptide.mass())              # 849.3426002717299
print(peptide.mz(charge=2))        # 425.6785766024859

# Chained edits return a modified annotation
print(peptide.set_charge(2).set_peptide_name("Peptacular").serialize())
# (>Peptacular)PEM[Oxidation]TIDE/2
```

## What else it can do

Digest a protein and generate fragment ions that round-trip through
[paftacular](https://github.com/tacular-omics/paftacular)'s mzPAF parser:

```python
import peptacular as pt

peptides = pt.digest("MKVLATSAGERTIDEK", enzyme="trypsin", missed_cleavages=1)
print([seq for seq, _ in peptides])
# ['MK', 'MKVLATSAGER', 'VLATSAGER', 'VLATSAGERTIDEK', 'TIDEK']

fragments = pt.fragment("PEPTIDE", ion_types=("b", "y"), charges=[1])
print(fragments[1].to_mzpaf())  # b2{PE}
```

The functional API operates on lists directly, auto-parallelizing for larger
batches:

```python
import peptacular as pt

peptides = ["[Acetyl]-PEPTIDES", "<13C>ARE", "SICK/2"]
print(pt.mass(peptides))               # [928.4025574375299, 388.23835027296, 451.24535797517194]
print(pt.mz(peptides, charge=2))       # [465.20855518538593, 195.12645160310103, 225.62267898758597]
```

For streaming input and per-item error collection instead of a raised
exception, see the [streaming guide](https://peptacular.readthedocs.io/en/latest/streaming.html):

```python
import peptacular as pt

results = pt.batch("mass", ["PEPTIDE", "PEP[UnknownModification]TIDE"], errors="collect")
print(results[0].value)                # 799.3599640328299
print(results[1].error.code)           # unresolved_modification
```

Raised errors are typed and all subclass `pt.PeptacularError` (a `ValueError`):
invalid ProForma raises `pt.ProFormaFormatError`, an unresolved modification
`pt.UnknownModificationError`, and so on. See the
[streaming guide](https://peptacular.readthedocs.io/en/latest/streaming.html) for the full list.

Expand an ambiguous modification into its localization isomers and find the fragment
ions that tell them apart. Candidate sites come from the ProForma string alone; peptacular
has no built-in list of which residues a mod can sit on:

```python
import peptacular as pt

isomers = pt.localization_isomers("PEP(ST)[Phospho]IDE")
print([a.serialize() for a in isomers])  # ['PEPS[Phospho]TIDE', 'PEPST[Phospho]IDE']
ions = pt.site_determining_ions(isomers, ion_types=("b", "y"), charges=(1,))
print([[f"{f.ion_type}{f.position}" for f in frags] for frags in ions])  # [['b4', 'y4'], ['b4', 'y4']]
```

| Area | Entry points |
|---|---|
| Digestion | `pt.digest`, `pt.semi_digest`, `pt.nonspecific_digest` |
| Fragmentation | `pt.fragment`, `pt.fast_fragment` |
| Localization | `pt.localization_isomers`, `pt.candidate_sites`, `pt.site_determining_ions`, `pt.pairwise_site_determining_ions` ([guide](https://peptacular.readthedocs.io/en/latest/localization.html)) |
| Isotopes | `pt.isotopic_distribution`, `pt.brain_isotopic_distribution` |
| Tables | `pt.digest_records`, `pt.fragment_records` (plain dicts for pandas or polars), `pt.fragment_arrays` (numpy columns) |
| Batch / streaming | `pt.batch`, `pt.iter_batch`, `pt.diagnose` (read FASTA with fastatacular) |
| JSON interchange | see the [JSON serialization guide](https://peptacular.readthedocs.io/en/latest/json_serialization.html) |

## Tables with pandas or polars

peptacular does not ship pandas or polars. `pt.digest_records` and `pt.fragment_records`
return a list of plain dicts (strings, numbers, booleans, `None`), one per peptide or ion,
which either library turns into a table. Column names are listed in
`pt.DIGEST_RECORD_KEYS` and `pt.FRAGMENT_RECORD_KEYS`:

```python
import peptacular as pt

rows = pt.digest_records("MKVLATSAGERTIDEK", "trypsin", missed_cleavages=1)
print(rows[0])  # {'peptide': 'MK', 'stripped_sequence': 'MK', 'start': 0, 'end': 2, 'missed_cleavages': 0, 'semi': False, 'accession': None}
ions = pt.fragment_records(pt.fragment("PEPTIDE/2", ion_types=("b", "y"), charges=(1, 2)))
print(ions[1]["ion_type"], ions[1]["position"], ions[1]["charge_state"], ions[1]["mzpaf"])  # b 2 1 b2{PE}
# pandas.DataFrame(rows) or polars.DataFrame(ions) gives a table
```

A FASTA entry's `accession` (or a PEFF entry's `db_unique_id`) is copied into each digest row.

For many peptides, `pt.fragment_arrays` (needs `pip install "peptacular[numpy]"`) returns the
same ions as `pt.fragment` as a dict of numpy columns, one row per ion, with a `peptide_index`
column. It is about 10x faster than building a table from `Fragment` objects:

```python
cols = pt.fragment_arrays(["PEPTIDE/2", "PEM[Oxidation]K"], ion_types=("b", "y"), charges=(1, 2))
print(cols["peptide_index"][:3].tolist(), cols["mz"][:3].round(4).tolist())  # [0, 0, 0] [98.06, 227.1026, 324.1554]
# polars.DataFrame(cols) or pyarrow.table(cols) takes the dict as is
```

See the [tables guide](https://peptacular.readthedocs.io/en/latest/records.html).

## Local MCP integration

Peptacular includes 12 optional MCP tools for agents to inspect annotations,
calculate theoretical properties, digest protein sequences, and transform
annotations. Calls accept small inline batches and return results directly,
with no stored data or job setup. Install with `pip install "peptacular[mcp]"`,
then check the installation:

```text
peptacular-mcp --check
```

See the [local MCP guide](https://github.com/tacular-omics/peptacular/blob/main/docs/mcp.rst)
for client setup, tool examples, and limits.

## Documentation

- Full docs: [peptacular.readthedocs.io](https://peptacular.readthedocs.io/en/latest/index.html)
- Changelog: [CHANGELOG.md](https://github.com/tacular-omics/peptacular/blob/main/CHANGELOG.md)
- ProForma 2.1 compliance status: [PROFORMA_COMPLIANCE.md](https://github.com/tacular-omics/peptacular/blob/main/PROFORMA_COMPLIANCE.md)
- Contributing: [CONTRIBUTING.md](https://github.com/tacular-omics/peptacular/blob/main/CONTRIBUTING.md)

## License

MIT

## Citation

Working on a JOSS submission, but in the meantime use:

https://doi.org/10.5281/zenodo.15054278
