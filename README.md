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
  automatically, with streaming FASTA/gzip input and per-item error collection.
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
print(peptide.mz(charge=2))        # 425.67857658818554

# Chained edits return a modified annotation
print(peptide.set_charge(2).set_peptide_name("Peptacular").serialize())
# (>Peptacular)PEM[Oxidation]TIDE/2
```

## What else it can do

Digest a protein and generate fragment ions that round-trip through
[paftacular](https://github.com/tacular-omics/paftacular)'s mzPAF parser:

```python
import peptacular as pt

trypsin = pt.PROTEASE_LOOKUP["trypsin"]
peptides = pt.digest("MKVLATSAGERTIDEK", enzyme_regex=trypsin.regex, missed_cleavages=1)
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
print(pt.mass(peptides))               # [928.4025574375299, 388.23835027296, 451.245357946571]
print(pt.mz(peptides, charge=2))       # [465.20855517108555, 195.12645158880056, 225.6226789732855]
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

| Area | Entry points |
|---|---|
| Digestion | `pt.digest`, `pt.semi_digest`, `pt.nonspecific_digest` |
| Fragmentation | `pt.fragment`, `pt.fast_fragment` |
| Isotopes | `pt.isotopic_distribution`, `pt.brain_isotopic_distribution` |
| FASTA / streaming | `pt.parse_fasta`, `pt.iter_fasta`, `pt.batch`, `pt.iter_batch` |
| JSON interchange | see the [JSON serialization guide](https://peptacular.readthedocs.io/en/latest/json_serialization.html) |

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
