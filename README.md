# Peptacular

<div align="center">
  <img src="peptacular_logo.png" alt="Peptacular Logo" width="400" style="margin: 20px;"/>
  
  
  A Python package for peptide sequence analysis built around **ProForma 2.1 notation**. Calculate masses, generate fragments, predict isotopic patterns, and more. Peptacular uses type annotations extensively, so it is type safe.
  
  [![Python package](https://github.com/tacular-omics/peptacular/actions/workflows/python-package.yml/badge.svg)](https://github.com/tacular-omics/peptacular/actions/workflows/python-package.yml)
  [![codecov](https://codecov.io/github/tacular-omics/peptacular/graph/badge.svg?token=1CTVZVFXF7)](https://codecov.io/github/tacular-omics/peptacular)
  [![PyPI version](https://badge.fury.io/py/peptacular.svg)](https://badge.fury.io/py/peptacular)
  [![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
  [![License: MIT](https://img.shields.io/badge/License-MIT-g.svg)](https://opensource.org/licenses/MIT)
  
</div>

## Features

- **Nearly Complete ProForma 2.1 Parsing**
- **Modifiable ProFormaAnnotation Objects (Factory Pattern)**
- **Mass/Mz/Composition Calculations**
- **Predicted Isotopic Distributions**
- **Enzymatic Protein Digestion** 
- **Fragment Ion Generation** 
- **Physiochemical Property Calculations** 
- **Built-in Parallel Processing** 

## Installation

```bash
pip install peptacular
```

## Quick Start (Object Based)

See docs for more detail.

```python
import peptacular as pt

# Parse a sequence into a ProFormaAnnotation
peptide: pt.ProFormaAnnotation = pt.parse("PEM[Oxidation]TIDE")

# Calculate mass and m/z
mass: float = peptide.mass() # 849.342
mz: float = peptide.mz(charge=2) # 425.678

# Factory pattern
print(peptide.set_charge(2).set_peptide_name("Peptacular").serialize())
# (>Peptacular)PEM[Oxidation]TIDE/2
```


## Quick Start (Functional Based)

When more than one item is passed to the functional API methods, it is automatically parallelized.

```python
import peptacular as pt

peptides = ['[Acetyl]-PEPTIDES', '<C13>ARE', 'SICK/2']

# Calculate mass and m/z for all peptides
masses: list[float] = pt.mass(peptides) # [928.4026, 374.1914, 451.2454]
mzs: list[float] = pt.mz(peptides, charge=2) # [465.2086, 188.103, 225.6227]
```


## ProForma 2.1 Compliance

See [PROFORMA_COMPLIANCE.md](PROFORMA_COMPLIANCE.md) for detailed compliance status.

## Contributing

Contributions welcome! Check the examples directory for code style and documentation patterns.

## License

MIT

## Citation

Working on a JOSS submission, but in the meantime use:

https://doi.org/10.5281/zenodo.15054278

