"""Optional interoperability helpers for third-party proteomics packages.

Integration modules import their third-party dependency only when a conversion
function is called. Importing :mod:`peptacular.interop` therefore never requires
an optional dependency.
"""

from ._errors import (
    InteropConversionError,
    InteropError,
    LossyConversionWarning,
    MissingOptionalDependencyError,
)
from ._policy import LossPolicy
from .alphabase import AlphaBasePeptide, from_alphabase, to_alphabase
from .psm_utils import from_psm_utils, to_psm_utils
from .pyteomics import (
    from_pyteomics,
    from_pyteomics_composition,
    to_pyteomics,
    to_pyteomics_composition,
)

__all__ = [
    "AlphaBasePeptide",
    "InteropConversionError",
    "InteropError",
    "LossPolicy",
    "LossyConversionWarning",
    "MissingOptionalDependencyError",
    "from_alphabase",
    "from_psm_utils",
    "from_pyteomics",
    "from_pyteomics_composition",
    "to_alphabase",
    "to_psm_utils",
    "to_pyteomics",
    "to_pyteomics_composition",
]
