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
from .alphabase import (
    AlphaBaseRow,
    from_alphabase_dataframe,
    from_alphabase_row,
    to_alphabase_dataframe,
    to_alphabase_row,
)
from .psm_utils import from_psm_utils, to_psm_utils
from .pyteomics import (
    from_pyteomics,
    from_pyteomics_composition,
    to_pyteomics,
    to_pyteomics_composition,
)

__all__ = [
    "AlphaBaseRow",
    "InteropConversionError",
    "InteropError",
    "LossPolicy",
    "LossyConversionWarning",
    "MissingOptionalDependencyError",
    "from_alphabase_dataframe",
    "from_alphabase_row",
    "from_psm_utils",
    "from_pyteomics",
    "from_pyteomics_composition",
    "to_alphabase_dataframe",
    "to_alphabase_row",
    "to_psm_utils",
    "to_pyteomics",
    "to_pyteomics_composition",
]
