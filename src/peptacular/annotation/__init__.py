"""
The annotation module contains classes and functions for handling peptide annotations.

All methods/properties are available on the ProFormaAnnotation class. The one exception is the AnnotationProperties,
which is accessible via the annot.prop attribute.
"""

from .annotation import (
    AnnotationProperties,
    CrossLinkFragmentMode,
    MultiProFormaAnnotation,
    ProFormaAnnotation,
    validate_crosslink_labels,
)
from .mod import Mod, Mods
from .parser import Interval
from .utils import Fragment

__all__ = [
    "Interval",
    "ProFormaAnnotation",
    "MultiProFormaAnnotation",
    "CrossLinkFragmentMode",
    "AnnotationProperties",
    "validate_crosslink_labels",
    "Mod",
    "Mods",
    "Fragment",
]
