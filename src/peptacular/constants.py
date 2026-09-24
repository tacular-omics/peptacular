from enum import StrEnum
from typing import Final, Literal

from tacular.constants import C13_C12_MASS_DIFF, ELECTRON_MASS, HYDROGEN_MASS, NEUTRON_MASS, PROTON_MASS

from .diagnostics import PeptacularError

__all__ = [
    "PROTON_MASS",
    "PROTON_CARRIER_MASS",
    "ELECTRON_MASS",
    "NEUTRON_MASS",
    "C13_NEUTRON_MASS",
    "PEPTIDE_AVERAGINE_NEUTRON_MASS",
    "CV",
    "Terminal",
    "ModType",
    "ModTypeLiteral",
    "ParallelMethod",
    "ParallelMethodLiteral",
]


# Physical constants come from tacular.constants (CODATA 2018; 13C-12C from tacular's
# isotope table) so every tacular-omics package uses the same values.
C13_NEUTRON_MASS: Final[float] = C13_C12_MASS_DIFF  # 13C - 12C, the isotope-peak spacing

# The monoisotopic mass one default (protonated) charge adds: a hydrogen atom minus an electron.
# This is 1.4e-8 Da below CODATA PROTON_MASS (the H 1s binding energy). peptacular uses the
# hydrogen-atom form so a charged mass agrees with the ion's elemental composition, which counts
# one H atom per charge. Average masses use the average H mass minus an electron instead.
PROTON_CARRIER_MASS: Final[float] = HYDROGEN_MASS - ELECTRON_MASS
PEPTIDE_AVERAGINE_NEUTRON_MASS: Final[float] = 1.002856


class CV(StrEnum):
    """One of the five supported controlled vocabularies"""

    UNIMOD = "UNIMOD"
    PSI_MOD = "MOD"
    RESID = "RESID"
    GNOME = "GNO"
    XL_MOD = "XLMOD"
    CUSTOM = "CUSTOM"
    OBSERVED = "OBSERVED"


_CV_TO_NAME_PREFIX: Final[dict[CV, str]] = {
    CV.UNIMOD: "",
    CV.PSI_MOD: "",
    CV.RESID: "R:",
    CV.GNOME: "G:",
    CV.XL_MOD: "X:",
    CV.CUSTOM: "C:",
}

_CV_TO_ACCESSION_PREFIX: Final[dict[CV, str]] = {
    CV.UNIMOD: "UNIMOD:",
    CV.PSI_MOD: "MOD:",
    CV.RESID: "RESID:",
    CV.GNOME: "GNO:",
    CV.XL_MOD: "XLMOD:",
}

_CV_TO_MASS_PREFIX: Final[dict[CV, str]] = {
    CV.UNIMOD: "U:",
    CV.PSI_MOD: "M:",
    CV.RESID: "R:",
    CV.GNOME: "G:",
    CV.XL_MOD: "X:",
    CV.CUSTOM: "C:",
    CV.OBSERVED: "Obs:",
}


class Terminal(StrEnum):
    """Terminal position specification"""

    ANYWHERE = "Anywhere"
    N_TERM = "N-term"
    C_TERM = "C-term"

    @classmethod
    def from_str(cls, term: str) -> "Terminal":
        """Get Terminal enum from string"""
        term_upper = term.upper()
        if term_upper == "N-TERM":
            return cls.N_TERM
        elif term_upper == "C-TERM":
            return cls.C_TERM
        raise PeptacularError(f"Unknown terminal type: {term}")


class ModType(StrEnum):
    """The kinds of modification a ProForma annotation stores.

    Used as keys by ``get_mods``, ``set_mods``, ``append_mods``, ``extend_mods`` and
    ``remove_mods``/``strip_mods``; the string values (``"nterm"``, ``"internal"``, ...)
    are accepted wherever a ``ModType`` is.

    - ``NTERM``/``CTERM``: terminal modifications (``[Acetyl]-PEP``, ``PEP-[Amidated]``)
    - ``INTERNAL``: residue modifications (``PEM[Oxidation]``)
    - ``ISOTOPE``: global isotope labels (``<13C>``)
    - ``STATIC``: fixed modifications (``<[Carbamidomethyl]@C>``)
    - ``LABILE``: labile modifications (``{Glycan:Hex}``)
    - ``UNKNOWN``: modifications of unknown position (``[Phospho]?PEP``)
    - ``INTERVAL``: modifications on a residue range (``P(EP)[Phospho]TIDE``)
    - ``CHARGE``: the charge state or adducts (``/2``, ``/[Na:z+1]``)
    """

    NTERM = "nterm"
    CTERM = "cterm"
    ISOTOPE = "isotope"
    STATIC = "static"
    LABILE = "labile"
    UNKNOWN = "unknown"
    INTERVAL = "interval"
    INTERNAL = "internal"
    CHARGE = "charge"


ModTypeLiteral = Literal[
    "nterm",
    "cterm",
    "isotope",
    "static",
    "labile",
    "unknown",
    "interval",
    "internal",
    "charge",
]


class ParallelMethod(StrEnum):
    """Backend used by the functional API for list input (the ``method=`` keyword).

    - ``PROCESS``: a ``multiprocessing`` process pool
    - ``THREAD``: a thread pool (the default on free-threaded Python)
    - ``SEQUENTIAL``: no pool, items run one after another
    """

    PROCESS = "process"
    THREAD = "thread"
    SEQUENTIAL = "sequential"


ParallelMethodLiteral = Literal["process", "thread", "sequential"]
