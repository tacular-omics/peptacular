from enum import StrEnum
from typing import Final, Literal

PROTON_MASS: Final[float] = 1.00727646688
ELECTRON_MASS: Final[float] = 0.00054857990946
NEUTRON_MASS: Final[float] = 1.00866491597
C13_NEUTRON_MASS: Final[float] = 1.003350
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


CV_TO_NAME_PREFIX: Final[dict[CV, str]] = {
    CV.UNIMOD: "",
    CV.PSI_MOD: "",
    CV.RESID: "R:",
    CV.GNOME: "G:",
    CV.XL_MOD: "X:",
    CV.CUSTOM: "C:",
}

CV_TO_ACCESSION_PREFIX: Final[dict[CV, str]] = {
    CV.UNIMOD: "UNIMOD:",
    CV.PSI_MOD: "MOD:",
    CV.RESID: "RESID:",
    CV.GNOME: "GNO:",
    CV.XL_MOD: "XLMOD:",
}

CV_TO_MASS_PREFIX: Final[dict[CV, str]] = {
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
        raise ValueError(f"Unknown terminal type: {term}")


class ModType(StrEnum):
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


class parallelMethod(StrEnum):
    PROCESS = "process"
    THREAD = "thread"
    SEQUENTIAL = "sequential"


parallelMethodLiteral = Literal["process", "thread", "sequential"]
