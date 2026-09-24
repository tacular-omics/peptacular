# Copyright 2003 Yair Benita.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Indices to be used with ProtParam."""

from enum import StrEnum
from typing import Final

__all__ = [
    "NEGATIVE_AMINO_ACIDS",
    "POSITIVE_AMINO_ACIDS",
    "CHARGED_AMINO_ACIDS",
    "AROMATIC_AMINO_ACIDS",
    "hphob_kyte_doolittle",
    "hphob_adoberin",
    "hphob_abraham_leo",
    "hphob_agros",
    "hphob_rao_argos",
    "hphob_black_mould",
    "hphob_bull_breese",
    "hphob_casari_sippl",
    "hphob_cid",
    "hphob_cowan_3_4",
    "hphob_cowan_7_5",
    "hphob_eisenberg",
    "hphob_engelman",
    "hphob_fasman",
    "hphob_fauchere",
    "hphob_goldsack",
    "hphob_guy",
    "hphob_jones",
    "hphob_juretic",
    "hphob_kidera",
    "hphob_miyazawa",
    "hphob_parker",
    "hphob_ponnuswamy",
    "hphob_manavalan",
    "hphob_rose",
    "hphob_roseman",
    "hphob_sweet",
    "hphob_tanford",
    "hphob_wilson",
    "hphob_zimmerman",
    "hphob_chothia",
    "hphob_janin",
    "hphob_wolfenden",
    "hphob_welling",
    "flexibility_vihinen",
    "hydrophilicity_hopp_wood",
    "surface_accessibility_vergoten",
    "surface_accessiblility_janin",
    "DIWV",
    "pk_nterminal",
    "pk_cterminal",
    "pk_sidechain",
    "deleage_roux_alpha_helix",
    "deleage_roux_beta_sheet",
    "deleage_roux_beta_turn",
    "deleage_roux_coil",
    "levitt_alpha_helix",
    "levitt_beta_sheet",
    "levitt_beta_turn",
    "chou_fasman_alpha_helix",
    "chou_fasman_beta_sheet",
    "chou_fasman_beta_turn",
    "beta_strand_parallel",
    "beta_strand_antiparallel",
    "beta_strand_total",
    "ratioside",
    "polarity_grantham",
    "polarity_zimmerman",
    "mutability",
    "codons",
    "refractivity",
    "bulkiness",
    "recognition_factors",
    "aa_composition_mccaldron",
    "aa_composition_swissprot",
    "aa_composition_scales",
    "transmembrane_tendency",
    "accessible_residues",
    "average_buried_area",
    "molecular_weights",
    "hplc_meek_2_1",
    "hplc_browne",
    "hplc_meek_7_4",
    "hplc_browne_tfa",
    "SecondaryStructureMethod",
    "SecondaryStructureType",
    "PropertyScale",
    "HydrophobicityScale",
    "SecondaryStructureScale",
    "SurfaceAccessibilityScale",
    "ChargeScale",
    "PolarityScale",
    "HPLCScale",
    "BetaStrandScale",
    "PhysicalPropertyScale",
    "CompositionScale",
    "secondary_structure_scales",
    "secondary_structure_scales_by_name",
    "PROPERTY_SCALES",
    "HYDROPHOBICITY_SCALES",
    "SURFACE_ACCESSIBILITY_SCALES",
    "HPLC_SCALES",
    "HYDROPHILICITY_SCALES",
    "FLEXIBILITY_SCALES",
    "POLARITY_SCALES",
    "COMPOSITION_SCALES",
    "PHYSICAL_PROPERTY_SCALES",
]


class _AA(StrEnum):
    A = "A"
    R = "R"
    N = "N"
    D = "D"
    C = "C"
    Q = "Q"
    E = "E"
    G = "G"
    H = "H"
    I = "I"
    L = "L"
    K = "K"
    M = "M"
    F = "F"
    P = "P"
    S = "S"
    T = "T"
    W = "W"
    Y = "Y"
    V = "V"


# Turn black code style off
# fmt: off

# Hydrophobicity

# Kyte & Doolittle index of hydrophobicity
# J. Mol. Biol. 157:105-132(1982).
# "KyteDoolittle"

NEGATIVE_AMINO_ACIDS: Final[set[str]] = {_AA.D, _AA.E, _AA.C, _AA.Y}
POSITIVE_AMINO_ACIDS: Final[set[str]] = {_AA.R, _AA.H, _AA.K}
CHARGED_AMINO_ACIDS: Final[set[str]] = {_AA.K, _AA.R, _AA.H, _AA.D, _AA.E, _AA.C, _AA.Y}
AROMATIC_AMINO_ACIDS: Final[set[str]] = {_AA.F, _AA.W, _AA.Y}

hphob_kyte_doolittle: Final[dict[str, float]] = {_AA.A: 1.8, _AA.R: -4.5, _AA.N: -3.5, _AA.D: -3.5, _AA.C: 2.5,
      _AA.Q: -3.5, _AA.E: -3.5, _AA.G: -0.4, _AA.H: -3.2, _AA.I: 4.5,
      _AA.L: 3.8, _AA.K: -3.9, _AA.M: 1.9, _AA.F: 2.8, _AA.P: -1.6,
      _AA.S: -0.8, _AA.T: -0.7, _AA.W: -0.9, _AA.Y: -1.3, _AA.V: 4.2}

# Aboderin hydrophobicity index
# International J. of Biochemistry, 2(11), 537-544.
# "Aboderin"
hphob_adoberin: Final[dict[str, float]] = {_AA.A: 5.1, _AA.R: 2.0, _AA.N: 0.6, _AA.D: 0.7, _AA.C: 0.0,
      _AA.Q: 1.4, _AA.E: 1.8, _AA.G: 4.1, _AA.H: 1.6, _AA.I: 9.3,
      _AA.L: 10.0, _AA.K: 1.3, _AA.M: 8.7, _AA.F: 9.6, _AA.P: 4.9,
      _AA.S: 3.1, _AA.T: 3.5, _AA.W: 9.2, _AA.Y: 8.0, _AA.V: 8.5}

# Abraham & Leo hydrophobicity index
# Proteins: Structure, Function and Genetics 2:130-152(1987).
# "AbrahamLeo"
hphob_abraham_leo: Final[dict[str, float]] = {_AA.A: 0.44, _AA.R: -2.42, _AA.N: -1.32, _AA.D: -0.31, _AA.C: 0.58,
      _AA.Q: -0.71, _AA.E: -0.34, _AA.G: 0.0, _AA.H: -0.01, _AA.I: 2.46,
      _AA.L: 2.46, _AA.K: -2.45, _AA.M: 1.1, _AA.F: 2.54, _AA.P: 1.29,
      _AA.S: -0.84, _AA.T: -0.41, _AA.W: 2.56, _AA.Y: 1.63, _AA.V: 1.73}

# Argos hydrophobicity index
# European Journal of Biochemistry, 128(2-3), 565-575.
# "Argos"
hphob_agros: Final[dict[str, float]] = {_AA.A: 0.61, _AA.R: 0.6, _AA.N: 0.06, _AA.D: 0.46, _AA.C: 1.07,
      _AA.Q: 0.0, _AA.E: 0.47, _AA.G: 0.07, _AA.H: 0.61, _AA.I: 2.22,
      _AA.L: 1.53, _AA.K: 1.15, _AA.M: 1.18, _AA.F: 2.02, _AA.P: 1.95,
      _AA.S: 0.05, _AA.T: 0.05, _AA.W: 2.65, _AA.Y: 1.88, _AA.V: 1.32}

"""
Amino acid scale: Membrane buried helix parameter.
Author(s):
Rao M.J.K., Argos P.
Reference: Biochim. Biophys. Acta 869:197-214(1986).
https://web.expasy.org/protscale/pscale/Hphob.Argos.html
"""
hphob_rao_argos: Final[dict[str, float]] = {_AA.A: 1.360, _AA.R: 0.150, _AA.N: 0.330, _AA.D: 0.110, _AA.C: 1.270,
      _AA.Q: 0.330, _AA.E: 0.250, _AA.G: 1.090, _AA.H: 0.680, _AA.I: 1.440,
      _AA.L: 1.470, _AA.K: 0.090, _AA.M: 1.420, _AA.F: 1.570, _AA.P: 0.540,
      _AA.S: 0.970, _AA.T: 1.080, _AA.W: 1.000, _AA.Y: 0.830, _AA.V: 1.370}

# Black & Mould hydrophobicity index
# Anal. Biochem. 193:72-82(1991).
# "BlackMould"
hphob_black_mould: Final[dict[str, float]] = {_AA.A: 0.616, _AA.R: 0.0, _AA.N: 0.236, _AA.D: 0.028, _AA.C: 0.68,
      _AA.Q: 0.251, _AA.E: 0.043, _AA.G: 0.501, _AA.H: 0.165, _AA.I: 0.943,
      _AA.L: 0.943, _AA.K: 0.283, _AA.M: 0.738, _AA.F: 1.0, _AA.P: 0.711,
      _AA.S: 0.359, _AA.T: 0.45, _AA.W: 0.878, _AA.Y: 0.88, _AA.V: 0.825}

# Bull & Breese hydrophobicity index
# Arch. Biochem. Biophys. 161:665-670(1974)
# "BullBreese"
hphob_bull_breese: Final[dict[str, float]] = {_AA.A: 0.61, _AA.R: 0.69, _AA.N: 0.89, _AA.D: 0.61, _AA.C: 0.36,
      _AA.Q: 0.97, _AA.E: 0.51, _AA.G: 0.81, _AA.H: 0.69, _AA.I: -1.45,
      _AA.L: -1.65, _AA.K: 0.46, _AA.M: -0.66, _AA.F: -1.52, _AA.P: -0.17,
      _AA.S: 0.42, _AA.T: 0.29, _AA.W: -1.2, _AA.Y: -1.43, _AA.V: -0.75}

# Casari & Sippl hydrophobic potential
# Journal of molecular biology, 224(3), 725-732.
# "Casari"
hphob_casari_sippl: Final[dict[str, float]] = {_AA.A: 0.2, _AA.R: -0.7, _AA.N: -0.5, _AA.D: -1.4, _AA.C: 1.9,
      _AA.Q: -1.1, _AA.E: -1.3, _AA.G: -0.1, _AA.H: 0.4, _AA.I: 1.4,
      _AA.L: 0.5, _AA.K: -1.6, _AA.M: 0.5, _AA.F: 1.0, _AA.P: -1.0,
      _AA.S: -0.7, _AA.T: -0.4, _AA.W: 1.6, _AA.Y: 0.5, _AA.V: 0.7}

# Cid hydrophobicity index
# Protein engineering, 5(5), 373-375.
# "Cid"
hphob_cid: Final[dict[str, float]] = {_AA.A: 0.02, _AA.R: -0.42, _AA.N: -0.77, _AA.D: -1.04, _AA.C: 0.77,
      _AA.Q: -1.1, _AA.E: -1.14, _AA.G: -0.8, _AA.H: 0.26, _AA.I: 1.81,
      _AA.L: 1.14, _AA.K: -0.41, _AA.M: 1.0, _AA.F: 1.35, _AA.P: -0.09,
      _AA.S: -0.97, _AA.T: -0.77, _AA.W: 1.71, _AA.Y: 1.11, _AA.V: 1.13}

# Cowan hydrophobicity indices at ph 3.4 and 7.5
# Peptide Research 3:75-80(1990).
# "Cowan3.4" "Conwan7.5"
_hphob_cowan: Final[dict[float, dict[str, float]]] = {3.4 : {_AA.A: 0.42, _AA.R: -1.56, _AA.N: -1.03, _AA.D: -0.51, _AA.C: 0.84,
             _AA.Q: -0.96, _AA.E: -0.37, _AA.G: 0.0, _AA.H: -2.28, _AA.I: 1.81,
             _AA.L: 1.8, _AA.K: -2.03, _AA.M: 1.18, _AA.F: 1.74, _AA.P: 0.86,
             _AA.S: -0.64, _AA.T: -0.26, _AA.W: 1.46, _AA.Y: 0.51, _AA.V: 1.34},
      7.5 : {_AA.A: 0.35, _AA.R: -1.5, _AA.N: -0.99, _AA.D: -2.15, _AA.C: 0.76,
             _AA.Q: -0.93, _AA.E: -1.95, _AA.G: 0.0, _AA.H: -0.65, _AA.I: 1.83,
             _AA.L: 1.8, _AA.K: -1.54, _AA.M: 1.1, _AA.F: 1.69, _AA.P: 0.84,
             _AA.S: -0.63, _AA.T: -0.27, _AA.W: 1.35, _AA.Y: 0.39, _AA.V: 1.32}
      }

hphob_cowan_3_4: Final[dict[str, int | float]] = _hphob_cowan[3.4]
hphob_cowan_7_5: Final[dict[str, int | float]] = _hphob_cowan[7.5]

# Eisenberg Normalized consensus hydrophobicity scale
# J. Mol. Biol. 179:125-142(1984)
# "Eisenberg"
hphob_eisenberg: Final[dict[str, float]] = {_AA.A: 0.62, _AA.R: -2.53, _AA.N: -0.78, _AA.D: -0.9, _AA.C: 0.29,
      _AA.Q: -0.85, _AA.E: -0.74, _AA.G: 0.48, _AA.H: -0.4, _AA.I: 1.38,
      _AA.L: 1.06, _AA.K: -1.5, _AA.M: 0.64, _AA.F: 1.19, _AA.P: 0.12,
      _AA.S: -0.18, _AA.T: -0.05, _AA.W: 0.81, _AA.Y: 0.26, _AA.V: 1.08}

# Engelman Hydrophobic Transfer Free Energies
# Annual review of biophysics and biophysical chemistry, 15(1), 321-353.
# "Engelman"
hphob_engelman: Final[dict[str, float]] = {_AA.A: -1.6, _AA.R: 12.3, _AA.N: 4.8, _AA.D: 9.2, _AA.C: -2,
      _AA.Q: 4.1, _AA.E: 8.2, _AA.G: -1, _AA.H: 3, _AA.I: -3.1,
      _AA.L: -2.8, _AA.K: 8.8, _AA.M: -3.4, _AA.F: -3.7, _AA.P: 0.2,
      _AA.S: -0.6, _AA.T: -1.2, _AA.W: -1.9, _AA.Y: 0.7, _AA.V: -2.6}

# Fasman hydrophobicity index
# (1989). Prediction of protein structure and the principles of protein conformation. Springer.
# "Fasman"
hphob_fasman: Final[dict[str, float]] = {_AA.A: -0.21, _AA.R: 2.11, _AA.N: 0.96, _AA.D: 1.36, _AA.C: -6.04,
      _AA.Q: 1.52, _AA.E: 2.3, _AA.G: 0, _AA.H: -1.23, _AA.I: -4.81,
      _AA.L: -4.68, _AA.K: 3.88, _AA.M: -3.66, _AA.F: -4.65, _AA.P: 0.75,
      _AA.S: 1.74, _AA.T: 0.78, _AA.W: -3.32, _AA.Y: -1.01, _AA.V: -3.5}

# Fauchere Hydrophobicity scale
# Eur. J. Med. Chem. 18:369-375(1983).
# "Fauchere"
hphob_fauchere: Final[dict[str, float]] = {_AA.A: 0.31, _AA.R: -1.01, _AA.N: -0.6, _AA.D: -0.77, _AA.C: 1.54,
      _AA.Q: -0.22, _AA.E: -0.64, _AA.G: 0, _AA.H: 0.13, _AA.I: 1.8,
      _AA.L: 1.7, _AA.K: -0.99, _AA.M: 1.23, _AA.F: 1.79, _AA.P: 0.72,
      _AA.S: -0.04, _AA.T: 0.26, _AA.W: 2.25, _AA.Y: 0.96, _AA.V: 1.22}

# Goldsack & Chalifoux Free Energy of Mixing of the Hydrophobic Side Chains
# Journal of theoretical biology, 39(3), 645-651.
# "Goldsack"
hphob_goldsack: Final[dict[str, float]] = {_AA.A: 0.75, _AA.R: 0.75, _AA.N: 0.69, _AA.D: 0, _AA.C: 1,
      _AA.Q: 0.59, _AA.E: 0, _AA.G: 0, _AA.H: 0, _AA.I: 2.95,
      _AA.L: 2.4, _AA.K: 1.5, _AA.M: 1.3, _AA.F: 2.65, _AA.P: 2.6,
      _AA.S: 0, _AA.T: 0.45, _AA.W: 3, _AA.Y: 2.85, _AA.V: 1.7}

# Guy Hydrophobicity scale based on free energy of transfer (kcal/mole).
# Biophys J. 47:61-70(1985)
# "Guy"
hphob_guy: Final[dict[str, float]] = {_AA.A: 0.1, _AA.R: 1.91, _AA.N: 0.48, _AA.D: 0.78, _AA.C: -1.42,
      _AA.Q: 0.95, _AA.E: 0.83, _AA.G: 0.33, _AA.H: -0.5, _AA.I: -1.13,
      _AA.L: -1.18, _AA.K: 1.4, _AA.M: -1.59, _AA.F: -2.12, _AA.P: 0.73,
      _AA.S: 0.52, _AA.T: 0.07, _AA.W: -0.51, _AA.Y: -0.21, _AA.V: -1.27}

# Jones Hydrophobicity scale
# Journal of theoretical biology, 50(1), 167-183.
# "Jones"
hphob_jones: Final[dict[str, float]] = {_AA.A: 0.87, _AA.R: 0.85, _AA.N: 0.09, _AA.D: 0.66, _AA.C: 1.52,
      _AA.Q: 0, _AA.E: 0.67, _AA.G: 0.1, _AA.H: 0.87, _AA.I: 3.15,
      _AA.L: 2.17, _AA.K: 1.64, _AA.M: 1.67, _AA.F: 2.87, _AA.P: 2.77,
      _AA.S: 0.07, _AA.T: 0.07, _AA.W: 3.77, _AA.Y: 2.67, _AA.V: 1.87}

# Juretic Hydrophobicity scale
# Theoretical and computational chemistry, 5, 405-445.
# "Juretic"
hphob_juretic: Final[dict[str, float]] = {_AA.A: 1.1, _AA.R: -5.1, _AA.N: -3.5, _AA.D: -3.6, _AA.C: 2.5,
      _AA.Q: -3.68, _AA.E: -3.2, _AA.G: -0.64, _AA.H: -3.2, _AA.I: 4.5,
      _AA.L: 3.8, _AA.K: -4.11, _AA.M: 1.9, _AA.F: 2.8, _AA.P: -1.9,
      _AA.S: -0.5, _AA.T: -0.7, _AA.W: -0.46, _AA.Y: -1.3, _AA.V: 4.2}

# Kidera Hydrophobicity Factors
# Journal of Protein Chemistry, 4(1), 23-55.
# "Kidera"
hphob_kidera: Final[dict[str, float]] = {_AA.A: -0.27, _AA.R: 1.87, _AA.N: 0.81, _AA.D: 0.81, _AA.C: -1.05,
      _AA.Q: 1.1, _AA.E: 1.17, _AA.G: -0.16, _AA.H: 0.28, _AA.I: -0.77,
      _AA.L: -1.1, _AA.K: 1.7, _AA.M: -0.73, _AA.F: -1.43, _AA.P: -0.75,
      _AA.S: 0.42, _AA.T: 0.63, _AA.W: -1.57, _AA.Y: -0.56, _AA.V: -0.4}

# Miyazawa Hydrophobicity scale (contact energy derived from 3D data)
# Macromolecules 18:534-552(1985)
# "Miyazawa"
hphob_miyazawa: Final[dict[str, float]] = {_AA.A: 5.33, _AA.R: 4.18, _AA.N: 3.71, _AA.D: 3.59, _AA.C: 7.93,
      _AA.Q: 3.87, _AA.E: 3.65, _AA.G: 4.48, _AA.H: 5.1, _AA.I: 8.83,
      _AA.L: 8.47, _AA.K: 2.95, _AA.M: 8.95, _AA.F: 9.03, _AA.P: 3.87,
      _AA.S: 4.09, _AA.T: 4.49, _AA.W: 7.66, _AA.Y: 5.89, _AA.V: 7.63}

# Parker Hydrophilicity scale derived from HPLC peptide retention times
# Biochemistry 25:5425-5431(1986)
# "Parker"
hphob_parker: Final[dict[str, float]] = {_AA.A: 2.1, _AA.R: 4.2, _AA.N: 7, _AA.D: 10, _AA.C: 1.4,
      _AA.Q: 6, _AA.E: 7.8, _AA.G: 5.7, _AA.H: 2.1, _AA.I: -8,
      _AA.L: -9.2, _AA.K: 5.7, _AA.M: -4.2, _AA.F: -9.2, _AA.P: 2.1,
      _AA.S: 6.5, _AA.T: 5.2, _AA.W: -10, _AA.Y: -1.9, _AA.V: -3.7}

# Ponnuswamy Hydrophobic characteristics of folded proteins
# Progress in biophysics and molecular biology, 59(1), 57-103.
# "Ponnuswamy"
hphob_ponnuswamy: Final[dict[str, float]] = {_AA.A: 0.85, _AA.R: 0.2, _AA.N: -0.48, _AA.D: -1.1, _AA.C: 2.1,
      _AA.Q: -0.42, _AA.E: -0.79, _AA.G: 0, _AA.H: 0.22, _AA.I: 3.14,
      _AA.L: 1.99, _AA.K: -1.19, _AA.M: 1.42, _AA.F: 1.69, _AA.P: -1.14,
      _AA.S: -0.52, _AA.T: -0.08, _AA.W: 1.76, _AA.Y: 1.37, _AA.V: 2.53}

"""
Amino acid scale: Average surrounding hydrophobicity.
Author(s):
Manavalan P., Ponnuswamy P.K.
Reference: Nature 275:673-674(1978).
https://web.expasy.org/protscale/pscale/Hphob.Manavalan.html
"""
hphob_manavalan: Final[dict[str, float]] = {_AA.A: 12.97, _AA.R: 11.72, _AA.N: 11.42, _AA.D: 10.85, _AA.C: 14.63,
      _AA.Q: 11.76, _AA.E: 11.89, _AA.G: 12.43, _AA.H: 12.16, _AA.I: 15.67,
      _AA.L: 14.9, _AA.K: 11.36, _AA.M: 14.39, _AA.F: 14.0, _AA.P: 11.37,
      _AA.S: 11.23, _AA.T: 11.69, _AA.W: 13.93, _AA.Y: 13.42, _AA.V: 15.71}

# Rose Hydrophobicity scale
# Science 229:834-838(1985)
# "Rose"
hphob_rose: Final[dict[str, float]] = {_AA.A: 0.74, _AA.R: 0.64, _AA.N: 0.63, _AA.D: 0.62, _AA.C: 0.91,
      _AA.Q: 0.62, _AA.E: 0.62, _AA.G: 0.72, _AA.H: 0.78, _AA.I: 0.88,
      _AA.L: 0.85, _AA.K: 0.52, _AA.M: 0.85, _AA.F: 0.88, _AA.P: 0.64,
      _AA.S: 0.66, _AA.T: 0.7, _AA.W: 0.85, _AA.Y: 0.76, _AA.V: 0.86}

# Roseman Hydrophobicity scale
# J. Mol. Biol. 200:513-522(1988)
# "Roseman"
hphob_roseman: Final[dict[str, float]] = {_AA.A: 0.39, _AA.R: -3.95, _AA.N: -1.91, _AA.D: -3.81, _AA.C: 0.25,
      _AA.Q: -1.3, _AA.E: -2.91, _AA.G: 0, _AA.H: -0.64, _AA.I: 1.82,
      _AA.L: 1.82, _AA.K: -2.77, _AA.M: 0.96, _AA.F: 2.27, _AA.P: 0.99,
      _AA.S: -1.24, _AA.T: -1, _AA.W: 2.13, _AA.Y: 1.47, _AA.V: 1.3}

# Sweet Optimized Matchig Hydrophobicity (OMH)
# J. Mol. Biol. 171:479-488(1983).
# "Sweet
hphob_sweet: Final[dict[str, float]] = {_AA.A: -0.4, _AA.R: -0.59, _AA.N: -0.92, _AA.D: -1.31, _AA.C: 0.17,
      _AA.Q: -0.91, _AA.E: -1.22, _AA.G: -0.67, _AA.H: -0.64, _AA.I: 1.25,
      _AA.L: 1.22, _AA.K: -0.67, _AA.M: 1.02, _AA.F: 1.92, _AA.P: -0.49,
      _AA.S: -0.55, _AA.T: -0.28, _AA.W: 0.5, _AA.Y: 1.67, _AA.V: 0.91}

# Tanford Hydrophobicity scale
# J. Am. Chem. Soc. 84:4240-4274(1962)
# "Tanford"
hphob_tanford: Final[dict[str, float]] = {_AA.A: 0.62, _AA.R: -2.53, _AA.N: -0.78, _AA.D: -0.09, _AA.C: 0.29,
      _AA.Q: -0.85, _AA.E: -0.74, _AA.G: 0.48, _AA.H: -0.4, _AA.I: 1.38,
      _AA.L: 1.53, _AA.K: -1.5, _AA.M: 0.64, _AA.F: 1.19, _AA.P: 0.12,
      _AA.S: -0.18, _AA.T: -0.05, _AA.W: 0.81, _AA.Y: 0.26, _AA.V: 1.8}

# Wilson Hydrophobic constants derived from HPLC peptide retention times
# Biochem. J. 199:31-41(1981)
# "Wilson"
hphob_wilson: Final[dict[str, float]] = {_AA.A: -0.3, _AA.R: -1.1, _AA.N: -0.2, _AA.D: -1.4, _AA.C: 6.3,
      _AA.Q: -0.2, _AA.E: 0, _AA.G: 1.2, _AA.H: -1.3, _AA.I: 4.3,
      _AA.L: 6.6, _AA.K: -3.6, _AA.M: 2.5, _AA.F: 7.5, _AA.P: 2.2,
      _AA.S: -0.6, _AA.T: -2.2, _AA.W: 7.9, _AA.Y: 7.1, _AA.V: 5.9}

# Zimmerman Hydrophobicity scale
# Journal of theoretical biology, 21(2), 170-201.
# "Zimmerman"
hphob_zimmerman: Final[dict[str, float]] = {_AA.A: 0.83, _AA.R: 0.83, _AA.N: 0.09, _AA.D: 0.64, _AA.C: 1.48,
      _AA.Q: 0, _AA.E: 0.65, _AA.G: 0.1, _AA.H: 1.1, _AA.I: 3.07,
      _AA.L: 2.52, _AA.K: 1.6, _AA.M: 1.4, _AA.F: 2.75, _AA.P: 2.7,
      _AA.S: 0.14, _AA.T: 0.54, _AA.W: 0.31, _AA.Y: 2.97, _AA.V: 1.79}

"""
Amino acid scale: Proportion of residues 95% buried (in 12 proteins).
Author(s):
Chothia C.
Reference: J. Mol. Biol. 105:1-14(1976).
https://web.expasy.org/cgi-bin/protscale/protscale.pl
"""

hphob_chothia: Final[dict[str, float]] = {_AA.A: 0.38, _AA.R: 0.01, _AA.N: 0.12, _AA.D: 0.15, _AA.C: 0.5,
      _AA.Q: 0.07, _AA.E: 0.18, _AA.G: 0.36, _AA.H: 0.17, _AA.I: 0.6,
      _AA.L: 0.45, _AA.K: 0.03, _AA.M: 0.4, _AA.F: 0.5, _AA.P: 0.18,
      _AA.S: 0.22, _AA.T: 0.23, _AA.W: 0.27, _AA.Y: 0.15, _AA.V: 0.54}


""""
Amino acid scale: Free energy of transfer from inside to outside of a globular protein.
Author(s):
Janin J.
Reference: Nature 277:491-492(1979).
https://web.expasy.org/protscale/pscale/Hphob.Janin.html
"""

hphob_janin: Final[dict[str, float]] = {_AA.A: 0.3, _AA.R: -1.4, _AA.N: -0.5, _AA.D: -0.6, _AA.C: 0.9,
      _AA.Q: -0.7, _AA.E: -0.7, _AA.G: 0.3, _AA.H: -0.1, _AA.I: 0.7,
      _AA.L: 0.5, _AA.K: -1.8, _AA.M: 0.4, _AA.F: 0.5, _AA.P: -0.3,
      _AA.S: -0.1, _AA.T: -0.2, _AA.W: 0.3, _AA.Y: -0.4, _AA.V: 0.6}

"""
Amino acid scale: Hydration potential (kcal/mole) at 25�C.
Author(s):
Wolfenden R.V., Andersson L., Cullis P.M., Southgate C.C.F.
Reference: Biochemistry 20:849-855(1981).
https://web.expasy.org/protscale/pscale/Hphob.Wolfenden.html
"""

hphob_wolfenden: Final[dict[str, float]] = {_AA.A: 1.940, _AA.R: -19.920, _AA.N: -9.680, _AA.D: -10.950, _AA.C: -1.240,
      _AA.Q: -9.380, _AA.E: -10.200, _AA.G: 2.390, _AA.H: -10.270, _AA.I: 2.150,
      _AA.L: 2.280, _AA.K: -9.520, _AA.M: -1.480, _AA.F: -0.760, _AA.P: 0.000,
      _AA.S: -5.060, _AA.T: -4.880, _AA.W: -5.880, _AA.Y: -6.110, _AA.V: 1.990}



"""
Amino acid scale: Antigenicity value X 10.
Author(s):
Welling G.W., Weijer W.J., Van der Zee R., Welling-Wester S.
Reference: FEBS Lett. 188:215-218(1985).
https://web.expasy.org/protscale/pscale/Hphob.Welling.html
"""

# Welling antigenicity scale
hphob_welling: Final[dict[str, float]] = {_AA.A: 1.150, _AA.R: 0.580, _AA.N: -0.770, _AA.D: 0.650, _AA.C: -1.200,
      _AA.Q: -0.110, _AA.E: -0.710, _AA.G: -1.840, _AA.H: 3.120, _AA.I: -2.920,
      _AA.L: 0.750, _AA.K: 2.060, _AA.M: -3.850, _AA.F: -1.410, _AA.P: -0.530,
      _AA.S: -0.260, _AA.T: -0.450, _AA.W: -1.140, _AA.Y: 0.130, _AA.V: -0.130}


# Flexibility
# Normalized flexibility parameters (B-values), average
# Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).
flexibility_vihinen: Final[dict[str, float]] = {_AA.A: 0.984, _AA.C: 0.906, _AA.E: 1.094, _AA.D: 1.068,
        _AA.G: 1.031, _AA.F: 0.915, _AA.I: 0.927, _AA.H: 0.950,
        _AA.K: 1.102, _AA.M: 0.952, _AA.L: 0.935, _AA.N: 1.048,
        _AA.Q: 1.037, _AA.P: 1.049, _AA.S: 1.046, _AA.R: 1.008,
        _AA.T: 0.997, _AA.W: 0.904, _AA.V: 0.931, _AA.Y: 0.929}

# Hydrophilicity
# 1 Hopp & Wood
# Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).
hydrophilicity_hopp_wood: Final[dict[str, float]] = {_AA.A: -0.5, _AA.R: 3.0, _AA.N: 0.2, _AA.D: 3.0, _AA.C: -1.0,
      _AA.Q: 0.2, _AA.E: 3.0, _AA.G: 0.0, _AA.H: -0.5, _AA.I: -1.8,
      _AA.L: -1.8, _AA.K: 3.0, _AA.M: -1.3, _AA.F: -2.5, _AA.P: 0.0,
      _AA.S: 0.3, _AA.T: -0.4, _AA.W: -3.4, _AA.Y: -2.3, _AA.V: -1.5}

# Surface accessibility
# Vergoten G & Theophanides T, Biomolecular Structure and Dynamics,
# pg.138 (1997).
# 1 Emini Surface fractional probability
surface_accessibility_vergoten: Final[dict[str, float]] = {_AA.A: 0.815, _AA.R: 1.475, _AA.N: 1.296, _AA.D: 1.283, _AA.C: 0.394,
      _AA.Q: 1.348, _AA.E: 1.445, _AA.G: 0.714, _AA.H: 1.180, _AA.I: 0.603,
      _AA.L: 0.603, _AA.K: 1.545, _AA.M: 0.714, _AA.F: 0.695, _AA.P: 1.236,
      _AA.S: 1.115, _AA.T: 1.184, _AA.W: 0.808, _AA.Y: 1.089, _AA.V: 0.606}

# 2 Janin Interior to surface transfer energy scale
surface_accessiblility_janin: Final[dict[str, float]] = {_AA.A: 0.28, _AA.R: -1.14, _AA.N: -0.55, _AA.D: -0.52, _AA.C: 0.97,
      _AA.Q: -0.69, _AA.E: -1.01, _AA.G: 0.43, _AA.H: -0.31, _AA.I: 0.60,
      _AA.L: 0.60, _AA.K: -1.62, _AA.M: 0.43, _AA.F: 0.46, _AA.P: -0.42,
      _AA.S: -0.19, _AA.T: -0.32, _AA.W: 0.29, _AA.Y: -0.15, _AA.V: 0.60}


# A two dimensional dictionary for calculating the instability index.
# Guruprasad K., Reddy B.V.B., Pandit M.W. Protein Engineering 4:155-161(1990).
# It is based on dipeptide values; therefore, the value for the dipeptide DG
# is DIWV['D']['G'].
DIWV: Final[dict[_AA, dict[str, float]]] = {
        _AA.A: {_AA.A: 1.0, _AA.C: 44.94, _AA.E: 1.0, _AA.D: -7.49,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 1.0, _AA.H: -7.49,
              _AA.K: 1.0, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 1.0, _AA.P: 20.26, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: 1.0, _AA.W: 1.0, _AA.V: 1.0, _AA.Y: 1.0},
        _AA.C: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 20.26,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 33.60,
              _AA.K: 1.0, _AA.M: 33.60, _AA.L: 20.26, _AA.N: 1.0,
              _AA.Q: -6.54, _AA.P: 20.26, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: 33.60, _AA.W: 24.68, _AA.V: -6.54, _AA.Y: 1.0},
        _AA.E: {_AA.A: 1.0, _AA.C: 44.94, _AA.E: 33.60, _AA.D: 20.26,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 20.26, _AA.H: -6.54,
              _AA.K: 1.0, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 20.26, _AA.P: 20.26, _AA.S: 20.26, _AA.R: 1.0,
              _AA.T: 1.0, _AA.W: -14.03, _AA.V: 1.0, _AA.Y: 1.0},
        _AA.D: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: 1.0, _AA.F: -6.54, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: -7.49, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 1.0, _AA.P: 1.0, _AA.S: 20.26, _AA.R: -6.54,
              _AA.T: -14.03, _AA.W: 1.0, _AA.V: 1.0, _AA.Y: 1.0},
        _AA.G: {_AA.A: -7.49, _AA.C: 1.0, _AA.E: -6.54, _AA.D: 1.0,
              _AA.G: 13.34, _AA.F: 1.0, _AA.I: -7.49, _AA.H: 1.0,
              _AA.K: -7.49, _AA.M: 1.0, _AA.L: 1.0, _AA.N: -7.49,
              _AA.Q: 1.0, _AA.P: 1.0, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: -7.49, _AA.W: 13.34, _AA.V: 1.0, _AA.Y: -7.49},
        _AA.F: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 13.34,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: -14.03, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 1.0, _AA.P: 20.26, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: 1.0, _AA.W: 1.0, _AA.V: 1.0, _AA.Y: 33.601},
        _AA.I: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 44.94, _AA.D: 1.0,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 13.34,
              _AA.K: -7.49, _AA.M: 1.0, _AA.L: 20.26, _AA.N: 1.0,
              _AA.Q: 1.0, _AA.P: -1.88, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: 1.0, _AA.W: 1.0, _AA.V: -7.49, _AA.Y: 1.0},
        _AA.H: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: -9.37, _AA.F: -9.37, _AA.I: 44.94, _AA.H: 1.0,
              _AA.K: 24.68, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 24.68,
              _AA.Q: 1.0, _AA.P: -1.88, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: -6.54, _AA.W: -1.88, _AA.V: 1.0, _AA.Y: 44.94},
        _AA.K: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: -7.49, _AA.F: 1.0, _AA.I: -7.49, _AA.H: 1.0,
              _AA.K: 1.0, _AA.M: 33.60, _AA.L: -7.49, _AA.N: 1.0,
              _AA.Q: 24.64, _AA.P: -6.54, _AA.S: 1.0, _AA.R: 33.60,
              _AA.T: 1.0, _AA.W: 1.0, _AA.V: -7.49, _AA.Y: 1.0},
        _AA.M: {_AA.A: 13.34, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 58.28,
              _AA.K: 1.0, _AA.M: -1.88, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: -6.54, _AA.P: 44.94, _AA.S: 44.94, _AA.R: -6.54,
              _AA.T: -1.88, _AA.W: 1.0, _AA.V: 1.0, _AA.Y: 24.68},
        _AA.L: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: -7.49, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 33.60, _AA.P: 20.26, _AA.S: 1.0, _AA.R: 20.26,
              _AA.T: 1.0, _AA.W: 24.68, _AA.V: 1.0, _AA.Y: 1.0},
        _AA.N: {_AA.A: 1.0, _AA.C: -1.88, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: -14.03, _AA.F: -14.03, _AA.I: 44.94, _AA.H: 1.0,
              _AA.K: 24.68, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: -6.54, _AA.P: -1.88, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: -7.49, _AA.W: -9.37, _AA.V: 1.0, _AA.Y: 1.0},
        _AA.Q: {_AA.A: 1.0, _AA.C: -6.54, _AA.E: 20.26, _AA.D: 20.26,
              _AA.G: 1.0, _AA.F: -6.54, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: 1.0, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 20.26, _AA.P: 20.26, _AA.S: 44.94, _AA.R: 1.0,
              _AA.T: 1.0, _AA.W: 1.0, _AA.V: -6.54, _AA.Y: -6.54},
        _AA.P: {_AA.A: 20.26, _AA.C: -6.54, _AA.E: 18.38, _AA.D: -6.54,
              _AA.G: 1.0, _AA.F: 20.26, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: 1.0, _AA.M: -6.54, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 20.26, _AA.P: 20.26, _AA.S: 20.26, _AA.R: -6.54,
              _AA.T: 1.0, _AA.W: -1.88, _AA.V: 20.26, _AA.Y: 1.0},
        _AA.S: {_AA.A: 1.0, _AA.C: 33.60, _AA.E: 20.26, _AA.D: 1.0,
              _AA.G: 1.0, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: 1.0, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 20.26, _AA.P: 44.94, _AA.S: 20.26, _AA.R: 20.26,
              _AA.T: 1.0, _AA.W: 1.0, _AA.V: 1.0, _AA.Y: 1.0},
        _AA.R: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: -7.49, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 20.26,
              _AA.K: 1.0, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 13.34,
              _AA.Q: 20.26, _AA.P: 20.26, _AA.S: 44.94, _AA.R: 58.28,
              _AA.T: 1.0, _AA.W: 58.28, _AA.V: 1.0, _AA.Y: -6.54},
        _AA.T: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 20.26, _AA.D: 1.0,
              _AA.G: -7.49, _AA.F: 13.34, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: 1.0, _AA.M: 1.0, _AA.L: 1.0, _AA.N: -14.03,
              _AA.Q: -6.54, _AA.P: 1.0, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: 1.0, _AA.W: -14.03, _AA.V: 1.0, _AA.Y: 1.0},
        _AA.W: {_AA.A: -14.03, _AA.C: 1.0, _AA.E: 1.0, _AA.D: 1.0,
              _AA.G: -9.37, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 24.68,
              _AA.K: 1.0, _AA.M: 24.68, _AA.L: 13.34, _AA.N: 13.34,
              _AA.Q: 1.0, _AA.P: 1.0, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: -14.03, _AA.W: 1.0, _AA.V: -7.49, _AA.Y: 1.0},
        _AA.V: {_AA.A: 1.0, _AA.C: 1.0, _AA.E: 1.0, _AA.D: -14.03,
              _AA.G: -7.49, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 1.0,
              _AA.K: -1.88, _AA.M: 1.0, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 1.0, _AA.P: 20.26, _AA.S: 1.0, _AA.R: 1.0,
              _AA.T: -7.49, _AA.W: 1.0, _AA.V: 1.0, _AA.Y: -6.54},
        _AA.Y: {_AA.A: 24.68, _AA.C: 1.0, _AA.E: -6.54, _AA.D: 24.68,
              _AA.G: -7.49, _AA.F: 1.0, _AA.I: 1.0, _AA.H: 13.34,
              _AA.K: 1.0, _AA.M: 44.94, _AA.L: 1.0, _AA.N: 1.0,
              _AA.Q: 1.0, _AA.P: 13.34, _AA.S: 1.0, _AA.R: -15.91,
              _AA.T: -7.49, _AA.W: -9.37, _AA.V: 1.0, _AA.Y: 13.34},
        }



"""
https://www.peptideweb.com/images/pdf/pKa-and-pI-values-of-amino-acids.pdf
"""

pk_nterminal: Final[dict[str, float]] = {
      _AA.A: 9.69, _AA.R: 9.04, _AA.N: 8.80, _AA.D: 9.82, _AA.C: 10.78,
      _AA.Q: 9.13, _AA.E: 9.67, _AA.G: 9.60, _AA.H: 9.17, _AA.I: 9.60,
      _AA.L: 9.60, _AA.K: 8.95, _AA.M: 9.21, _AA.F: 9.13, _AA.P: 10.60,
      _AA.S: 9.15, _AA.T: 9.10, _AA.W: 9.44, _AA.Y: 9.11, _AA.V: 9.62
}
pk_cterminal: Final[dict[str, float]] = {
      _AA.A: 2.34, _AA.R: 2.17, _AA.N: 2.02, _AA.D: 2.09, _AA.C: 1.71,
      _AA.Q: 2.17, _AA.E: 2.19, _AA.G: 2.34, _AA.H: 1.82, _AA.I: 2.36,
      _AA.L: 2.36, _AA.K: 2.18, _AA.M: 2.28, _AA.F: 1.83, _AA.P: 1.99,
      _AA.S: 2.21, _AA.T: 2.09, _AA.W: 2.43, _AA.Y: 2.20, _AA.V: 2.32
}

pk_sidechain: Final[dict[str, float]] = {
      _AA.A: 0.0, _AA.R: 12.48, _AA.N: 0.0, _AA.D: 3.86, _AA.C: 8.33,
      _AA.Q: 0.0, _AA.E: 4.25, _AA.G: 0.0, _AA.H: 6.00, _AA.I: 0.0,
      _AA.L: 0.0, _AA.K: 10.79, _AA.M: 0.0, _AA.F: 0.0, _AA.P: 0.0,
      _AA.S: 0.0, _AA.T: 0.0, _AA.W: 0.0, _AA.Y: 10.07, _AA.V: 0.0
}


"""
Amino acid scale: Conformational parameter for alpha helix.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/alpha-helixRoux.html
"""

deleage_roux_alpha_helix: Final[dict[str, float]] = {
    _AA.A: 1.489, _AA.R: 1.224, _AA.N: 0.772, _AA.D: 0.924, _AA.C: 0.966,
    _AA.Q: 1.164, _AA.E: 1.504, _AA.G: 0.510, _AA.H: 1.003, _AA.I: 1.003,
    _AA.L: 1.236, _AA.K: 1.172, _AA.M: 1.363, _AA.F: 1.195, _AA.P: 0.492,
    _AA.S: 0.739, _AA.T: 0.785, _AA.W: 1.090, _AA.Y: 0.787, _AA.V: 0.990
}

"""
Amino acid scale: Conformational parameter for beta-sheet.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/beta-sheetRoux.html
"""

deleage_roux_beta_sheet: Final[dict[str, float]] = {
    _AA.A: 0.709, _AA.R: 0.920, _AA.N: 0.604, _AA.D: 0.541, _AA.C: 1.191,
    _AA.Q: 0.840, _AA.E: 0.567, _AA.G: 0.657, _AA.H: 0.863, _AA.I: 1.799,
    _AA.L: 1.261, _AA.K: 0.721, _AA.M: 1.210, _AA.F: 1.393, _AA.P: 0.354,
    _AA.S: 0.928, _AA.T: 1.221, _AA.W: 1.306, _AA.Y: 1.266, _AA.V: 1.965
}

"""
Amino acid scale: Conformational parameter for beta-turn.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/beta-turnRoux.html
"""

deleage_roux_beta_turn: Final[dict[str, float]] = {
    _AA.A: 0.788, _AA.R: 0.912, _AA.N: 1.572, _AA.D: 1.197, _AA.C: 0.965,
    _AA.Q: 0.997, _AA.E: 1.149, _AA.G: 1.860, _AA.H: 0.970, _AA.I: 0.240,
    _AA.L: 0.670, _AA.K: 1.302, _AA.M: 0.436, _AA.F: 0.624, _AA.P: 1.415,
    _AA.S: 1.316, _AA.T: 0.739, _AA.W: 0.546, _AA.Y: 0.795, _AA.V: 0.387
}

"""
Amino acid scale: Conformational parameter for coil.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/CoilRoux.html
"""

deleage_roux_coil: Final[dict[str, float]] = {
    _AA.A: 0.824, _AA.R: 0.893, _AA.N: 1.167, _AA.D: 1.197, _AA.C: 0.953,
    _AA.Q: 0.947, _AA.E: 0.761, _AA.G: 1.251, _AA.H: 1.068, _AA.I: 0.886,
    _AA.L: 0.810, _AA.K: 0.897, _AA.M: 0.810, _AA.F: 0.797, _AA.P: 1.540,
    _AA.S: 1.130, _AA.T: 1.148, _AA.W: 0.941, _AA.Y: 1.109, _AA.V: 0.772
}


"""
Amino acid scale: Normalized frequency for alpha helix.
Author(s):
Levitt M.
Reference: Biochemistry 17:4277-4285(1978).
https://web.expasy.org/protscale/pscale/alpha-helixLevitt.html
"""

levitt_alpha_helix: Final[dict[str, float]] = {
    _AA.A: 1.290, _AA.R: 0.960, _AA.N: 0.900, _AA.D: 1.040, _AA.C: 1.110,
    _AA.Q: 1.270, _AA.E: 1.440, _AA.G: 0.560, _AA.H: 1.220, _AA.I: 0.970,
    _AA.L: 1.300, _AA.K: 1.230, _AA.M: 1.470, _AA.F: 1.070, _AA.P: 0.520,
    _AA.S: 0.820, _AA.T: 0.820, _AA.W: 0.990, _AA.Y: 0.720, _AA.V: 0.910
}

"""
Amino acid scale: Normalized frequency for beta-sheet.
Author(s):
Levitt M.
Reference: Biochemistry 17:4277-4285(1978).
https://web.expasy.org/protscale/pscale/beta-sheetLevitt.html
"""

levitt_beta_sheet: Final[dict[str, float]] = {
    _AA.A: 0.900, _AA.R: 0.990, _AA.N: 0.760, _AA.D: 0.720, _AA.C: 0.740,
    _AA.Q: 0.800, _AA.E: 0.750, _AA.G: 0.920, _AA.H: 1.080, _AA.I: 1.450,
    _AA.L: 1.020, _AA.K: 0.770, _AA.M: 0.970, _AA.F: 1.320, _AA.P: 0.640,
    _AA.S: 0.950, _AA.T: 1.210, _AA.W: 1.140, _AA.Y: 1.250, _AA.V: 1.490
}

"""
Amino acid scale: Normalized frequency for beta-turn.
Author(s):
Levitt M.
Reference: Biochemistry 17:4277-4285(1978).
https://web.expasy.org/protscale/pscale/beta-turnLevitt.html
"""

levitt_beta_turn: Final[dict[str, float]] = {
    _AA.A: 0.770, _AA.R: 0.880, _AA.N: 1.280, _AA.D: 1.410, _AA.C: 0.810,
    _AA.Q: 0.980, _AA.E: 0.990, _AA.G: 1.640, _AA.H: 0.680, _AA.I: 0.510,
    _AA.L: 0.580, _AA.K: 0.960, _AA.M: 0.410, _AA.F: 0.590, _AA.P: 1.910,
    _AA.S: 1.320, _AA.T: 1.040, _AA.W: 0.760, _AA.Y: 1.050, _AA.V: 0.470
}

"""
Amino acid scale: Conformational parameter for alpha helix (computed from 29 proteins).
Author(s):
Chou P.Y., Fasman G.D.
Reference: Adv. Enzym. 47:45-148(1978).
https://web.expasy.org/protscale/pscale/alpha-helixFasman.html
"""

chou_fasman_alpha_helix: Final[dict[str, float]] = {
    _AA.A: 1.420, _AA.R: 0.980, _AA.N: 0.670, _AA.D: 1.010, _AA.C: 0.700,
    _AA.Q: 1.110, _AA.E: 1.510, _AA.G: 0.570, _AA.H: 1.000, _AA.I: 1.080,
    _AA.L: 1.210, _AA.K: 1.160, _AA.M: 1.450, _AA.F: 1.130, _AA.P: 0.570,
    _AA.S: 0.770, _AA.T: 0.830, _AA.W: 1.080, _AA.Y: 0.690, _AA.V: 1.060
}

"""
Amino acid scale: Conformational parameter for beta-sheet (computed from 29 proteins).
Author(s):
Chou P.Y., Fasman G.D.
Reference: Adv. Enzym. 47:45-148(1978).
  https://web.expasy.org/protscale/pscale/beta-sheetFasman.html
"""

chou_fasman_beta_sheet: Final[dict[str, float]] = {
    _AA.A: 0.830, _AA.R: 0.930, _AA.N: 0.890, _AA.D: 0.540, _AA.C: 1.190,
    _AA.Q: 1.100, _AA.E: 0.370, _AA.G: 0.750, _AA.H: 0.870, _AA.I: 1.600,
    _AA.L: 1.300, _AA.K: 0.740, _AA.M: 1.050, _AA.F: 1.380, _AA.P: 0.550,
    _AA.S: 0.750, _AA.T: 1.190, _AA.W: 1.370, _AA.Y: 1.470, _AA.V: 1.700
}

"""
Amino acid scale: Conformational parameter for beta-turn (computed from 29 proteins).
Author(s):
Chou P.Y., Fasman G.D.
Reference: Adv. Enzym. 47:45-148(1978).
https://web.expasy.org/protscale/pscale/beta-turnFasman.html
"""

chou_fasman_beta_turn: Final[dict[str, float]] = {
    _AA.A: 0.660, _AA.R: 0.950, _AA.N: 1.560, _AA.D: 1.460, _AA.C: 1.190,
    _AA.Q: 0.980, _AA.E: 0.740, _AA.G: 1.560, _AA.H: 0.950, _AA.I: 0.470,
    _AA.L: 0.590, _AA.K: 1.010, _AA.M: 0.600, _AA.F: 0.600, _AA.P: 1.520,
    _AA.S: 1.430, _AA.T: 0.960, _AA.W: 0.960, _AA.Y: 1.140, _AA.V: 0.500
}





"""
Amino acid scale: Conformational preference for parallel beta strand.
Author(s):
Lifson S., Sander C.
Reference: Nature 282:109-111(1979).
https://web.expasy.org/protscale/pscale/Parallelbeta-strand.html
"""

beta_strand_parallel: Final[dict[str, float]] = {
    _AA.A: 1.000, _AA.R: 0.680, _AA.N: 0.540, _AA.D: 0.500, _AA.C: 0.910,
    _AA.Q: 0.280, _AA.E: 0.590, _AA.G: 0.790, _AA.H: 0.380, _AA.I: 2.600,
    _AA.L: 1.420, _AA.K: 0.590, _AA.M: 1.490, _AA.F: 1.300, _AA.P: 0.350,
    _AA.S: 0.700, _AA.T: 0.590, _AA.W: 0.890, _AA.Y: 1.080, _AA.V: 2.630
}

"""
Amino acid scale: Conformational preference for antiparallel beta strand.
Author(s):
Lifson S., Sander C.
Reference: Nature 282:109-111(1979).
https://web.expasy.org/protscale/pscale/Antiparallelbeta-strand.html
"""

beta_strand_antiparallel: Final[dict[str, float]] = {
    _AA.A: 0.900, _AA.R: 1.020, _AA.N: 0.620, _AA.D: 0.470, _AA.C: 1.240,
    _AA.Q: 1.180, _AA.E: 0.620, _AA.G: 0.560, _AA.H: 1.120, _AA.I: 1.540,
    _AA.L: 1.260, _AA.K: 0.740, _AA.M: 1.090, _AA.F: 1.230, _AA.P: 0.420,
    _AA.S: 0.870, _AA.T: 1.300, _AA.W: 1.750, _AA.Y: 1.680, _AA.V: 1.530
}

"""
Amino acid scale: Conformational preference for total beta strand (antiparallel+parallel).
Author(s):
Lifson S., Sander C.
Reference: Nature 282:109-111(1979).
https://web.expasy.org/protscale/pscale/Totalbeta-strand.html
"""

beta_strand_total: Final[dict[str, float]] = {
    _AA.A: 0.920, _AA.R: 0.930, _AA.N: 0.600, _AA.D: 0.480, _AA.C: 1.160,
    _AA.Q: 0.950, _AA.E: 0.610, _AA.G: 0.610, _AA.H: 0.930, _AA.I: 1.810,
    _AA.L: 1.300, _AA.K: 0.700, _AA.M: 1.190, _AA.F: 1.250, _AA.P: 0.400,
    _AA.S: 0.820, _AA.T: 1.120, _AA.W: 1.540, _AA.Y: 1.530, _AA.V: 1.810
}

"""
Amino acid scale: Atomic weight ratio of hetero elements in end group to C in side chain.
Author(s):
Grantham R.
Reference: Science 185:862-864(1974).
https://web.expasy.org/protscale/pscale/Ratioside.html
"""

ratioside: Final[dict[str, float]] = {
    _AA.A: 0.000, _AA.R: 0.650, _AA.N: 1.330, _AA.D: 1.380, _AA.C: 2.750,
    _AA.Q: 0.890, _AA.E: 0.920, _AA.G: 0.740, _AA.H: 0.580, _AA.I: 0.000,
    _AA.L: 0.000, _AA.K: 0.330, _AA.M: 0.000, _AA.F: 0.000, _AA.P: 0.390,
    _AA.S: 1.420, _AA.T: 0.710, _AA.W: 0.130, _AA.Y: 0.200, _AA.V: 0.000
}


"""
Amino acid scale: Polarity (p).
Author(s):
Grantham R.
Reference: Science 185:862-864(1974).
https://web.expasy.org/protscale/pscale/PolarityGrantham.html
"""

polarity_grantham: Final[dict[str, float]] = {
      _AA.A: 8.100, _AA.R: 10.500, _AA.N: 11.600, _AA.D: 13.000, _AA.C: 5.500,
      _AA.Q: 10.500, _AA.E: 12.300, _AA.G: 9.000, _AA.H: 10.400, _AA.I: 5.200,
      _AA.L: 4.900, _AA.K: 11.300, _AA.M: 5.700, _AA.F: 5.200, _AA.P: 8.000,
      _AA.S: 9.200, _AA.T: 8.600, _AA.W: 5.400, _AA.Y: 6.200, _AA.V: 5.900
      }

"""
Amino acid scale: Polarity.
Author(s):
Zimmerman J.M., Eliezer N., Simha R.
Reference: J. Theor. Biol. 21:170-201(1968).
https://web.expasy.org/protscale/pscale/PolarityZimmerman.html
"""

polarity_zimmerman: Final[dict[str, float]] = {
    _AA.A: 0.000, _AA.R: 52.000, _AA.N: 3.380, _AA.D: 49.700, _AA.C: 1.480,
    _AA.Q: 3.530, _AA.E: 49.900, _AA.G: 0.000, _AA.H: 51.600, _AA.I: 0.130,
    _AA.L: 0.130, _AA.K: 49.500, _AA.M: 1.430, _AA.F: 0.350, _AA.P: 1.580,
    _AA.S: 1.670, _AA.T: 1.660, _AA.W: 2.100, _AA.Y: 1.610, _AA.V: 0.130
}

"""
Amino acid scale: Relative mutability of amino acids (Ala=100).
Author(s):
Dayhoff M.O., Schwartz R.M., Orcutt B.C.
Reference: In "Atlas of Protein Sequence and Structure", Vol.5, Suppl.3 (1978).
https://web.expasy.org/protscale/pscale/Relativemutability.html
"""

mutability: Final[dict[str, float]] = {
      _AA.A: 100.000, _AA.R: 65.000, _AA.N: 134.000, _AA.D: 106.000, _AA.C: 20.000,
      _AA.Q: 93.000, _AA.E: 102.000, _AA.G: 49.000, _AA.H: 66.000, _AA.I: 96.000,
      _AA.L: 40.000, _AA.K: 56.000, _AA.M: 94.000, _AA.F: 41.000, _AA.P: 56.000,
      _AA.S: 120.000, _AA.T: 97.000, _AA.W: 18.000, _AA.Y: 41.000, _AA.V: 74.000
}

"""
Amino acid scale: Number of codon(s) coding for each amino acid in universal genetic code.
Author(s):
-
Reference: Most textbooks.
https://web.expasy.org/protscale/pscale/Numbercodons.html
"""

codons: Final[dict[str, float]] = {
    _AA.A: 4.0, _AA.R: 6.0, _AA.N: 2.0, _AA.D: 2.0, _AA.C: 1.0,
    _AA.Q: 2.0, _AA.E: 2.0, _AA.G: 4.0, _AA.H: 2.0, _AA.I: 3.0,
    _AA.L: 6.0, _AA.K: 2.0, _AA.M: 1.0, _AA.F: 2.0, _AA.P: 4.0,
    _AA.S: 6.0, _AA.T: 4.0, _AA.W: 1.0, _AA.Y: 2.0, _AA.V: 4.0
}

"""
Amino acid scale: Refractivity.
Author(s):
Jones. D.D.
Reference: J. Theor. Biol. 50:167-184(1975).
https://web.expasy.org/protscale/pscale/Refractivity.html
"""

refractivity: Final[dict[str, float]] = {
    _AA.A: 4.340, _AA.R: 26.660, _AA.N: 13.280, _AA.D: 12.000, _AA.C: 35.770,
    _AA.Q: 17.560, _AA.E: 17.260, _AA.G: 0.000, _AA.H: 21.810, _AA.I: 19.060,
    _AA.L: 18.780, _AA.K: 21.290, _AA.M: 21.640, _AA.F: 29.400, _AA.P: 10.930,
    _AA.S: 6.350, _AA.T: 11.010, _AA.W: 42.530, _AA.Y: 31.530, _AA.V: 13.920
}

"""
Amino acid scale: Bulkiness.
Author(s):
Zimmerman J.M., Eliezer N., Simha R.
Reference: J. Theor. Biol. 21:170-201(1968).

https://web.expasy.org/protscale/pscale/Bulkiness.html
"""

bulkiness: Final[dict[str, float]] = {
    _AA.A: 11.500, _AA.R: 14.280, _AA.N: 12.820, _AA.D: 11.680, _AA.C: 13.460,
    _AA.Q: 14.450, _AA.E: 13.570, _AA.G: 3.400, _AA.H: 13.690, _AA.I: 21.400,
    _AA.L: 21.400, _AA.K: 15.710, _AA.M: 16.250, _AA.F: 19.800, _AA.P: 17.430,
    _AA.S: 9.470, _AA.T: 15.770, _AA.W: 21.670, _AA.Y: 18.030, _AA.V: 21.570
}

"""
Amino acid scale: Recognition factors.
Author(s):
Fraga S.
Reference: Can. J. Chem. 60:2606-2610(1982).
https://web.expasy.org/protscale/pscale/Recognitionfactors.html
"""

recognition_factors: Final[dict[str, float]] = {
    _AA.A: 78.000, _AA.R: 95.000, _AA.N: 94.000, _AA.D: 81.000, _AA.C: 89.000,
    _AA.Q: 87.000, _AA.E: 78.000, _AA.G: 84.000, _AA.H: 84.000, _AA.I: 88.000,
    _AA.L: 85.000, _AA.K: 87.000, _AA.M: 80.000, _AA.F: 81.000, _AA.P: 91.000,
    _AA.S: 107.000, _AA.T: 93.000, _AA.W: 104.000, _AA.Y: 84.000, _AA.V: 89.000
}

"""
Amino acid scale: Overall amino acid composition (%).
Author(s):
McCaldon P., Argos P.
Reference: Proteins: Structure, Function and Genetics 4:99-122(1988).
https://web.expasy.org/protscale/pscale/A.A.composition.html
"""

aa_composition_mccaldron: Final[dict[str, float]] = {
      _AA.A: 8.300, _AA.R: 5.700, _AA.N: 4.400, _AA.D: 5.300, _AA.C: 1.700,
      _AA.Q: 4.000, _AA.E: 6.200, _AA.G: 7.200, _AA.H: 2.200, _AA.I: 5.200,
      _AA.L: 9.000, _AA.K: 5.700, _AA.M: 2.400, _AA.F: 3.900, _AA.P: 5.100,
      _AA.S: 6.900, _AA.T: 5.800, _AA.W: 1.300, _AA.Y: 3.200, _AA.V: 6.600
      }

"""
Amino acid scale: Amino acid composition (%) in the UniProtKB/Swiss-Prot data bank.
Author(s):
Bairoch A.
Reference: Release notes for UniProtKB/Swiss-Prot release 2013_04 - April 2013.
https://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html
"""

aa_composition_swissprot: Final[dict[str, float]] = {
    _AA.A: 8.25, _AA.R: 5.53, _AA.N: 4.06, _AA.D: 5.45, _AA.C: 1.37,
    _AA.Q: 3.93, _AA.E: 6.75, _AA.G: 7.07, _AA.H: 2.27, _AA.I: 5.96,
    _AA.L: 9.66, _AA.K: 5.84, _AA.M: 2.42, _AA.F: 3.86, _AA.P: 4.70,
    _AA.S: 6.56, _AA.T: 5.34, _AA.W: 1.08, _AA.Y: 2.92, _AA.V: 6.87
}

aa_composition_scales: Final[dict[str, dict[str, float]]] = {
    "McCaldron": aa_composition_mccaldron,
    "Swiss-Prot": aa_composition_swissprot
}

"""
Amino acid scale: Transmembrane tendency
Author(s):
Zhao, G., London E.
Reference: Protein Sci. 15:1987-2001(2006).
https://web.expasy.org/protscale/pscale/Transmembranetendency.html
"""

transmembrane_tendency: Final[dict[str, float]] = {
    _AA.A: 0.380, _AA.R: -2.570, _AA.N: -1.620, _AA.D: -3.270, _AA.C: -0.300,
    _AA.Q: -1.840, _AA.E: -2.900, _AA.G: -0.190, _AA.H: -1.440, _AA.I: 1.970,
    _AA.L: 1.820, _AA.K: -3.460, _AA.M: 1.400, _AA.F: 1.980, _AA.P: -1.440,
    _AA.S: -0.530, _AA.T: -0.320, _AA.W: 1.530, _AA.Y: 0.490, _AA.V: 1.460
}

"""
Amino acid scale: Molar fraction (%) of 3220 accessible residues.
Author(s):
Janin J.
Reference: Nature 277:491-492(1979).
https://web.expasy.org/protscale/pscale/accessibleresidues.html
"""

accessible_residues: Final[dict[str, float]] = {
    _AA.A: 6.600, _AA.R: 4.500, _AA.N: 6.700, _AA.D: 7.700, _AA.C: 0.900,
    _AA.Q: 5.200, _AA.E: 5.700, _AA.G: 6.700, _AA.H: 2.500, _AA.I: 2.800,
    _AA.L: 4.800, _AA.K: 10.300, _AA.M: 1.000, _AA.F: 2.400, _AA.P: 4.800,
    _AA.S: 9.400, _AA.T: 7.000, _AA.W: 1.400, _AA.Y: 5.100, _AA.V: 4.500
}

"""
Amino acid scale: Average area buried on transfer from standard state to folded protein.
Author(s):
Rose G.D., Geselowitz A.R., Lesser G.J., Lee R.H., Zehfus M.H.
Reference: Science 229:834-838(1985).
https://web.expasy.org/protscale/pscale/Averageburied.html
"""

average_buried_area: Final[dict[str, float]] = {
    _AA.A: 86.600, _AA.R: 162.200, _AA.N: 103.300, _AA.D: 97.800, _AA.C: 132.300,
    _AA.Q: 119.200, _AA.E: 113.900, _AA.G: 62.900, _AA.H: 155.800, _AA.I: 158.000,
    _AA.L: 164.100, _AA.K: 115.500, _AA.M: 172.900, _AA.F: 194.100, _AA.P: 92.900,
    _AA.S: 85.600, _AA.T: 106.500, _AA.W: 224.600, _AA.Y: 177.700, _AA.V: 141.000
}

"""
Amino acid scale: Molecular weight of each amino acid.
Author(s):
-
Reference: Most textbooks.
https://web.expasy.org/protscale/pscale/Molecularweight.html
"""

molecular_weights: Final[dict[str, float]] = {
    _AA.A: 89.000, _AA.R: 174.000, _AA.N: 132.000, _AA.D: 133.000, _AA.C: 121.000,
    _AA.Q: 146.000, _AA.E: 147.000, _AA.G: 75.000, _AA.H: 155.000, _AA.I: 131.000,
    _AA.L: 131.000, _AA.K: 146.000, _AA.M: 149.000, _AA.F: 165.000, _AA.P: 115.000,
    _AA.S: 105.000, _AA.T: 119.000, _AA.W: 204.000, _AA.Y: 181.000, _AA.V: 117.000
}

"""
Amino acid scale: Retention coefficient in HPLC, pH 2.1.
Author(s):
Meek J.L.
Reference: Proc. Natl. Acad. Sci. USA 77:1632-1636(1980).
https://web.expasy.org/protscale/pscale/HPLC2.1.html
"""

hplc_meek_2_1: Final[dict[str, float]] = {
    _AA.A: -0.100, _AA.R: -4.500, _AA.N: -1.600, _AA.D: -2.800, _AA.C: -2.200,
    _AA.Q: -2.500, _AA.E: -7.500, _AA.G: -0.500, _AA.H: 0.800, _AA.I: 11.800,
    _AA.L: 10.000, _AA.K: -3.200, _AA.M: 7.100, _AA.F: 13.900, _AA.P: 8.000,
    _AA.S: -3.700, _AA.T: 1.500, _AA.W: 18.100, _AA.Y: 8.200, _AA.V: 3.300
}

"""
Amino acid scale: Retention coefficient in HFBA.
Author(s):
Browne C.A., Bennett H.P.J., Solomon S.
Reference: Anal. Biochem. 124:201-208(1982).
https://web.expasy.org/protscale/pscale/HPLCHFBA.html

"""

hplc_browne: Final[dict[str, float]] = {
    _AA.A: 3.900, _AA.R: 3.200, _AA.N: -2.800, _AA.D: -2.800, _AA.C: -14.300,
    _AA.Q: 1.800, _AA.E: -7.500, _AA.G: -2.300, _AA.H: 2.000, _AA.I: 11.000,
    _AA.L: 15.000, _AA.K: -2.500, _AA.M: 4.100, _AA.F: 14.700, _AA.P: 5.600,
    _AA.S: -3.500, _AA.T: 1.100, _AA.W: 17.800, _AA.Y: 3.800, _AA.V: 2.100
}

"""
Amino acid scale: Retention coefficient in HPLC, pH 7.4.
Author(s):
Meek J.L.
Reference: Proc. Natl. Acad. Sci. USA 77:1632-1636(1980).
https://web.expasy.org/protscale/pscale/HPLC7.4.html
"""

hplc_meek_7_4: Final[dict[str, float]] = {
    _AA.A: 0.500, _AA.R: 0.800, _AA.N: 0.800, _AA.D: -8.200, _AA.C: -6.800,
    _AA.Q: -4.800, _AA.E: -16.900, _AA.G: 0.000, _AA.H: -3.500, _AA.I: 13.900,
    _AA.L: 8.800, _AA.K: 0.100, _AA.M: 4.800, _AA.F: 13.200, _AA.P: 6.100,
    _AA.S: 1.200, _AA.T: 2.700, _AA.W: 14.900, _AA.Y: 6.100, _AA.V: 2.700
}

"""
Amino acid scale: Retention coefficient in TFA.
Author(s):
Browne C.A., Bennett H.P.J., Solomon S.
Reference: Anal. Biochem. 124:201-208(1982).
https://web.expasy.org/protscale/pscale/HPLCTFA.html
"""

hplc_browne_tfa: Final[dict[str, float]] = {
    _AA.A: 7.300, _AA.R: -3.600, _AA.N: -5.700, _AA.D: -2.900, _AA.C: -9.200,
    _AA.Q: -0.300, _AA.E: -7.100, _AA.G: -1.200, _AA.H: -2.100, _AA.I: 6.600,
    _AA.L: 20.000, _AA.K: -3.700, _AA.M: 5.600, _AA.F: 19.200, _AA.P: 5.100,
    _AA.S: -4.100, _AA.T: 0.800, _AA.W: 16.300, _AA.Y: 5.900, _AA.V: 3.500
}



class SecondaryStructureMethod(StrEnum):
    """Enum for secondary structure prediction methods."""

    DELEAGE_ROUX = "DeleageRoux"
    LEVITT = "Levitt"
    CHOU_FASMAN = "ChouFasman"


class SecondaryStructureType(StrEnum):
    """Enum for secondary structure types."""

    ALPHA_HELIX = "alpha_helix"
    BETA_SHEET = "beta_sheet"
    BETA_TURN = "beta_turn"
    COIL = "coil"


class PropertyScale(StrEnum):
    """Base enum for all amino acid property scales."""

    pass


class HydrophobicityScale(PropertyScale):
    """Hydrophobicity scales."""

    KYTE_DOOLITTLE = "hphob_kyte_doolittle"
    ADOBERIN = "hphob_adoberin"
    ABRAHAM_LEO = "hphob_abraham_leo"
    AGROS = "hphob_agros"
    RAO_ARGOS = "hphob_rao_argos"
    BLACK_MOULD = "hphob_black_mould"
    BULL_BREESE = "hphob_bull_breese"
    CASARI_SIPPL = "hphob_casari_sippl"
    CID = "hphob_cid"
    COWAN_3_4 = "hphob_cowan_3_4"
    COWAN_7_5 = "hphob_cowan_7_5"
    EISENBERG = "hphob_eisenberg"
    ENGELMAN = "hphob_engelman"
    FASMAN = "hphob_fasman"
    FAUCHERE = "hphob_fauchere"
    GOLDSACK = "hphob_goldsack"
    GUY = "hphob_guy"
    JONES = "hphob_jones"
    JURETIC = "hphob_juretic"
    KIDERA = "hphob_kidera"
    MIYAZAWA = "hphob_miyazawa"
    PARKER = "hphob_parker"
    PONNUSWAMY = "hphob_ponnuswamy"
    MANAVALAN = "hphob_manavalan"
    ROSE = "hphob_rose"
    ROSEMAN = "hphob_roseman"
    SWEET = "hphob_sweet"
    TANFORD = "hphob_tanford"
    WILSON = "hphob_wilson"
    ZIMMERMAN = "hphob_zimmerman"
    CHOTHIA = "hphob_chothia"
    JANIN = "hphob_janin"
    WOLFENDEN = "hphob_wolfenden"
    WELLING = "hphob_welling"


class SecondaryStructureScale(PropertyScale):
    """Secondary structure prediction scales."""

    DELEAGE_ROUX_ALPHA_HELIX = "deleage_roux_alpha_helix"
    DELEAGE_ROUX_BETA_SHEET = "deleage_roux_beta_sheet"
    DELEAGE_ROUX_BETA_TURN = "deleage_roux_beta_turn"
    DELEAGE_ROUX_COIL = "deleage_roux_coil"
    LEVITT_ALPHA_HELIX = "levitt_alpha_helix"
    LEVITT_BETA_SHEET = "levitt_beta_sheet"
    LEVITT_BETA_TURN = "levitt_beta_turn"
    CHOU_FASMAN_ALPHA_HELIX = "chou_fasman_alpha_helix"
    CHOU_FASMAN_BETA_SHEET = "chou_fasman_beta_sheet"
    CHOU_FASMAN_BETA_TURN = "chou_fasman_beta_turn"


class SurfaceAccessibilityScale(PropertyScale):
    """Surface accessibility and area scales."""

    VERGOTEN = "surface_accessibility_vergoten"
    JANIN = "surface_accessiblility_janin"  # Note: keeping original typo
    ACCESSIBLE_RESIDUES = "accessible_residues"
    AVERAGE_BURIED_AREA = "average_buried_area"


class ChargeScale(PropertyScale):
    """Charge and pK related scales."""

    PK_NTERMINAL = "pk_nterminal"
    PK_CTERMINAL = "pk_cterminal"
    PK_SIDECHAIN = "pk_sidechain"


class PolarityScale(PropertyScale):
    """Polarity scales."""

    GRANTHAM = "polarity_grantham"
    ZIMMERMAN = "polarity_zimmerman"


class HPLCScale(PropertyScale):
    """HPLC retention time scales."""

    MEEK_2_1 = "hplc_meek_2_1"
    BROWNE = "hplc_browne"
    MEEK_7_4 = "hplc_meek_7_4"
    BROWNE_TFA = "hplc_browne_tfa"


class BetaStrandScale(PropertyScale):
    """Beta strand propensity scales."""

    PARALLEL = "beta_strand_parallel"
    ANTIPARALLEL = "beta_strand_antiparallel"
    TOTAL = "beta_strand_total"


class PhysicalPropertyScale(PropertyScale):
    """Physical and chemical property scales."""

    MOLECULAR_WEIGHTS = "molecular_weights"
    BULKINESS = "bulkiness"
    REFRACTIVITY = "refractivity"
    FLEXIBILITY_VIHINEN = "flexibility_vihinen"
    HYDROPHILICITY_HOP_WOOD = "hydrophilicity_hop_wood"
    RATIOSIDE = "ratioside"
    MUTABILITY = "mutability"
    CODONS = "codons"
    RECOGNITION_FACTORS = "recognition_factors"
    TRANSMEMBRANE_TENDENCY = "transmembrane_tendency"


class CompositionScale(PropertyScale):
    """Amino acid composition scales."""

    MCCALDRON = "aa_composition_mccaldron"
    SWISSPROT = "aa_composition_swissprot"


# Updated dictionaries using the enums
secondary_structure_scales: dict[str, dict[str, float]] = {
    SecondaryStructureScale.DELEAGE_ROUX_ALPHA_HELIX: deleage_roux_alpha_helix,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_SHEET: deleage_roux_beta_sheet,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_TURN: deleage_roux_beta_turn,
    SecondaryStructureScale.DELEAGE_ROUX_COIL: deleage_roux_coil,
    SecondaryStructureScale.LEVITT_ALPHA_HELIX: levitt_alpha_helix,
    SecondaryStructureScale.LEVITT_BETA_SHEET: levitt_beta_sheet,
    SecondaryStructureScale.LEVITT_BETA_TURN: levitt_beta_turn,
    SecondaryStructureScale.CHOU_FASMAN_ALPHA_HELIX: chou_fasman_alpha_helix,
    SecondaryStructureScale.CHOU_FASMAN_BETA_SHEET: chou_fasman_beta_sheet,
    SecondaryStructureScale.CHOU_FASMAN_BETA_TURN: chou_fasman_beta_turn,
}

secondary_structure_scales_by_name: dict[str, dict[str, dict[str, float]]] = {
    SecondaryStructureMethod.DELEAGE_ROUX: {
        SecondaryStructureType.ALPHA_HELIX: deleage_roux_alpha_helix,
        SecondaryStructureType.BETA_SHEET: deleage_roux_beta_sheet,
        SecondaryStructureType.BETA_TURN: deleage_roux_beta_turn,
        SecondaryStructureType.COIL: deleage_roux_coil,
    },
    SecondaryStructureMethod.LEVITT: {
        SecondaryStructureType.ALPHA_HELIX: levitt_alpha_helix,
        SecondaryStructureType.BETA_SHEET: levitt_beta_sheet,
        SecondaryStructureType.BETA_TURN: levitt_beta_turn,
    },
    SecondaryStructureMethod.CHOU_FASMAN: {
        SecondaryStructureType.ALPHA_HELIX: chou_fasman_alpha_helix,
        SecondaryStructureType.BETA_SHEET: chou_fasman_beta_sheet,
        SecondaryStructureType.BETA_TURN: chou_fasman_beta_turn,
    },
}



PROPERTY_SCALES: Final[dict[str, dict[str, float]]] = {
    # Hydrophobicity scales
    HydrophobicityScale.KYTE_DOOLITTLE: hphob_kyte_doolittle,
    HydrophobicityScale.ADOBERIN: hphob_adoberin,
    HydrophobicityScale.ABRAHAM_LEO: hphob_abraham_leo,
    HydrophobicityScale.AGROS: hphob_agros,
    HydrophobicityScale.RAO_ARGOS: hphob_rao_argos,
    HydrophobicityScale.BLACK_MOULD: hphob_black_mould,
    HydrophobicityScale.BULL_BREESE: hphob_bull_breese,
    HydrophobicityScale.CASARI_SIPPL: hphob_casari_sippl,
    HydrophobicityScale.CID: hphob_cid,
    HydrophobicityScale.COWAN_3_4: hphob_cowan_3_4,
    HydrophobicityScale.COWAN_7_5: hphob_cowan_7_5,
    HydrophobicityScale.EISENBERG: hphob_eisenberg,
    HydrophobicityScale.ENGELMAN: hphob_engelman,
    HydrophobicityScale.FASMAN: hphob_fasman,
    HydrophobicityScale.FAUCHERE: hphob_fauchere,
    HydrophobicityScale.GOLDSACK: hphob_goldsack,
    HydrophobicityScale.GUY: hphob_guy,
    HydrophobicityScale.JONES: hphob_jones,
    HydrophobicityScale.JURETIC: hphob_juretic,
    HydrophobicityScale.KIDERA: hphob_kidera,
    HydrophobicityScale.MIYAZAWA: hphob_miyazawa,
    HydrophobicityScale.PARKER: hphob_parker,
    HydrophobicityScale.PONNUSWAMY: hphob_ponnuswamy,
    HydrophobicityScale.MANAVALAN: hphob_manavalan,
    HydrophobicityScale.ROSE: hphob_rose,
    HydrophobicityScale.ROSEMAN: hphob_roseman,
    HydrophobicityScale.SWEET: hphob_sweet,
    HydrophobicityScale.TANFORD: hphob_tanford,
    HydrophobicityScale.WILSON: hphob_wilson,
    HydrophobicityScale.ZIMMERMAN: hphob_zimmerman,
    HydrophobicityScale.CHOTHIA: hphob_chothia,
    HydrophobicityScale.JANIN: hphob_janin,
    HydrophobicityScale.WOLFENDEN: hphob_wolfenden,
    HydrophobicityScale.WELLING: hphob_welling,
    # Secondary structure scales
    SecondaryStructureScale.DELEAGE_ROUX_ALPHA_HELIX: deleage_roux_alpha_helix,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_SHEET: deleage_roux_beta_sheet,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_TURN: deleage_roux_beta_turn,
    SecondaryStructureScale.DELEAGE_ROUX_COIL: deleage_roux_coil,
    SecondaryStructureScale.LEVITT_ALPHA_HELIX: levitt_alpha_helix,
    SecondaryStructureScale.LEVITT_BETA_SHEET: levitt_beta_sheet,
    SecondaryStructureScale.LEVITT_BETA_TURN: levitt_beta_turn,
    SecondaryStructureScale.CHOU_FASMAN_ALPHA_HELIX: chou_fasman_alpha_helix,
    SecondaryStructureScale.CHOU_FASMAN_BETA_SHEET: chou_fasman_beta_sheet,
    SecondaryStructureScale.CHOU_FASMAN_BETA_TURN: chou_fasman_beta_turn,
    # Surface accessibility scales
    SurfaceAccessibilityScale.VERGOTEN: surface_accessibility_vergoten,
    SurfaceAccessibilityScale.JANIN: surface_accessiblility_janin,
    SurfaceAccessibilityScale.ACCESSIBLE_RESIDUES: accessible_residues,
    SurfaceAccessibilityScale.AVERAGE_BURIED_AREA: average_buried_area,
    # Charge scales
    ChargeScale.PK_NTERMINAL: pk_nterminal,
    ChargeScale.PK_CTERMINAL: pk_cterminal,
    ChargeScale.PK_SIDECHAIN: pk_sidechain,
    # Polarity scales
    PolarityScale.GRANTHAM: polarity_grantham,
    PolarityScale.ZIMMERMAN: polarity_zimmerman,
    # HPLC scales
    HPLCScale.MEEK_2_1: hplc_meek_2_1,
    HPLCScale.BROWNE: hplc_browne,
    HPLCScale.MEEK_7_4: hplc_meek_7_4,
    HPLCScale.BROWNE_TFA: hplc_browne_tfa,
    # Beta strand scales
    BetaStrandScale.PARALLEL: beta_strand_parallel,
    BetaStrandScale.ANTIPARALLEL: beta_strand_antiparallel,
    BetaStrandScale.TOTAL: beta_strand_total,
    # Physical property scales
    PhysicalPropertyScale.MOLECULAR_WEIGHTS: molecular_weights,
    PhysicalPropertyScale.BULKINESS: bulkiness,
    PhysicalPropertyScale.REFRACTIVITY: refractivity,
    PhysicalPropertyScale.FLEXIBILITY_VIHINEN: flexibility_vihinen,
    PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD: hydrophilicity_hopp_wood,
    PhysicalPropertyScale.RATIOSIDE: ratioside,
    PhysicalPropertyScale.MUTABILITY: mutability,
    PhysicalPropertyScale.CODONS: codons,
    PhysicalPropertyScale.RECOGNITION_FACTORS: recognition_factors,
    PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY: transmembrane_tendency,
    # Composition scales
    CompositionScale.MCCALDRON: aa_composition_mccaldron,
    CompositionScale.SWISSPROT: aa_composition_swissprot,
}

HYDROPHOBICITY_SCALES: Final[dict[str, dict[str, float]]] = {
    HydrophobicityScale.KYTE_DOOLITTLE: hphob_kyte_doolittle,
    HydrophobicityScale.ADOBERIN: hphob_adoberin,
    HydrophobicityScale.ABRAHAM_LEO: hphob_abraham_leo,
    HydrophobicityScale.AGROS: hphob_agros,
    HydrophobicityScale.BLACK_MOULD: hphob_black_mould,
    HydrophobicityScale.BULL_BREESE: hphob_bull_breese,
    HydrophobicityScale.CASARI_SIPPL: hphob_casari_sippl,
    HydrophobicityScale.CID: hphob_cid,
    HydrophobicityScale.COWAN_3_4: hphob_cowan_3_4,
    HydrophobicityScale.COWAN_7_5: hphob_cowan_7_5,
    HydrophobicityScale.EISENBERG: hphob_eisenberg,
    HydrophobicityScale.ENGELMAN: hphob_engelman,
    HydrophobicityScale.FASMAN: hphob_fasman,
    HydrophobicityScale.FAUCHERE: hphob_fauchere,
    HydrophobicityScale.GOLDSACK: hphob_goldsack,
    HydrophobicityScale.GUY: hphob_guy,
    HydrophobicityScale.JONES: hphob_jones,
    HydrophobicityScale.JURETIC: hphob_juretic,
    HydrophobicityScale.KIDERA: hphob_kidera,
    HydrophobicityScale.MIYAZAWA: hphob_miyazawa,
    HydrophobicityScale.PARKER: hphob_parker,
    HydrophobicityScale.PONNUSWAMY: hphob_ponnuswamy,
    HydrophobicityScale.ROSE: hphob_rose,
    HydrophobicityScale.ROSEMAN: hphob_roseman,
    HydrophobicityScale.SWEET: hphob_sweet,
    HydrophobicityScale.TANFORD: hphob_tanford,
    HydrophobicityScale.WILSON: hphob_wilson,
    HydrophobicityScale.ZIMMERMAN: hphob_zimmerman,
    HydrophobicityScale.RAO_ARGOS: hphob_rao_argos,
    HydrophobicityScale.MANAVALAN: hphob_manavalan,
    HydrophobicityScale.CHOTHIA: hphob_chothia,
    HydrophobicityScale.JANIN: hphob_janin,
    HydrophobicityScale.WOLFENDEN: hphob_wolfenden,
    HydrophobicityScale.WELLING: hphob_welling,
}


# Surface accessibility scales
SURFACE_ACCESSIBILITY_SCALES: Final[dict[str, dict[str, float]]] = {
    SurfaceAccessibilityScale.VERGOTEN: surface_accessibility_vergoten,
    SurfaceAccessibilityScale.JANIN: surface_accessiblility_janin,
    SurfaceAccessibilityScale.ACCESSIBLE_RESIDUES: accessible_residues,
    SurfaceAccessibilityScale.AVERAGE_BURIED_AREA: average_buried_area,
}


HPLC_SCALES: Final[dict[str, dict[str, float]]] = {
    HPLCScale.MEEK_2_1: hplc_meek_2_1,
    HPLCScale.BROWNE: hplc_browne,
    HPLCScale.MEEK_7_4: hplc_meek_7_4,
    HPLCScale.BROWNE_TFA: hplc_browne_tfa,
}


# Hydrophilicity scales
HYDROPHILICITY_SCALES: Final[dict[str, dict[str, float]]] = {
    PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD: hydrophilicity_hopp_wood
}

# Flexibility scales
FLEXIBILITY_SCALES: Final[dict[str, dict[str, float]]] = {
    PhysicalPropertyScale.FLEXIBILITY_VIHINEN: flexibility_vihinen
}

# Polarity scales (enum-keyed)
POLARITY_SCALES: Final[dict[str, dict[str, float]]] = {
    PolarityScale.GRANTHAM: polarity_grantham,
    PolarityScale.ZIMMERMAN: polarity_zimmerman,
}

# Composition scales (enum-keyed)
COMPOSITION_SCALES: Final[dict[str, dict[str, float]]] = {
    CompositionScale.MCCALDRON: aa_composition_mccaldron,
    CompositionScale.SWISSPROT: aa_composition_swissprot,
}

# Physical property scales — exhaustive mapping using PhysicalPropertyScale
PHYSICAL_PROPERTY_SCALES: Final[dict[str, dict[str, float]]] = {
    PhysicalPropertyScale.MOLECULAR_WEIGHTS: molecular_weights,
    PhysicalPropertyScale.BULKINESS: bulkiness,
    PhysicalPropertyScale.REFRACTIVITY: refractivity,
    PhysicalPropertyScale.FLEXIBILITY_VIHINEN: flexibility_vihinen,
    PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD: hydrophilicity_hopp_wood,
    PhysicalPropertyScale.RATIOSIDE: ratioside,
    PhysicalPropertyScale.MUTABILITY: mutability,
    PhysicalPropertyScale.CODONS: codons,
    PhysicalPropertyScale.RECOGNITION_FACTORS: recognition_factors,
    PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY: transmembrane_tendency,
}
