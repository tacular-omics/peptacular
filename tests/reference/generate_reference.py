"""Generate ``reference_data.json``: values from tools independent of peptacular.

peptacular is never imported here. Every value comes from one of:

* pyteomics 5.0.1 (``mass.calculate_mass``, ``mass.Composition``, ``mass.std_ion_comp``,
  ``parser.cleave`` with ``expasy_rules``/``psims_rules``, ``electrochem.charge``/``pI``,
  ``proforma.ProForma``), which uses NIST element masses and abundances.
* An exact isotope-cluster convolution written below. Abundances are pyteomics' NIST
  abundances; isotope masses are AME2020 values typed in below, because pyteomics 5.0.1
  ships older masses (e.g. 18O 17.999161 vs 17.99915961286) that shift cluster
  centres by up to 3e-6 Da.
* The peptideweb pKa table, typed in from
  https://www.peptideweb.com/images/pdf/pKa-and-pI-values-of-amino-acids.pdf
* ExPASy ProtScale amino-acid scale pages (https://web.expasy.org/protscale/pscale/),
  fetched at generation time, and Biopython 1.88 ``Bio.SeqUtils.ProtParamData``.
* Worked examples from the ProForma 2.0 specification (HUPO-PSI, 2022), with masses
  from pyteomics.proforma where it resolves them, or from an explicit composition.

Produced with (network needed for ProtScale; biopython is not a peptacular dependency):

    uv run --with biopython==1.88 python tests/reference/generate_reference.py

Generated 2026-09-23 with pyteomics 5.0.1, biopython 1.88, Python 3.12.
"""

from __future__ import annotations

import importlib.metadata
import json
import re
import urllib.request
import warnings
from pathlib import Path

from pyteomics import electrochem as pe
from pyteomics import mass as pm
from pyteomics import parser as pp
from pyteomics import proforma as ppf

warnings.filterwarnings("ignore")

OUT = Path(__file__).with_name("reference_data.json")

PEPTIDES = ["G", "GA", "PEPTIDE", "SAMPLER", "KRKRHH", "WYWYCC", "MVIMSEFSADPAGQGQGQQK", "ACDEFGHIKLMNPQRSTVWY"]

# Unimod compositions (Unimod ids 1, 4, 35, 21, 737), typed in from unimod.org.
MOD_COMPS = {
    "Acetyl": pm.Composition(formula="C2H2O"),
    "Carbamidomethyl": pm.Composition(formula="C2H3NO"),
    "Oxidation": pm.Composition(formula="O"),
    "Phospho": pm.Composition(formula="HPO3"),
    "TMT6plex": pm.Composition({"C": 8, "C[13]": 4, "H": 20, "N": 1, "N[15]": 1, "O": 2}),
}

# (proforma, unmodified sequence, mods added to the composition)
MODIFIED = [
    ("PEPTC[Carbamidomethyl]IDE", "PEPTCIDE", ["Carbamidomethyl"]),
    ("[Acetyl]-PEPTIDE", "PEPTIDE", ["Acetyl"]),
    ("PEPTIDE-[Oxidation]", "PEPTIDE", ["Oxidation"]),
    ("PEPS[Phospho]T[Phospho]IDE", "PEPSTIDE", ["Phospho", "Phospho"]),
    ("[TMT6plex]-PEPTK[TMT6plex]IDE", "PEPTKIDE", ["TMT6plex", "TMT6plex"]),
    ("M[Oxidation]", "M", ["Oxidation"]),
    ("{Phospho}PEPSTIDE", "PEPSTIDE", ["Phospho"]),
    ("[Phospho]?PEPSTIDE", "PEPSTIDE", ["Phospho"]),
    ("<[Carbamidomethyl]@C>PEPTCCIDE", "PEPTCCIDE", ["Carbamidomethyl"] * 2),
    ("<[Acetyl]@N-term>PEPTIDE", "PEPTIDE", ["Acetyl"]),
    ("PEPT[U:Phospho]IDE", "PEPTIDE", ["Phospho"]),
    ("PEPT[UNIMOD:21]IDE", "PEPTIDE", ["Phospho"]),
    ("PEPT[MOD:00046]IDE", "PEPTIDE", ["Phospho"]),
    ("PEPT[Formula:HPO3]IDE", "PEPTIDE", ["Phospho"]),
]

# peptacular ion name -> pyteomics std_ion_comp name
TERMINAL_IONS = {"a": "a", "b": "b", "c": "c", "x": "x", "y": "y", "z": "z", "z.": "z-dot", "z+H": "z+2", "c-H": "c-1"}

# Internal ions: residues - H2O + N-side offset (a: -CO, b: 0, c: +NH3) + C-side offset
# (x: +CO-H2, y: 0, z: -NH3). Same offsets as the terminal ions of the same letters.
INTERNAL_N = {"a": "C-1O-1", "b": "", "c": "NH3"}
INTERNAL_C = {"x": "COH-2", "y": "", "z": "N-1H-3"}
INTERNAL_TYPES = ["by", "ax", "cz", "ay", "az", "bx", "bz", "cx", "cy"]

PROTEIN = (
    "MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTEYAKNGEMF"
    "DYLTSHGHLSEKEARKKFWQILSAVEYCHDHKIVHRDLKPENLLLDDNMNIKIADFGFSNLFTPGQLLKTWCGSPPYAAPELFEGKEYDGPKVDIWSLGVVLYVLVCGALPF"
    "DGSTLQNLRARVLSGKFRIPFFMSTECEHLIRHMLVLDPNKRLSMEQICKHKWMKLGDADPNFDRLIAKNVPGPNRLG"
)
# peptacular protease name -> (pyteomics rule table, rule name)
PROTEASES = {
    "trypsin": ("expasy_rules", "trypsin"),
    "trypsin_full": ("psims_rules", "Trypsin/P"),
    "lys_c": ("expasy_rules", "lysc"),
    "asp_n": ("expasy_rules", "asp-n"),
    "chymotrypsin": ("psims_rules", "Chymotrypsin"),
    "arg_c": ("expasy_rules", "arg-c"),
    "glu_c": ("expasy_rules", "glutamyl endopeptidase"),
    "proteinase_k": ("expasy_rules", "proteinase k"),
}

# peptideweb table: residue -> (pKa COOH, pKa NH3+, pKa side chain or None)
PKA_TABLE = {
    "A": (2.34, 9.69, None), "R": (2.17, 9.04, 12.48), "N": (2.02, 8.80, None), "D": (2.09, 9.82, 3.86),
    "C": (1.71, 10.78, 8.33), "E": (2.19, 9.67, 4.25), "Q": (2.17, 9.13, None), "G": (2.34, 9.60, None),
    "H": (1.82, 9.17, 6.00), "I": (2.36, 9.60, None), "L": (2.36, 9.60, None), "K": (2.18, 8.95, 10.79),
    "M": (2.28, 9.21, None), "F": (1.83, 9.13, None), "P": (1.99, 10.60, None), "S": (2.21, 9.15, None),
    "T": (2.09, 9.10, None), "W": (2.43, 9.44, None), "Y": (2.20, 9.11, 10.07), "V": (2.32, 9.62, None),
}  # fmt: skip
PI_PEPTIDES = ["PEPTIDE", "KRKRHH", "ACDEFGHIKLMNPQRSTVWY", "DDDEEE", "CCYY", "G", "HHHH", "SAMPLEQ", "PEPTIDEE"]

# AME2020 isotope masses (Wang et al., Chinese Phys. C 45, 030003 (2021)).
AME2020 = {
    "H": {1: 1.00782503223, 2: 2.01410177812},
    "C": {12: 12.0, 13: 13.00335483507},
    "N": {14: 14.00307400443, 15: 15.00010889888},
    "O": {16: 15.99491461957, 17: 16.99913175650, 18: 17.99915961286},
    "S": {32: 31.9720711744, 33: 32.9714589098, 34: 33.967867004, 36: 35.96708071},
}
ISOTOPE_PEPTIDES = ["G", "PEPTIDE", "SAMPLER", "MCMCM", "ACDEFGHIK", "ACDEFGHIKLMNPQRSTVWY" * 2]
ISOTOPE_PEAKS = 6

# peptacular scale id -> (source, table name). ProtScale names are page names under
# https://web.expasy.org/protscale/pscale/<name>.html; Biopython names are attributes of
# Bio.SeqUtils.ProtParamData. The pKa scales are checked against PKA_TABLE instead.
SCALES = {
    "hphob_kyte_doolittle": ("protscale", "Hydropath.Doolittle"),
    "hphob_aboderin": ("protscale", "Hphob.mobility"),
    "hphob_abraham_leo": ("protscale", "Hphob.Leo"),
    "hphob_argos": ("biopython", "ag"),
    "hphob_rao_argos": ("protscale", "Hphob.Argos"),
    "hphob_black_mould": ("protscale", "Hphob.Black"),
    "hphob_bull_breese": ("protscale", "Hphob.Breese"),
    "hphob_casari_sippl": ("biopython", "cs"),
    "hphob_cid": ("biopython", "ci"),
    "hphob_cowan_3_4": ("protscale", "Hphob.pH3.4"),
    "hphob_cowan_7_5": ("protscale", "Hphob.pH7.5"),
    "hphob_eisenberg": ("protscale", "Hphob.Eisenberg"),
    "hphob_engelman": ("biopython", "eg"),
    "hphob_fasman": ("biopython", "fs"),
    "hphob_fauchere": ("protscale", "Hphob.Fauchere"),
    "hphob_goldsack": ("biopython", "gd"),
    "hphob_guy": ("protscale", "Hphob.Guy"),
    "hphob_jones": ("biopython", "jo"),
    "hphob_juretic": ("biopython", "ju"),
    "hphob_kidera": ("biopython", "ki"),
    "hphob_miyazawa": ("protscale", "Hphob.Miyazawa"),
    "hphob_parker": ("protscale", "Hphob.Parker"),
    "hphob_ponnuswamy": ("biopython", "po"),
    "hphob_manavalan": ("protscale", "Hphob.Manavalan"),
    "hphob_rose": ("protscale", "Hphob.Rose"),
    "hphob_roseman": ("protscale", "Hphob.Roseman"),
    "hphob_sweet": ("protscale", "Hphob.Sweet"),
    "hphob_tanford": ("protscale", "Hphob.Tanford"),
    "hphob_wilson": ("protscale", "Hphob.Wilson"),
    "hphob_zimmerman": ("biopython", "zi"),
    "hphob_chothia": ("protscale", "Hphob.Chothia"),
    "hphob_janin": ("protscale", "Hphob.Janin"),
    "hphob_wolfenden": ("protscale", "Hphob.Wolfenden"),
    "hphob_welling": ("protscale", "Hphob.Welling"),
    "deleage_roux_alpha_helix": ("protscale", "alpha-helixRoux"),
    "deleage_roux_beta_sheet": ("protscale", "beta-sheetRoux"),
    "deleage_roux_beta_turn": ("protscale", "beta-turnRoux"),
    "deleage_roux_coil": ("protscale", "CoilRoux"),
    "levitt_alpha_helix": ("protscale", "alpha-helixLevitt"),
    "levitt_beta_sheet": ("protscale", "beta-sheetLevitt"),
    "levitt_beta_turn": ("protscale", "beta-turnLevitt"),
    "chou_fasman_alpha_helix": ("protscale", "alpha-helixFasman"),
    "chou_fasman_beta_sheet": ("protscale", "beta-sheetFasman"),
    "chou_fasman_beta_turn": ("protscale", "beta-turnFasman"),
    "surface_accessibility_vergoten": ("biopython", "em"),
    "surface_accessibility_janin": ("biopython", "ja"),
    "accessible_residues": ("protscale", "accessibleresidues"),
    "average_buried_area": ("protscale", "Averageburied"),
    "polarity_grantham": ("protscale", "PolarityGrantham"),
    "polarity_zimmerman": ("protscale", "PolarityZimmerman"),
    "hplc_meek_2_1": ("protscale", "HPLC2.1"),
    "hplc_browne": ("protscale", "HPLCHFBA"),
    "hplc_meek_7_4": ("protscale", "HPLC7.4"),
    "hplc_browne_tfa": ("protscale", "HPLCTFA"),
    "beta_strand_parallel": ("protscale", "Parallelbeta-strand"),
    "beta_strand_antiparallel": ("protscale", "Antiparallelbeta-strand"),
    "beta_strand_total": ("protscale", "Totalbeta-strand"),
    "molecular_weights": ("protscale", "Molecularweight"),
    "bulkiness": ("protscale", "Bulkiness"),
    "refractivity": ("protscale", "Refractivity"),
    "flexibility_vihinen": ("biopython", "Flex"),
    "hydrophilicity_hop_wood": ("protscale", "Hphob.Woods"),
    "ratioside": ("protscale", "Ratioside"),
    "mutability": ("protscale", "Relativemutability"),
    "codons": ("protscale", "Numbercodons"),
    "recognition_factors": ("protscale", "Recognitionfactors"),
    "transmembrane_tendency": ("protscale", "Transmembranetendency"),
    "aa_composition_mccaldron": ("protscale", "A.A.composition"),
    "aa_composition_swissprot": ("protscale", "A.A.Swiss-Prot"),
}

# Worked examples from the ProForma 2.0 specification, one per line, in spec order.
SPEC_EXAMPLES = Path(__file__).with_name("proforma_spec_examples.txt")

PROTON = pm.nist_mass["H+"][0][0]
THREE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}  # fmt: skip


def comp_mass(comp: pm.Composition, charge: int = 0) -> float:
    return pm.calculate_mass(composition=comp, charge=charge or None)


def masses() -> dict:
    out: dict = {"unmodified": {}, "modified": {}}
    for s in PEPTIDES:
        out["unmodified"][s] = {
            "mono": pm.calculate_mass(sequence=s),
            "average": pm.calculate_mass(sequence=s, average=True),
            "mz": {z: pm.calculate_mass(sequence=s, charge=z) for z in (1, 2, 3, 4)},
        }
    for proforma, seq, mods in MODIFIED:
        comp = pm.Composition(sequence=seq)
        for m in mods:
            comp += MOD_COMPS[m]
        out["modified"][proforma] = {"mono": comp_mass(comp), "mz": {z: comp_mass(comp, z) for z in (1, 2, 3, 4)}}
    return out


# Modified peptides whose b/y ions are checked: proforma -> (sequence, {0-based index: mod}).
MODIFIED_FRAGMENTS = {
    "PEPTC[Carbamidomethyl]IDE": ("PEPTCIDE", {4: "Carbamidomethyl"}),
    "PEPS[Phospho]TM[Oxidation]IDEK": ("PEPSTMIDEK", {3: "Phospho", 5: "Oxidation"}),
}


def fragments() -> dict:
    out: dict = {"terminal": {}, "immonium": {}, "internal": {}, "neutral_loss": {}, "modified": {}}
    for s in PEPTIDES[1:]:
        n = len(s)
        out["terminal"][s] = {}
        for ion, py_ion in TERMINAL_IONS.items():
            for z in (1, 2, 3):
                series = []
                for i in range(1, n + 1):
                    sub = s[:i] if ion[0] in "abc" else s[n - i :]
                    series.append(pm.calculate_mass(sequence=sub, ion_type=py_ion, charge=z))
                out["terminal"][s][f"{ion}^{z}"] = series
        # Immonium: the residue minus CO, protonated (the a1 ion of that residue).
        out["immonium"][s] = [pm.calculate_mass(sequence=aa, ion_type="a", charge=1) for aa in s]
        out["internal"][s] = {}
        for it in INTERNAL_TYPES:
            nterm = pm.Composition(formula=INTERNAL_N[it[0]]) if INTERNAL_N[it[0]] else pm.Composition()
            cterm = pm.Composition(formula=INTERNAL_C[it[1]]) if INTERNAL_C[it[1]] else pm.Composition()
            rows = []
            for start in range(1, n + 1):
                for end in range(start, n + 1):
                    comp = pm.Composition(sequence=s[start - 1 : end]) - pm.Composition(formula="H2O") + nterm + cterm
                    rows.append([start, end, comp_mass(comp, 1)])
            out["internal"][s][it] = rows
    for s in ["PEPTIDE", "SAMPLER", "KRNQST"]:
        n = len(s)
        out["neutral_loss"][s] = {}
        for loss in ("H2O", "NH3"):
            for ion in ("b", "y"):
                out["neutral_loss"][s][f"{ion}-{loss}"] = [
                    pm.calculate_mass(sequence=s[:i] if ion == "b" else s[n - i :], ion_type=ion, charge=1) - pm.calculate_mass(formula=loss)
                    for i in range(1, n + 1)
                ]
    for proforma, (seq, mods) in MODIFIED_FRAGMENTS.items():
        n = len(seq)
        out["modified"][proforma] = {}
        for ion in ("b", "y"):
            for z in (1, 2):
                series = []
                for i in range(1, n + 1):
                    lo, hi = (0, i) if ion == "b" else (n - i, n)
                    comp = pm.Composition(sequence=seq[lo:hi], ion_type=ion)
                    for idx, mod in mods.items():
                        if lo <= idx < hi:
                            comp += MOD_COMPS[mod]
                    series.append(comp_mass(comp, z))
                out["modified"][proforma][f"{ion}^{z}"] = series
    return out


def exact_distribution(comp: pm.Composition, nmax: int) -> list[list[float]]:
    """Aggregated isotope cluster by exact polynomial convolution.

    Each peak is [neutron offset, abundance relative to the most abundant peak,
    abundance-weighted mean mass]. Abundances: pyteomics NIST. Masses: AME2020.
    """
    dist = {0: (1.0, 0.0)}
    for el, cnt in comp.items():
        isos = [(AME2020[el][k], a) for k, (_, a) in pm.nist_mass[el].items() if k and a > 0]
        m0 = min(m for m, _ in isos)
        single = {round(m - m0): (a, a * m) for m, a in isos}
        for _ in range(cnt):
            new: dict[int, tuple[float, float]] = {}
            for o1, (p1, w1) in dist.items():
                for o2, (p2, w2) in single.items():
                    o = o1 + o2
                    if o > nmax:
                        continue
                    p, w = new.get(o, (0.0, 0.0))
                    new[o] = (p + p1 * p2, w + w1 * p2 + p1 * w2)
            dist = new
    top = max(p for p, _ in dist.values())
    return [[o, p / top, w / p] for o, (p, w) in sorted(dist.items())]


def isotopes() -> dict:
    return {s: exact_distribution(pm.Composition(sequence=s), ISOTOPE_PEAKS) for s in ISOTOPE_PEPTIDES}


def digestion() -> dict:
    out: dict = {"protein": PROTEIN, "proteases": {}}
    for name, (table, rule) in PROTEASES.items():
        regex = getattr(pp, table)[rule]
        cases = {}
        for mc in (0, 1, 2):
            for semi in (False, True):
                if semi and mc == 2:
                    continue
                cases[f"mc{mc}{'_semi' if semi else ''}"] = sorted(pp.cleave(PROTEIN, regex, missed_cleavages=mc, semi=semi))
        out["proteases"][name] = {"source": f"pyteomics.parser.{table}[{rule!r}]", "regex": regex, "peptides": cases}
    return out


def charge_and_pi() -> dict:
    out: dict = {"pka_table": PKA_TABLE, "peptides": {}}
    side = {aa: [(v[2], 1 if aa in "RHK" else -1)] for aa, v in PKA_TABLE.items() if v[2] is not None}
    for s in PI_PEPTIDES:
        pk = dict(side)
        pk["H-"] = [(PKA_TABLE[s[0]][1], 1)]
        pk["-OH"] = [(PKA_TABLE[s[-1]][0], -1)]
        out["peptides"][s] = {
            "charge": {ph: pe.charge(s, ph, pK=pk) for ph in (2.0, 4.0, 7.0, 10.0, 12.0)},
            "pi": pe.pI(s, pK=pk, precision_pI=1e-8),
        }
    return out


def fetch_protscale(name: str) -> dict[str, float]:
    url = f"https://web.expasy.org/protscale/pscale/{name}.html"
    with urllib.request.urlopen(url, timeout=30) as resp:
        text = re.sub(r"<[^>]*>", "", resp.read().decode("latin-1"))
    pattern = r"\b(" + "|".join(THREE) + r"):\s*(-?[\d.]+)"
    values = {THREE[m.group(1)]: float(m.group(2)) for m in re.finditer(pattern, text)}
    if len(values) != 20:
        raise ValueError(f"{url}: parsed {len(values)} residues")
    return values


def scales() -> dict:
    from Bio.SeqUtils import ProtParamData

    out = {}
    for scale_id, (source, name) in SCALES.items():
        if source == "protscale":
            values = fetch_protscale(name)
            ref = f"https://web.expasy.org/protscale/pscale/{name}.html"
        else:
            table = getattr(ProtParamData, name, None)
            if table is None:
                table = ProtParamData.gravy_scales[name]
            values = {aa: float(v) for aa, v in table.items() if aa in THREE.values()}
            ref = f"Bio.SeqUtils.ProtParamData.{name} (biopython {importlib.metadata.version('biopython')})"
        out[scale_id] = {"source": ref, "values": values}
    return out


# Monosaccharide residue formulas (ProForma 2.0 spec, Table of monosaccharides).
MONOSACCHARIDES = {"HexNAc": "C8H13NO5", "NeuAc": "C11H17NO8", "Hex": "C6H10O5"}
GLYCAN_RE = re.compile(r"Glycan:((?:(?:HexNAc|NeuAc|Hex)\d*)+)\]")
ISOTOPE_LABELS = {"13C": ("C", 13), "15N": ("N", 15), "D": ("H", 2)}


def _pyteomics_monosaccharide(name: str) -> float:
    return ppf.ProForma.parse(f"G[Glycan:{name}]").mass - ppf.ProForma.parse("G").mass


def _glycan_correction(line: str) -> float:
    """pyteomics' monosaccharide masses differ from the formulas by up to 7e-5 Da; swap them."""
    delta = 0.0
    for block in GLYCAN_RE.findall(line):
        for name, count in re.findall(r"(HexNAc|NeuAc|Hex)(\d*)", block):
            n = int(count or 1)
            delta += n * (pm.calculate_mass(formula=MONOSACCHARIDES[name]) - _pyteomics_monosaccharide(name))
    for name in re.findall(r"\{Glycan:(HexNAc|NeuAc|Hex)\}", line):
        delta += pm.calculate_mass(formula=MONOSACCHARIDES[name]) - _pyteomics_monosaccharide(name)
    return delta


def _isotope_labelled_mass(line: str) -> float | None:
    """pyteomics.proforma ignores <13C>/<15N>/<D>; label every atom of that element."""
    labels = re.findall(r"<(13C|15N|D)>", line)
    if not labels:
        return None
    comp = pm.Composition(sequence=re.sub(r"<[^>]*>", "", line))
    swap = dict(ISOTOPE_LABELS[lab] for lab in labels)
    return sum(cnt * pm.nist_mass[el][swap[el] if el in swap else 0][0] for el, cnt in comp.items())


def spec_examples() -> list[dict]:
    """Each spec example with a neutral-mass reference, or pyteomics' error.

    ``reference_mass`` is pyteomics.proforma's mass, corrected for the two things it does
    differently from the spec: monosaccharide masses (typed in above from their formulas)
    and global isotope labels (which it ignores).
    """
    rows = []
    for line in SPEC_EXAMPLES.read_text().splitlines():
        if not line.strip():
            continue
        row: dict = {"proforma": line, "reference_mass": None, "source": None, "pyteomics_error": None}
        labelled = _isotope_labelled_mass(line)
        if labelled is not None:
            row["reference_mass"], row["source"] = labelled, "pyteomics composition, every atom of the labelled element heavy"
        else:
            try:
                mass = ppf.ProForma.parse(line).mass
            except Exception as exc:  # noqa: BLE001 - recorded, not raised
                row["pyteomics_error"] = f"{type(exc).__name__}: {str(exc)[:120]}"
            else:
                correction = _glycan_correction(line)
                row["reference_mass"] = mass + correction
                row["source"] = "pyteomics.proforma" + (" + monosaccharide formula masses" if correction else "")
        rows.append(row)
    return rows


def main() -> None:
    data = {
        "_about": {
            "generator": "tests/reference/generate_reference.py",
            "versions": {
                "pyteomics": importlib.metadata.version("pyteomics"),
                "biopython": importlib.metadata.version("biopython"),
            },
            "proton": PROTON,
        },
        "masses": masses(),
        "fragments": fragments(),
        "isotopes": isotopes(),
        "digestion": digestion(),
        "charge_pi": charge_and_pi(),
        "scales": scales(),
        "proforma_spec": spec_examples(),
    }
    OUT.write_text(json.dumps(data, indent=1, sort_keys=False) + "\n")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
