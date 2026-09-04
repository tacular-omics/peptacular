# Notes on Compliance levels

### Base Handling
- **Controlled vocabularies**: Unimod, PSI-MOD, RESID, XL-MOD, and GNO modifications
  - Only supports modifications with unambiguous mass values
- **Delta masses**: `PEM[+15.995]AT` 
    - Composition cannot be determined
- **Labile modifications**: `{Glycan:Hex}EM[U:Oxidation]EV`
  - Only applied to precursor ('p') or neutral ('n') fragment ion types
- **Information tags**: `ELV[INFO:AnyString]IS` 
    - Preserved but not used in calculations. Treated as 0.0 mass.
- **Glycan compositions**: `NEEYN[Glycan:Hex{H2O}5HexNAc4NeuAc1]K`
  - Mixed named, molecular-formula, and bare-mass components are supported. Bare-mass components support mass calculation but do not have an elemental composition.

### Ambiguity Handling
- **Mass gaps**: `PEX[+147.035]AT`
    - Composition cannot be determined
- **Ambiguous amino acids**: B, Z, J, X
  - **J**: Sequence-ambiguous but mass-unambiguous (all calculations work)
  - **B, Z**: Mass-ambiguous (will raise errors for mass/composition/fragment calculations)
- **Amino acid ambiguity sets**: `(?VCH)AT`
  - Mass/mz/composition calculations work, but fragmentation fails (no defined localization)

### Localization Variants
- **Unknown position**: `[Oxidation]?PEMAT`
  - Mass/mz/composition work, fragmentation fails
- **Position sets / scores**: `PEP[Oxidation#1]M[#1]AT` or `PEP[Oxidation#1(0.95)]M[#1(0.05)]AT`
  - Always uses the declared location (first occurrence with modification specified)
- **Position ranges**: `PRT(ESFRMS)[Oxidation]ISK`
  - Mass/mz/composition work (Assuming modification has valid composition), fragmentation fails


### Isotope Labeling
```
<13C>PEPTIDE
```
- Mass calculation flow:
  1. Collect amino acid and modification composition, then apply the ion formula offset.
  2. Apply user isotope substitutions and elemental deltas to the complete neutral ion.
  3. Apply global isotope substitution.
  4. Apply external charge carriers, which are unaffected by global isotope substitution.
  5. Correct electron mass using intrinsic plus external charge. Mass-only deltas remain additive.

Global isotope substitutions apply to modification and delta formulas as well as the sequence and ion offset.

- Calculating Composition:
    - All modifications must have a known/valid composition.

- Calcualting Mass/Mz:
    - Allows for DeltaMass Tags, despite not having a valid composition. Since these are typically user provided, or represent an observed mass shift, they are taken as is (The delta mass is assumed to already have been calculated with the isotopic shift in mind).
    - All other modification types must have a valid composition.

- Mass and composition calculation paths use the same elemental adjustment order.

### Fixed Modifications
```
<[Oxidation]@M>ATPEMILTCMGCLK
```
Automatically mapped to applicable residues during mass calculation and fragmentation.

### Chimeric Spectra
```
NEEYN+SEQUEN
```
- Use `pt.parse_chimeric()` → returns list of ProForma objects
- Use `pt.serialize_chimeric()` → converts list back to chimeric notation
  - **Requirement**: All entries must share compound name, isotope mods, and static mods

## Charge State Notation
```
SEQUEN/2
SEQUEN/[Na:z+1,H:z+1]
```

Explicit adducts contribute their atom counts and charge. Fragment objects distinguish external carrier charge from intrinsic modification charge. Impossible deprotonation or elemental loss counts raise errors.

### Charged Formulas
```
SEQUEN[Formula:Zn1:z+2]CE
```
- Annotation-level charge and adducts only reflect end-of-sequence tags
- Internal charged formulas affect mz and fragment ion calculations
- Fragment objects reflect internal charges, but adducts show only user-supplied / end-of-sequence tags values. 
- `fast_fragment()` uses regular calculations when modifications carry intrinsic charge

## Interpretation Hints

### Mass with Interpretation
```
PEM[+15.995|Oxidation]AT
```
Only the first modification tag is used - treated as delta mass.

### Prefixed Delta Masses
```
PEM[U:+15.995]AT
```
Currently treated as delta mass. Mass-based lookup was attempted but abandoned due to conflicts from modifications with identical masses.

### Placement Control
```
PTI(MERMERME)[+32|Position:E]PTIDE
```
Parses correctly but not utilized by peptacular. Accessible if needed for external tools.

## Not Supported

### Cross-Linking
- ❌ Intrachain cross-linkers: `EVTK[X:Aryl azide#XL1]LEK[#XL1]SEFD`
- ❌ Interchain cross-linkers: `EVTK[X:Aryl azide#XL1]L//EK[#XL1]SEFD`
- ❌ Branches: `ED[MOD:00093#BRANCH]//D[#BRANCH]ATR`

### Ion Notation
`SEQUEN-[b-type-ion]` supports parsing, serialization, mass, and composition calculations. Fragment objects also support mzPAF serialization.
