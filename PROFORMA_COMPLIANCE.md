# ProForma Parser Implementation Checklist

## 5.3.1 Base-ProForma Compliance

- [x] **Amino acids (+UO)** - `AAHCFKUOT` (§6.1)
- [x] **Unimod names** - `PEM[Oxidation]AT` (§6.2.1)
- [x] **PSI-MOD names** - `PEM[monohydroxylated residue]AT` (§6.2.1)
- [x] **Unimod numbers** - `PEM[UNIMOD:35]AT` (§6.2.2)
- [x] **PSI-MOD numbers** - `PEM[MOD:00425]AT` (§6.2.2)
- [x] **Delta masses** - `PEM[+15.995]AT` (§6.2.3)
- [x] **N-terminal modifications** - `[Carbamyl]-QPEPTIDE` (§6.3)
- [x] **C-terminal modifications** - `PEPTIDEG-[Methyl]` (§6.3)
- [x] **Labile modifications** - `{Glycan:Hex}EM[U:Oxidation]EV` (§6.4)
- [x] **Multiple modifications** - `MPGNW[Oxidation][Carboxymethyl]PESQE` (§6.5)
- [x] **Information tag** - `ELV[INFO:AnyString]IS` (§6.6)

## 5.3.2 Level 2-ProForma Compliance

- [x] **Ambiguous amino acids** - `BZJX` (§7.1)
- [x] **Prefixed delta masses** - `PEM[U:+15.995]AT` (§7.2)
- [x] **Mass gap** - `PEX[+147.035]AT` (§7.3)
- [x] **Formulas** - `PEM[Formula:O]AT`, `PEM[Formula:[17O1]]AT` (§7.4)
- [x] **Mass with interpretation** - `PEM[+15.995|Oxidation]AT` (§7.5)
- [x] **Unknown mod position** - `[Oxidation]?PEMAT` (§7.6.1)
- [x] **Set of positions** - `PEP[Oxidation#1]M[#1]AT` (§7.6.2)
- [x] **Range of positions** - `PRT(ESFRMS)[+19.0523]ISK` (§7.6.3)
- [x] **Position scores** - `PEP[Oxidation#1(0.95)]M[#1(0.05)]AT` (§7.6.4)
- [x] **Range position scores** - `(PEP)[Oxidation#1(0.95)]M[#1(0.05)]AT` (§7.6.5)
- [x] **Amino acid ambiguity** - `(?VCH)AT` (§7.7)
- [x] **Modification prefixes** - `PEPM[U:Oxidation]AS[M:O-phospho-L-serine]` (§7.8)

## 5.3.3 Level 2-ProForma + Top-Down

- [x] **RESID modifications** - `EM[R:L-methionine sulfone]EM[RESID:AA0581]` (§8.1)
- [x] **Names** - `(>Heavy chain)EVQLVESG` (§8.2)

## 5.3.4 Level 2-ProForma + Cross-Linking

- [x] **XL-MOD modifications** - `EVTK[X:Aryl azide]LEK[XLMOD:00114]SEFD` (§9.1)
- [ ] **Cross-linkers (intrachain)** - `EVTK[X:Aryl azide#XL1]LEK[#XL1]SEFD` (§9.2.1)
- [ ] **Cross-linkers (interchain)** - `EVTK[X:Aryl azide#XL1]L//EK[#XL1]SEFD` (§9.2.2)
- [ ] **Branches** - `ED[MOD:00093#BRANCH]//D[#BRANCH]ATR` (§9.3) 

## 5.3.5 Level 2-ProForma + Glycans

- [x] **GNO modifications** - `NEEYN[GNO:G59626AS]K` (§10.1)
- [x] **Glycan compositions** - `NEEYN[Glycan:Hex5HexNAc4NeuAc1]K` (§10.2)

## 5.3.6 Level 2-ProForma + Advanced Complexity

- [x] **Charged formulas** - `SEQUEN[Formula:Zn1:z+2]CE` (§11.1)
- [x] **Controlling placement** - `PTI(MERMERME)[+32|Position:E]PTIDE` (§11.2) 
- [x] **Global isotope** - `<13C>CARBON` (§11.3.1)
- [x] **Fixed modifications** - `<[Oxidation]@M>ATPEMILTCMGCLK` (§11.3.2)
- [x] **Chimeric spectra** - `NEEYN+SEQUEN` (§11.4)
- [x] **Charges** - `SEQUEN/2`, `SEQUEN/[Na:z+1,H:z+1]` (§11.5)
- [x] **Ion notation** - `SEQUEN-[b-type-ion]` (§11.6)