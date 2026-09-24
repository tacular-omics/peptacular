from collections import Counter

import tacular
from fastatacular import read_fasta

import peptacular as pt

if __name__ == "__main__":
    FASTA = "/home/patrick-garrett/Data/fasta/human_and_contaminants.fasta"
    fasta_entries: list[str] = [f.sequence for f in read_fasta(FASTA)]

    # remove any peptides with X, B, J, U, or Z
    fasta_entries = [seq for seq in fasta_entries if not any(aa in seq for aa in "XBJUZ")]

    masses: list[int | float] = pt.mass(fasta_entries, monoisotopic=False, ion_type='n')
    compositions: list[Counter[tacular.ElementInfo]] = pt.comp(fasta_entries)

    combined_masses = sum(masses)
    combined_composition = Counter[tacular.ElementInfo]()
    for comp in compositions:
        combined_composition += comp

    print(f"Average mass of all sequences: {combined_masses / len(masses)} Da")
    print("Average composition of all sequences:")
    for element, count in combined_composition.items():
        print(f"  {element.symbol}: {count}")

    averagine_ratios_per_dalton: dict[str, float] = {}

    # divide each element count by total mass to get ratio per dalton
    for element, count in combined_composition.items():
        averagine_ratios_per_dalton[str(element.symbol)] = count / combined_masses

    # roudn to 5 decimal places
    for element in averagine_ratios_per_dalton:
        averagine_ratios_per_dalton[element] = round(averagine_ratios_per_dalton[element], 7)

    print("Averagine composition ratios per dalton:")
    for element, ratio in averagine_ratios_per_dalton.items():
        print(f"  {element}: {ratio} atoms/Da")

    print(averagine_ratios_per_dalton)

