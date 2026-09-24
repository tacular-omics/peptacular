import time

if __name__ == "__main__":
    import peptacular as pt

    start_time = time.time()
    fasta_seqs: list[pt.FastaSequence] = pt.parse_fasta("/home/patrick-garrett/Data/fasta/human_and_contaminants.fasta")[:100]
    elapsed = time.time() - start_time
    print(f"Parsed {len(fasta_seqs):,} FASTA sequences in {elapsed:.3f} seconds")

    start_time = time.time()
    proteins = [fasta_seq.sequence for fasta_seq in fasta_seqs]
    peptides = set()
    for i, peptide_spans in enumerate(
        pt.digest(
            proteins,
            enzyme=pt.Protease.TRYPSIN,
            missed_cleavages=2,
            semi=False,
            min_len=6,
            max_len=50,
        )
    ):
        for pep_seq, span in peptide_spans:
            peptides.add(pep_seq)

    peptides = list(peptides)
    elapsed = time.time() - start_time
    print(f"Generated {len(peptides):,} unique peptides in {elapsed:.3f} seconds")

    # calculate the mass
    start_time = time.time()
    masses = pt.mass(peptides)
    elapsed = time.time() - start_time

    print(f"Calculated masses for {len(peptides):,} peptides in {elapsed:.3f} seconds")

    start_time = time.time()
    for frags in pt.fragment(
        peptides,
        ion_types=[pt.IonType.B, pt.IonType.Y],
        charges=[1, 2],
        monoisotopic=True,
        isotopes=[1],
    ):
        pass

    elapsed = time.time() - start_time
    print(f"Fragmented {len(peptides):,} peptides in {elapsed:.3f} seconds")
