"""Expand ambiguous modification sites and find ions that distinguish them."""

import peptacular as pt


def run():
    # The range names the candidate residues. No residue specificity is inferred.
    isomers = pt.localization_isomers("PEP(ST)[Phospho]IDE", max_isomers=100)
    assert [a.serialize() for a in isomers] == [
        "PEPS[Phospho]TIDE",
        "PEPST[Phospho]IDE",
    ]
    print("Localization isomers:")
    for isomer in isomers:
        print(isomer.serialize())

    # Choose eligible residues explicitly when adding a modification.
    candidates = pt.candidate_sites("PEPSTIDE", "Phospho", residues="ST")
    assert [position for position, _ in candidates] == [3, 4]
    print("Candidate sites (zero-based):")
    for position, isomer in candidates:
        print(position, isomer.serialize())

    # These fragments are unique to one isomer at the requested m/z tolerance.
    ions = pt.site_determining_ions(
        isomers,
        ion_types=("b", "y"),
        charges=(1, 2),
        tolerance=10,
        tolerance_unit="ppm",
    )
    assert len(ions) == 2 and all(ions)
    for isomer, fragments in zip(isomers, ions, strict=True):
        print(isomer.serialize())
        for fragment in fragments:
            print(f"  {fragment.ion_type}{fragment.position}+{fragment.charge_state}: {fragment.mz:.3f} m/z")

    # With three sites, a middle isomer may share all its ions with neighbours.
    # Pairwise comparisons show which ions distinguish each ordered pair.
    three_sites = pt.localization_isomers("PEP(STY)[Phospho]IDEK")
    pairs = pt.pairwise_site_determining_ions(three_sites)
    assert len(pairs) == 6
    print("Pairwise site-determining ions:")
    for (left, right), fragments in pairs.items():
        print(left, right, [f"{f.ion_type}{f.position}" for f in fragments])


if __name__ == "__main__":
    run()
