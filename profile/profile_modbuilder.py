import cProfile
import pstats
from io import StringIO

import peptacular as pt


def profile_single_peptide_mods():
    """Profile building mods for a single peptide."""
    peptide = "PEPTIDEWITHMETHIONINEPEPTIDEWITHMETHIONINEPEPTIDEWITHMETHIONINEPEPTIDEWITHMETHIONINEPEPTIDEWITHMETHIONINEPEPTIDEWITHMETHIONINEPEPTIDEWITHMETHIONINEPEPTIDEWITHMETHIONINE"

    print(f"Profiling modifications for single peptide: {peptide}")

    profiler = cProfile.Profile()
    profiler.enable()

    modified = list(
        pt.modify(
            peptide,
            internal_static={"C": [57.021]},
            internal_variable={"M": [15.995], "E": [79.966]},
            max_variable_mods=2,
        )
    )

    profiler.disable()

    print(f"Generated {len(modified)} variants")

    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
    ps.print_stats(40)
    print(s.getvalue())


if __name__ == "__main__":
    print("=" * 60)
    print("SINGLE PEPTIDE MOD BUILDING")
    print("=" * 60)
    profile_single_peptide_mods()
