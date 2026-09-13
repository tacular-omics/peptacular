"""Deterministic mass timings, including parsing and reused annotations.

Run with ``python profile/benchmark_mass.py`` from an installed checkout.
"""

import json
import platform
import random
import statistics
import timeit

import peptacular as pt


def main():
    rng = random.Random(42)
    peptides = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=20)) for _ in range(200)]
    modified = ["[Acetyl]-M[Oxidation]" + seq for seq in peptides]
    proteins = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=500)) for _ in range(200)]
    results = {}
    for label, sequences in [("peptide", peptides), ("modified", modified), ("protein", proteins)]:
        for representation, inputs in [("string", sequences), ("annotation", [pt.parse(seq) for seq in sequences])]:
            for operation, kwargs in [("mass", {}), ("mass_z2", {"charge": 2}), ("mz_z2", {"charge": 2})]:
                func = pt.mz if operation == "mz_z2" else pt.mass

                def run(func=func, inputs=inputs, kwargs=kwargs):
                    return func(inputs, method="sequential", **kwargs)

                run()
                timings = timeit.repeat(run, number=10, repeat=5)
                results[f"{label}/{representation}/{operation}"] = statistics.median(timings) * 1e6 / (10 * len(inputs))
    print(json.dumps({"python": platform.python_version(), "microseconds_per_item": results}, indent=2))


if __name__ == "__main__":
    main()
