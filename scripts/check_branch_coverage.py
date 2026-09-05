"""Enforce branch coverage separately from coverage.py's combined percentage."""

import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("report", nargs="?", type=Path, default=Path("coverage.json"))
parser.add_argument("--minimum", type=float, default=79.0)
args = parser.parse_args()
totals = json.loads(args.report.read_text())["totals"]
branches = totals["num_branches"]
if branches <= 0:
    raise SystemExit("No branch data found. Run pytest with --cov-branch.")
percentage = 100 * totals["covered_branches"] / branches
print(f"Branch coverage: {percentage:.2f}% (minimum {args.minimum:.2f}%)")
if percentage < args.minimum:
    raise SystemExit(1)
