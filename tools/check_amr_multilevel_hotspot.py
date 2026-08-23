#!/usr/bin/env python3
"""Check ordered, complete, positive three-level reactive AMR output."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--domain-length", type=float, default=0.012)
    args = parser.parse_args()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise AssertionError("multilevel AMR output is empty")
    levels = {int(row["level"]) for row in rows}
    if levels != {0, 1, 2}:
        raise AssertionError(f"expected three active levels, found {levels}")
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("multilevel AMR output contains non-finite data")
    spacings = sorted({row["cell_dx"] for row in rows}, reverse=True)
    if len(spacings) != 3:
        raise AssertionError(f"expected three cell spacings, found {spacings}")
    if any(abs(coarse / fine - 2.0) > 1.0e-12 for coarse, fine in zip(spacings, spacings[1:])):
        raise AssertionError(f"unexpected refinement ratios for {spacings}")
    coverage = math.fsum(row["cell_dx"] for row in rows)
    if abs(coverage - args.domain_length) > 3.0e-13:
        raise AssertionError(f"composite coverage is {coverage}")
    positions = [row["x"] for row in rows]
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise AssertionError("multilevel composite output is not ordered")
    if min(row["rho"] for row in rows) <= 0.0:
        raise AssertionError("multilevel density is not positive")
    if min(row["pressure"] for row in rows) <= 0.0:
        raise AssertionError("multilevel pressure is not positive")
    if min(row["temperature"] for row in rows) <= 0.0:
        raise AssertionError("multilevel temperature is not positive")
    closure = max(
        abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0)
        for row in rows
    )
    if closure > 3.0e-10:
        raise AssertionError(f"species closure error is {closure}")
    print(f"rows={len(rows)}")
    print(f"coverage={coverage:.16e}")
    print(f"maximum_closure_error={closure:.16e}")
    print("Multilevel AMR hotspot output: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
