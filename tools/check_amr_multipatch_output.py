#!/usr/bin/env python3
"""Check ordered, conservative-coverage separated multipatch AMR output."""
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
    parser.add_argument("--minimum-patches", type=int, default=2)
    args = parser.parse_args()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise AssertionError("multipatch AMR output is empty")
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("multipatch AMR output contains non-finite data")
    levels = {int(row["level"]) for row in rows}
    if levels != {0, 1}:
        raise AssertionError(f"expected parent and fine levels, found {levels}")
    spacings = sorted({row["cell_dx"] for row in rows}, reverse=True)
    if len(spacings) != 2 or abs(spacings[0] / spacings[1] - 2.0) > 1.0e-12:
        raise AssertionError(f"unexpected multipatch spacings {spacings}")
    fine_positions = [row["x"] for row in rows if int(row["level"]) == 1]
    patch_count = 1 + sum(
        right - left > 1.5 * spacings[1]
        for left, right in zip(fine_positions, fine_positions[1:])
    )
    if patch_count < args.minimum_patches:
        raise AssertionError(f"expected separated fine patches, found {patch_count}")
    coverage = math.fsum(row["cell_dx"] for row in rows)
    if abs(coverage - args.domain_length) > 3.0e-13:
        raise AssertionError(f"multipatch composite coverage is {coverage}")
    positions = [row["x"] for row in rows]
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise AssertionError("multipatch composite output is not ordered")
    if min(row["rho"] for row in rows) <= 0.0:
        raise AssertionError("multipatch density is not positive")
    if min(row["pressure"] for row in rows) <= 0.0:
        raise AssertionError("multipatch pressure is not positive")
    if min(row["temperature"] for row in rows) <= 0.0:
        raise AssertionError("multipatch temperature is not positive")
    closure = max(
        abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0)
        for row in rows
    )
    if closure > 3.0e-10:
        raise AssertionError(f"multipatch species closure error is {closure}")
    print(f"rows={len(rows)}")
    print(f"patches={patch_count}")
    print(f"coverage={coverage:.16e}")
    print(f"maximum_closure_error={closure:.16e}")
    print("Dynamic multipatch AMR output: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
