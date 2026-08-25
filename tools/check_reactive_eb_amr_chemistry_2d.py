#!/usr/bin/env python3
"""Check uniform-reactor parity on both active EB AMR levels."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_rows(path: Path, expected: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected:
        raise AssertionError(f"{path.name}: expected {expected} rows, got {len(rows)}")
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite output")
    return rows


def close(actual: float, expected: float, tolerance: float = 8.0e-10) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def check_level(
    rows: list[dict[str, str]], reference: dict[str, str], label: str
) -> None:
    species = [name for name in reference if name.startswith("Y_")]
    state = [
        "rho",
        "u",
        "v",
        "w",
        "pressure",
        "temperature",
        "rhoE",
        *species,
    ]
    cell_types: set[int] = set()
    active = 0
    for row in rows:
        if abs(float(row["time"]) - 2.0e-7) > 4.0e-20:
            raise AssertionError(f"{label}: incorrect final time")
        cell_type = int(row["cell_type"])
        cell_types.add(cell_type)
        if cell_type == 0:
            continue
        active += 1
        for name in state:
            if not close(float(row[name]), float(reference[name])):
                raise AssertionError(
                    f"{label}: active {name} mismatch: "
                    f"{row[name]} != {reference[name]}"
                )
        closure = sum(float(row[name]) for name in species)
        if abs(closure - 1.0) > 5.0e-13:
            raise AssertionError(f"{label}: species closure drift")
    if active == 0 or cell_types != {0, 1, 2}:
        raise AssertionError(f"{label}: incomplete EB cell coverage {cell_types}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine", required=True, type=Path)
    args = parser.parse_args()

    reference_rows = load_rows(args.reference, 8 * 8)
    coarse_rows = load_rows(args.coarse, 8 * 8)
    fine_rows = load_rows(args.fine, 10 * 10)
    reference = reference_rows[0]
    if abs(float(reference["temperature"]) - 1200.0) <= 1.0e-12:
        raise AssertionError("reference chemistry did not advance")
    check_level(coarse_rows, reference, "coarse")
    check_level(fine_rows, reference, "fine")
    print("check_reactive_eb_amr_chemistry_2d: PASS")


if __name__ == "__main__":
    main()
