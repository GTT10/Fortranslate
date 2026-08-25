#!/usr/bin/env python3
"""Check public static three-level reactive EB AMR output."""

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


def close(actual: float, expected: float, tolerance: float = 2.0e-8) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def check_level(
    rows: list[dict[str, str]], reference: dict[str, str], label: str
) -> None:
    species = [name for name in reference if name.startswith("Y_")]
    state = ["rho", "pressure", "temperature", "rhoE", *species]
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
        if abs(float(row["u"])) > 2.0e-7 or abs(float(row["v"])) > 2.0e-7:
            raise AssertionError(f"{label}: stationary velocity drift")
        closure = sum(float(row[name]) for name in species)
        if abs(closure - 1.0) > 8.0e-12:
            raise AssertionError(f"{label}: species closure drift")
    if active == 0 or cell_types != {0, 1, 2}:
        raise AssertionError(f"{label}: incomplete EB cell coverage {cell_types}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--middle", required=True, type=Path)
    parser.add_argument("--finest", required=True, type=Path)
    args = parser.parse_args()

    reference_rows = load_rows(args.reference, 8 * 8)
    root_rows = load_rows(args.root, 8 * 8)
    middle_rows = load_rows(args.middle, 12 * 12)
    finest_rows = load_rows(args.finest, 16 * 16)
    reference = reference_rows[0]
    if abs(float(reference["temperature"]) - 1200.0) <= 1.0e-12:
        raise AssertionError("reference chemistry did not advance")
    check_level(root_rows, reference, "root")
    check_level(middle_rows, reference, "middle")
    check_level(finest_rows, reference, "finest")
    print("check_reactive_eb_amr_three_level_2d: PASS")


if __name__ == "__main__":
    main()
