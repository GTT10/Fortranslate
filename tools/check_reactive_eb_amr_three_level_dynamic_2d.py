#!/usr/bin/env python3
"""Check tag-driven finest-patch regridding in public three-level EB AMR."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError(f"{path.name}: empty output")
    species = [name for name in rows[0] if name.startswith("Y_")]
    cell_types: set[int] = set()
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
        if abs(float(row["time"]) - 1.0e-7) > 2.0e-20:
            raise AssertionError(f"{path.name}: incorrect final time")
        cell_types.add(int(row["cell_type"]))
        if int(row["cell_type"]) != 0:
            if float(row["rho"]) <= 0.0 or float(row["temperature"]) <= 0.0:
                raise AssertionError(f"{path.name}: invalid active thermodynamics")
            closure = sum(float(row[name]) for name in species)
            if abs(closure - 1.0) > 8.0e-12:
                raise AssertionError(f"{path.name}: species closure drift")
    if cell_types != {0, 1, 2}:
        raise AssertionError(f"{path.name}: incomplete EB coverage {cell_types}")
    return rows


def dimensions(rows: list[dict[str, str]]) -> tuple[int, int]:
    x_coordinates = {float(row["x"]) for row in rows}
    y_coordinates = {float(row["y"]) for row in rows}
    if len(rows) != len(x_coordinates) * len(y_coordinates):
        raise AssertionError("output is not one complete rectangular level")
    return len(x_coordinates), len(y_coordinates)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--middle", required=True, type=Path)
    parser.add_argument("--finest", required=True, type=Path)
    args = parser.parse_args()

    root = load(args.root)
    middle = load(args.middle)
    finest = load(args.finest)
    if dimensions(root) != (12, 12):
        raise AssertionError("incorrect root dimensions")
    if dimensions(middle) == (20, 20):
        raise AssertionError("middle patch retained the configured seed")
    finest_dimensions = dimensions(finest)
    if min(finest_dimensions) < 8:
        raise AssertionError("finest patch is smaller than the tagged minimum")
    x_coordinates = sorted({float(row["x"]) for row in finest})
    y_coordinates = sorted({float(row["y"]) for row in finest})
    if not (x_coordinates[0] < 0.00437 < x_coordinates[-1]):
        raise AssertionError("finest patch does not cross the embedded boundary")
    print("check_reactive_eb_amr_three_level_dynamic_2d: PASS")


if __name__ == "__main__":
    main()
