#!/usr/bin/env python3
"""Check parent-first regridding in public fixed three-level EB AMR."""

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
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
        if abs(float(row["time"]) - 1.0e-7) > 2.0e-20:
            raise AssertionError(f"{path.name}: incorrect final time")
        if int(row["cell_type"]) != 0:
            if float(row["rho"]) <= 0.0 or float(row["temperature"]) <= 0.0:
                raise AssertionError(f"{path.name}: invalid active thermodynamics")
            closure = sum(float(row[name]) for name in species)
            if abs(closure - 1.0) > 8.0e-12:
                raise AssertionError(f"{path.name}: species closure drift")
    return rows


def dimensions(rows: list[dict[str, str]]) -> tuple[int, int]:
    x_coordinates = {float(row["x"]) for row in rows}
    y_coordinates = {float(row["y"]) for row in rows}
    if len(rows) != len(x_coordinates) * len(y_coordinates):
        raise AssertionError("output is not one complete rectangular level")
    return len(x_coordinates), len(y_coordinates)


def axis_extent(
    rows: list[dict[str, str]], coordinate: str
) -> tuple[float, float, float]:
    coordinates = sorted({float(row[coordinate]) for row in rows})
    spacings = [b - a for a, b in zip(coordinates, coordinates[1:])]
    if not spacings:
        raise AssertionError(f"{coordinate}: level has fewer than two cells")
    spacing = sum(spacings) / len(spacings)
    tolerance = 2.0e-12 * max(1.0, abs(spacing))
    if any(abs(value - spacing) > tolerance for value in spacings):
        raise AssertionError(f"{coordinate}: nonuniform level spacing")
    return (
        coordinates[0] - 0.5 * spacing,
        coordinates[-1] + 0.5 * spacing,
        spacing,
    )


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
    root_cell_types = {int(row["cell_type"]) for row in root}
    if root_cell_types != {0, 1, 2}:
        raise AssertionError(f"incomplete root EB coverage {root_cell_types}")
    if dimensions(middle) == (20, 20):
        raise AssertionError("middle patch retained the configured seed")
    middle_cell_types = {int(row["cell_type"]) for row in middle}
    if 1 not in middle_cell_types:
        raise AssertionError(f"middle: missing cut-cell coverage {middle_cell_types}")
    finest_cell_types = {int(row["cell_type"]) for row in finest}
    if finest_cell_types != {2}:
        raise AssertionError(f"finest: expected regular cells {finest_cell_types}")
    finest_dimensions = dimensions(finest)
    if min(finest_dimensions) < 8:
        raise AssertionError("finest patch is smaller than the tagged minimum")
    for coordinate in ("x", "y"):
        middle_lower, middle_upper, middle_spacing = axis_extent(
            middle, coordinate
        )
        finest_lower, finest_upper, finest_spacing = axis_extent(
            finest, coordinate
        )
        tolerance = 2.0e-12 * max(1.0, abs(middle_spacing))
        if abs(2.0 * finest_spacing - middle_spacing) > tolerance:
            raise AssertionError(f"{coordinate}: incorrect refinement ratio")
        if (
            finest_lower < middle_lower + 2.0 * middle_spacing - tolerance
            or finest_upper > middle_upper - 2.0 * middle_spacing + tolerance
        ):
            raise AssertionError(
                f"{coordinate}: finest patch violates the two-cell middle margin"
            )
    print("check_reactive_eb_amr_three_level_dynamic_2d: PASS")


if __name__ == "__main__":
    main()
