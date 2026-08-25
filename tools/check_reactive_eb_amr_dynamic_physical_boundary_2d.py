#!/usr/bin/env python3
"""Validate dynamic EB AMR patch planning on an outflow physical boundary."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


ROOT_CELLS = 14
ROOT_DX = 1.0 / ROOT_CELLS
FINE_DX = 0.5 * ROOT_DX
FINAL_TIME = 1.0e-8


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError(f"{path.name}: empty output")
    return rows


def check_physical(rows: list[dict[str, str]]) -> set[int]:
    species = [name for name in rows[0] if name.startswith("Y_")]
    classes: set[int] = set()
    for row in rows:
        values = [float(value) for value in row.values()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("nonfinite dynamic-boundary EB AMR output")
        if abs(float(row["time"]) - FINAL_TIME) > 2.0e-20:
            raise AssertionError("dynamic-boundary final-time mismatch")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError("nonpositive dynamic-boundary state")
        if not 0.0 <= float(row["volume_fraction"]) <= 1.0:
            raise AssertionError("invalid dynamic-boundary volume fraction")
        if abs(sum(float(row[name]) for name in species) - 1.0) > 8.0e-13:
            raise AssertionError("dynamic-boundary species closure drift")
        classes.add(int(row["cell_type"]))
    return classes


def unique_axis(rows: list[dict[str, str]], name: str) -> list[float]:
    return sorted({float(row[name]) for row in rows})


def infer_bounds(axis: list[float]) -> tuple[int, int]:
    for left, right in zip(axis, axis[1:]):
        if abs((right - left) - FINE_DX) > 3.0e-15:
            raise AssertionError("dynamic-boundary fine spacing mismatch")
    lower = round(axis[0] / ROOT_DX - 0.25) + 1
    upper = round(axis[-1] / ROOT_DX + 0.25)
    expected_lower = (lower - 1) * ROOT_DX + 0.5 * FINE_DX
    expected_upper = upper * ROOT_DX - 0.5 * FINE_DX
    if abs(axis[0] - expected_lower) > 3.0e-15 or abs(
        axis[-1] - expected_upper
    ) > 3.0e-15:
        raise AssertionError("dynamic-boundary fine alignment mismatch")
    if len(axis) != 2 * (upper - lower + 1):
        raise AssertionError("dynamic-boundary refinement-ratio mismatch")
    return lower, upper


def check_patch(
    rows: list[dict[str, str]], center_x: float, center_y: float
) -> tuple[int, int, int, int, set[int]]:
    xs = unique_axis(rows, "x")
    ys = unique_axis(rows, "y")
    if len(rows) != len(xs) * len(ys):
        raise AssertionError("dynamic-boundary fine output is not rectangular")
    i_lower, i_upper = infer_bounds(xs)
    j_lower, j_upper = infer_bounds(ys)
    if not (xs[0] < center_x < xs[-1] and ys[0] < center_y < ys[-1]):
        raise AssertionError("dynamic-boundary patch missed its hotspot")
    return i_lower, i_upper, j_lower, j_upper, check_physical(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine-one", required=True, type=Path)
    parser.add_argument("--fine-two", required=True, type=Path)
    args = parser.parse_args()

    coarse = read_rows(args.coarse)
    if len(coarse) != ROOT_CELLS * ROOT_CELLS:
        raise AssertionError("incorrect dynamic-boundary root dimensions")
    coarse_classes = check_physical(coarse)
    first_rows = read_rows(args.fine_one)
    second_rows = read_rows(args.fine_two)
    first = check_patch(first_rows, 0.03, 0.35)
    second = check_patch(second_rows, 0.75, 0.75)
    if first[0] != 1:
        raise AssertionError("first dynamic patch missed the x-lower boundary")
    separated = (
        first[1] + 2 < second[0]
        or second[1] + 2 < first[0]
        or first[3] + 2 < second[2]
        or second[3] + 2 < first[2]
    )
    if not separated:
        raise AssertionError("dynamic-boundary patches violate separation")
    if coarse_classes != {0, 1, 2}:
        raise AssertionError(f"incomplete root EB classes: {coarse_classes}")
    if set.union(first[4], second[4]) != {0, 1, 2}:
        raise AssertionError("dynamic fine patches do not span all EB classes")
    if max(float(row["temperature"]) for row in first_rows) <= 1200.0:
        raise AssertionError("physical-side patch lost its hotspot")
    if max(float(row["temperature"]) for row in second_rows) <= 1200.0:
        raise AssertionError("interior patch lost its hotspot")
    print("check_reactive_eb_amr_dynamic_physical_boundary_2d: PASS")


if __name__ == "__main__":
    main()
