#!/usr/bin/env python3
"""Validate the public two-level reactive EB multipatch application."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError(f"{path.name}: empty output")
    return rows


def coordinates(rows: list[dict[str, str]], name: str) -> list[float]:
    return sorted({float(row[name]) for row in rows})


def check_physical(rows: list[dict[str, str]]) -> set[int]:
    species_columns = [name for name in rows[0] if name.startswith("Y_")]
    cell_types: set[int] = set()
    for row in rows:
        values = [float(value) for value in row.values()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("nonfinite multipatch EB AMR output")
        if abs(float(row["time"]) - 1.0e-8) > 2.0e-20:
            raise AssertionError("incorrect multipatch final time")
        if not 0.0 <= float(row["volume_fraction"]) <= 1.0:
            raise AssertionError("invalid multipatch volume fraction")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError("nonpositive multipatch state")
        closure = sum(float(row[name]) for name in species_columns)
        if abs(closure - 1.0) > 8.0e-13:
            raise AssertionError("multipatch composition closure drift")
        cell_types.add(int(row["cell_type"]))
    return cell_types


def check_patch(
    rows: list[dict[str, str]], center_x: float, center_y: float
) -> tuple[list[float], list[float]]:
    xs = coordinates(rows, "x")
    ys = coordinates(rows, "y")
    if len(rows) != len(xs) * len(ys):
        raise AssertionError("fine patch output is not rectangular")
    if len(xs) != 10 or len(ys) != 10:
        raise AssertionError("fine patch is not the expected ratio-two 5x5 region")
    if not xs[0] < center_x < xs[-1] or not ys[0] < center_y < ys[-1]:
        raise AssertionError("fine patch missed its hotspot")
    check_physical(rows)
    return xs, ys


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine-one", required=True, type=Path)
    parser.add_argument("--fine-two", required=True, type=Path)
    args = parser.parse_args()

    coarse_rows = read_rows(args.coarse)
    if len(coarse_rows) != 14 * 14:
        raise AssertionError("incorrect multipatch coarse dimensions")
    coarse_types = check_physical(coarse_rows)
    first_rows = read_rows(args.fine_one)
    second_rows = read_rows(args.fine_two)
    first_x, first_y = check_patch(first_rows, 0.25, 0.55)
    second_x, second_y = check_patch(second_rows, 0.75, 0.75)
    if max(first_x) >= min(second_x) or max(first_y) >= min(second_y):
        raise AssertionError("fine patch outputs are not separated and ordered")
    if coarse_types != {0, 1, 2}:
        raise AssertionError(f"missing coarse EB cell class: {sorted(coarse_types)}")
    if max(float(row["temperature"]) for row in first_rows + second_rows) <= 1200.0:
        raise AssertionError("multipatch outputs lost the hotspot temperature")
    print("check_reactive_eb_amr_multipatch_2d: PASS")


if __name__ == "__main__":
    main()
