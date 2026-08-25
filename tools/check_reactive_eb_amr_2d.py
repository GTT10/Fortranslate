#!/usr/bin/env python3
"""Validate the runnable stationary two-level reactive EB AMR case."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_rows(path: Path, expected_count: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected_count:
        raise AssertionError(
            f"{path.name}: expected {expected_count} cells, found {len(rows)}"
        )
    return rows


def check_level(
    rows: list[dict[str, str]],
    nx: int,
    ny: int,
    x_lower: float,
    x_upper: float,
    y_lower: float,
    y_upper: float,
) -> None:
    species_columns = [name for name in rows[0] if name.startswith("Y_")]
    cell_types: set[int] = set()
    dx = (x_upper - x_lower) / nx
    dy = (y_upper - y_lower) / ny
    for index, row in enumerate(rows):
        values = [float(value) for value in row.values()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("nonfinite EB AMR output")
        if abs(float(row["time"]) - 1.0e-7) > 2.0e-20:
            raise AssertionError("incorrect EB AMR final time")
        expected_x = x_lower + (index % nx + 0.5) * dx
        expected_y = y_lower + (index // nx + 0.5) * dy
        if abs(float(row["x"]) - expected_x) > 3.0e-16 or abs(
            float(row["y"]) - expected_y
        ) > 3.0e-16:
            raise AssertionError("incorrect EB AMR coordinate order")
        cell_types.add(int(row["cell_type"]))
        if not 0.0 <= float(row["volume_fraction"]) <= 1.0:
            raise AssertionError("invalid EB AMR volume fraction")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError("nonpositive EB AMR state")
        if abs(float(row["pressure"]) - 101325.0) > 3.0e-6:
            raise AssertionError("stationary EB AMR pressure drift")
        if abs(float(row["temperature"]) - 1000.0) > 3.0e-8:
            raise AssertionError("stationary EB AMR temperature drift")
        if math.hypot(float(row["u"]), float(row["v"])) > 3.0e-8:
            raise AssertionError("stationary EB AMR velocity drift")
        closure = sum(float(row[name]) for name in species_columns)
        if abs(closure - 1.0) > 3.0e-13:
            raise AssertionError("EB AMR composition closure drift")
    if cell_types != {0, 1, 2}:
        raise AssertionError(f"missing EB cell class: {sorted(cell_types)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine", required=True, type=Path)
    args = parser.parse_args()

    coarse_rows = read_rows(args.coarse, 12 * 12)
    fine_rows = read_rows(args.fine, 18 * 18)
    check_level(coarse_rows, 12, 12, 0.0, 0.01, 0.0, 0.01)
    lower = 0.01 / 12.0
    upper = 10.0 * 0.01 / 12.0
    check_level(fine_rows, 18, 18, lower, upper, lower, upper)
    print("check_reactive_eb_amr_2d: PASS")


if __name__ == "__main__":
    main()
