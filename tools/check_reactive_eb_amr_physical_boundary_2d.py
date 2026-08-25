#!/usr/bin/env python3
"""Validate a stationary EB AMR patch on an outflow physical boundary."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_rows(path: Path, expected_rows: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected_rows:
        raise AssertionError(
            f"{path.name}: expected {expected_rows} rows, got {len(rows)}"
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
) -> set[int]:
    dx = (x_upper - x_lower) / nx
    dy = (y_upper - y_lower) / ny
    species = [name for name in rows[0] if name.startswith("Y_")]
    cell_types: set[int] = set()
    for index, row in enumerate(rows):
        values = [float(value) for value in row.values()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("nonfinite physical-boundary EB AMR output")
        expected_x = x_lower + (index % nx + 0.5) * dx
        expected_y = y_lower + (index // nx + 0.5) * dy
        if abs(float(row["x"]) - expected_x) > 3.0e-16 or abs(
            float(row["y"]) - expected_y
        ) > 3.0e-16:
            raise AssertionError("physical-boundary EB AMR coordinate mismatch")
        if abs(float(row["time"]) - 1.0e-7) > 2.0e-20:
            raise AssertionError("physical-boundary EB AMR final time mismatch")
        if not 0.0 <= float(row["volume_fraction"]) <= 1.0:
            raise AssertionError("invalid physical-boundary volume fraction")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError("nonpositive physical-boundary state")
        if abs(float(row["pressure"]) - 101325.0) > 3.0e-6:
            raise AssertionError("physical-boundary pressure drift")
        if abs(float(row["temperature"]) - 1000.0) > 3.0e-8:
            raise AssertionError("physical-boundary temperature drift")
        if math.hypot(float(row["u"]), float(row["v"])) > 3.0e-8:
            raise AssertionError("physical-boundary velocity drift")
        if abs(sum(float(row[name]) for name in species) - 1.0) > 3.0e-13:
            raise AssertionError("physical-boundary species closure drift")
        cell_types.add(int(row["cell_type"]))
    return cell_types


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine", required=True, type=Path)
    args = parser.parse_args()

    coarse = load_rows(args.coarse, 12 * 12)
    fine = load_rows(args.fine, 12 * 18)
    coarse_types = check_level(coarse, 12, 12, 0.0, 0.01, 0.0, 0.01)
    fine_types = check_level(fine, 12, 18, 0.0, 0.005, 0.01 / 12.0, 0.01 * 10 / 12.0)
    if coarse_types != {0, 1, 2} or fine_types != {0, 1, 2}:
        raise AssertionError(
            f"incomplete EB classes: coarse={coarse_types}, fine={fine_types}"
        )
    fine_x = min(float(row["x"]) for row in fine)
    if abs(fine_x - 0.5 * (0.005 / 12.0)) > 3.0e-16:
        raise AssertionError("fine patch does not begin at the x-lower boundary")
    print("check_reactive_eb_amr_physical_boundary_2d: PASS")


if __name__ == "__main__":
    main()
