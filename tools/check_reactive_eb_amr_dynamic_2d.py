#!/usr/bin/env python3
"""Validate the runnable temperature-tagged reactive EB AMR case."""

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


def unique_coordinates(rows: list[dict[str, str]], name: str) -> list[float]:
    return sorted({float(row[name]) for row in rows})


def check_physical_rows(rows: list[dict[str, str]]) -> None:
    species_columns = [name for name in rows[0] if name.startswith("Y_")]
    for row in rows:
        values = [float(value) for value in row.values()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("nonfinite dynamic EB AMR output")
        if abs(float(row["time"]) - 1.0e-8) > 2.0e-21:
            raise AssertionError("incorrect dynamic EB AMR final time")
        if not 0.0 <= float(row["volume_fraction"]) <= 1.0:
            raise AssertionError("invalid dynamic EB volume fraction")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError("nonpositive dynamic EB AMR state")
        closure = sum(float(row[name]) for name in species_columns)
        if abs(closure - 1.0) > 5.0e-13:
            raise AssertionError("dynamic EB AMR composition closure drift")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine", required=True, type=Path)
    args = parser.parse_args()

    coarse_rows = read_rows(args.coarse)
    fine_rows = read_rows(args.fine)
    if len(coarse_rows) != 12 * 12:
        raise AssertionError("incorrect dynamic EB AMR coarse dimensions")
    check_physical_rows(coarse_rows)
    check_physical_rows(fine_rows)

    fine_x = unique_coordinates(fine_rows, "x")
    fine_y = unique_coordinates(fine_rows, "y")
    if len(fine_rows) != len(fine_x) * len(fine_y):
        raise AssertionError("dynamic fine output is not rectangular")
    if len(fine_x) % 2 != 0 or len(fine_y) % 2 != 0:
        raise AssertionError("dynamic fine patch is not ratio-two aligned")
    if not (fine_x[0] < 0.0072 < fine_x[-1]):
        raise AssertionError("dynamic fine patch missed hotspot in x")
    if not (fine_y[0] < 0.0062 < fine_y[-1]):
        raise AssertionError("dynamic fine patch missed hotspot in y")
    static_upper = 5.0 * 0.01 / 12.0
    if fine_x[-1] <= static_upper or fine_y[-1] <= static_upper:
        raise AssertionError("dynamic fine patch retained static bounds")
    temperatures = [float(row["temperature"]) for row in fine_rows]
    if max(temperatures) < 1200.0 or min(temperatures) < 900.0:
        raise AssertionError("dynamic fine patch lost hotspot structure")
    print("check_reactive_eb_amr_dynamic_2d: PASS")


if __name__ == "__main__":
    main()
