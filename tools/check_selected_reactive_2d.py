#!/usr/bin/env python3
"""Validate dynamic-schema output from a selected regular 2D application."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


BASE_COLUMNS = (
    "time",
    "x",
    "y",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
)


def relative_span(values: list[float]) -> float:
    scale = max(1.0, max(abs(value) for value in values))
    return (max(values) - min(values)) / scale


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--x-lower", type=float, default=0.0)
    parser.add_argument("--x-upper", type=float, default=0.004)
    parser.add_argument("--y-lower", type=float, default=0.0)
    parser.add_argument("--y-upper", type=float, default=0.004)
    parser.add_argument("--require-uniform", action="store_true")
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    args = parser.parse_args()

    species = tuple(args.species)
    columns = (*BASE_COLUMNS, *(f"Y_{name}" for name in species))
    with args.input.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != columns:
            raise AssertionError(f"unexpected columns: {reader.fieldnames}")
        rows = [{name: float(row[name]) for name in columns} for row in reader]

    if len(rows) != args.nx * args.ny:
        raise AssertionError("selected reactive 2D row count is incorrect")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("selected reactive 2D output contains nonfinite data")

    dx = (args.x_upper - args.x_lower) / args.nx
    dy = (args.y_upper - args.y_lower) / args.ny
    maximum_closure_error = 0.0
    for index, row in enumerate(rows):
        i = index % args.nx
        j = index // args.nx
        expected_x = args.x_lower + (i + 0.5) * dx
        expected_y = args.y_lower + (j + 0.5) * dy
        if max(abs(row["x"] - expected_x), abs(row["y"] - expected_y)) > 5.0e-15:
            raise AssertionError(f"row {index} is not in x-fastest order")
        if abs(row["time"] - args.final_time) > 5.0e-15:
            raise AssertionError(f"row {index} has an unexpected time")
        if min(row["rho"], row["pressure"], row["temperature"], row["rhoE"]) <= 0.0:
            raise AssertionError("selected reactive 2D state is nonphysical")
        fractions = [row[f"Y_{name}"] for name in species]
        if min(fractions) < -2.0e-12:
            raise AssertionError("selected reactive 2D species is negative")
        maximum_closure_error = max(
            maximum_closure_error, abs(sum(fractions) - 1.0)
        )
    if maximum_closure_error > 2.0e-10:
        raise AssertionError(
            f"selected reactive 2D closure error {maximum_closure_error}"
        )

    state_columns = tuple(name for name in columns if name not in {"time", "x", "y"})
    maximum_uniform_error = max(
        relative_span([row[name] for row in rows]) for name in state_columns
    )
    if args.require_uniform and maximum_uniform_error > 2.0e-13:
        raise AssertionError(
            f"selected reactive 2D uniformity error {maximum_uniform_error}"
        )

    maximum_requested_species_change = 0.0
    activity_options = (args.activity_species, args.initial_mass_fraction)
    if any(value is not None for value in activity_options):
        if any(value is None for value in activity_options):
            raise AssertionError("both activity options are required")
        column = f"Y_{args.activity_species}"
        if column not in columns:
            raise AssertionError("activity species is not in the output schema")
        maximum_requested_species_change = max(
            abs(row[column] - args.initial_mass_fraction) for row in rows
        )
        if maximum_requested_species_change < args.minimum_change:
            raise AssertionError("selected reactive 2D chemistry was not active")

    print(f"minimum_density={min(row['rho'] for row in rows):.16e}")
    print(f"minimum_pressure={min(row['pressure'] for row in rows):.16e}")
    print(f"minimum_temperature={min(row['temperature'] for row in rows):.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    print(f"maximum_uniform_error={maximum_uniform_error:.16e}")
    print(
        "maximum_requested_species_change="
        f"{maximum_requested_species_change:.16e}"
    )
    print("selected reactive 2D regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
