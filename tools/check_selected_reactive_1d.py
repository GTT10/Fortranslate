#!/usr/bin/env python3
"""Validate a selected-mechanism reactive 1D CSV without fixed chemistry."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


BASE_COLUMNS = (
    "time",
    "x",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
)


def relative_span(values: list[float]) -> float:
    return (max(values) - min(values)) / max(1.0, max(abs(value) for value in values))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--require-uniform", action="store_true")
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected_columns = [*BASE_COLUMNS, *(f"Y_{name}" for name in args.species)]
        if reader.fieldnames != expected_columns:
            raise AssertionError(
                f"selected reactive 1D columns differ: {reader.fieldnames!r}"
            )
        rows = [{name: float(value) for name, value in row.items()} for row in reader]

    if len(rows) != args.nx:
        raise AssertionError(f"expected {args.nx} rows, found {len(rows)}")
    values = [value for row in rows for value in row.values()]
    if not all(math.isfinite(value) for value in values):
        raise AssertionError("selected reactive 1D output contains nonfinite data")
    time_error = max(abs(row["time"] - args.final_time) for row in rows)
    if time_error > 5.0e-15 * max(1.0, abs(args.final_time)):
        raise AssertionError(f"selected reactive 1D final-time error is {time_error}")
    coordinates = [row["x"] for row in rows]
    if coordinates != sorted(coordinates) or len(set(coordinates)) != args.nx:
        raise AssertionError("selected reactive 1D coordinates are not strictly ordered")

    minimum_density = min(row["rho"] for row in rows)
    minimum_pressure = min(row["pressure"] for row in rows)
    minimum_temperature = min(row["temperature"] for row in rows)
    minimum_fraction = min(
        row[f"Y_{name}"] for row in rows for name in args.species
    )
    closure_error = max(
        abs(sum(row[f"Y_{name}"] for name in args.species) - 1.0)
        for row in rows
    )
    if minimum_density <= 0.0 or minimum_pressure <= 0.0 or minimum_temperature <= 0.0:
        raise AssertionError("selected reactive 1D output lost thermodynamic positivity")
    if minimum_fraction < -2.0e-12 or closure_error > 2.0e-10:
        raise AssertionError("selected reactive 1D composition is invalid")

    uniform_error = 0.0
    if args.require_uniform:
        for name in (*BASE_COLUMNS[2:], *(f"Y_{item}" for item in args.species)):
            uniform_error = max(uniform_error, relative_span([row[name] for row in rows]))
        if uniform_error > 2.0e-13:
            raise AssertionError(
                f"selected reactive 1D uniform-state drift is {uniform_error}"
            )

    activity = 0.0
    activity_options = (args.activity_species, args.initial_mass_fraction)
    if any(value is not None for value in activity_options):
        if any(value is None for value in activity_options):
            raise AssertionError("both activity options are required")
        column = f"Y_{args.activity_species}"
        if column not in rows[0]:
            raise AssertionError(f"activity species is absent: {args.activity_species}")
        activity = max(
            abs(row[column] - args.initial_mass_fraction) for row in rows
        )
        if activity < args.minimum_change:
            raise AssertionError(
                f"selected reactive 1D chemistry activity is only {activity}"
            )

    print(f"minimum_density={minimum_density:.16e}")
    print(f"minimum_pressure={minimum_pressure:.16e}")
    print(f"minimum_temperature={minimum_temperature:.16e}")
    print(f"minimum_mass_fraction={minimum_fraction:.16e}")
    print(f"maximum_closure_error={closure_error:.16e}")
    print(f"maximum_uniform_error={uniform_error:.16e}")
    print(f"maximum_requested_species_change={activity:.16e}")
    print("selected reactive 1D regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
