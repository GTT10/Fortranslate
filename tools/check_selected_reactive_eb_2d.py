#!/usr/bin/env python3
"""Validate a dynamic-species selected reactive EB 2D CSV result."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


BASE_COLUMNS = [
    "time",
    "x",
    "y",
    "volume_fraction",
    "cell_type",
    "boundary_length",
    "boundary_normal_x",
    "boundary_normal_y",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--nx", required=True, type=int)
    parser.add_argument("--ny", required=True, type=int)
    parser.add_argument("--final-time", required=True, type=float)
    parser.add_argument("--expected-regular", required=True, type=int)
    parser.add_argument("--expected-cut", required=True, type=int)
    parser.add_argument("--expected-covered", required=True, type=int)
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(message)


def main() -> None:
    arguments = parse_arguments()
    expected_header = BASE_COLUMNS + [f"Y_{name}" for name in arguments.species]
    with arguments.input.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_header:
            fail(
                "selected reactive EB 2D CSV schema mismatch: "
                f"{reader.fieldnames!r}"
            )
        rows = list(reader)

    if arguments.activity_species:
        if arguments.activity_species not in arguments.species:
            fail("--activity-species must be present in --species")
        if arguments.initial_mass_fraction is None:
            fail("--activity-species requires --initial-mass-fraction")

    if len(rows) != arguments.nx * arguments.ny:
        fail("selected reactive EB 2D CSV row count mismatch")
    if len({float(row["x"]) for row in rows}) != arguments.nx:
        fail("selected reactive EB 2D x-coordinate count mismatch")
    if len({float(row["y"]) for row in rows}) != arguments.ny:
        fail("selected reactive EB 2D y-coordinate count mismatch")

    type_counts = {0: 0, 1: 0, 2: 0}
    maximum_closure_error = 0.0
    minimum_density = math.inf
    minimum_pressure = math.inf
    minimum_temperature = math.inf
    maximum_activity_change = 0.0
    for row in rows:
        try:
            values = [float(row[column]) for column in expected_header]
        except ValueError as error:
            fail(f"selected reactive EB 2D CSV contains nonnumeric data: {error}")
        if not all(math.isfinite(value) for value in values):
            fail("selected reactive EB 2D CSV contains nonfinite data")
        if not math.isclose(
            float(row["time"]), arguments.final_time, rel_tol=0.0, abs_tol=1.0e-14
        ):
            fail("selected reactive EB 2D final time mismatch")

        cell_type_value = float(row["cell_type"])
        cell_type = int(cell_type_value)
        if cell_type_value != cell_type or cell_type not in type_counts:
            fail("selected reactive EB 2D cell type is invalid")
        type_counts[cell_type] += 1
        volume_fraction = float(row["volume_fraction"])
        if volume_fraction < 0.0 or volume_fraction > 1.0:
            fail("selected reactive EB 2D volume fraction is outside [0,1]")
        if cell_type == 0 and volume_fraction != 0.0:
            fail("selected reactive EB 2D covered cell has nonzero volume")
        if cell_type == 1 and not 0.0 < volume_fraction < 1.0:
            fail("selected reactive EB 2D cut cell lacks a partial volume")
        if cell_type == 2 and volume_fraction != 1.0:
            fail("selected reactive EB 2D regular cell lacks unit volume")
        if float(row["boundary_length"]) < 0.0:
            fail("selected reactive EB 2D boundary length is negative")

        density = float(row["rho"])
        pressure = float(row["pressure"])
        temperature = float(row["temperature"])
        if density <= 0.0 or pressure <= 0.0 or temperature <= 0.0:
            fail("selected reactive EB 2D state is not positive")
        minimum_density = min(minimum_density, density)
        minimum_pressure = min(minimum_pressure, pressure)
        minimum_temperature = min(minimum_temperature, temperature)
        closure = sum(float(row[f"Y_{name}"]) for name in arguments.species)
        maximum_closure_error = max(maximum_closure_error, abs(closure - 1.0))
        if any(
            float(row[f"Y_{name}"]) < -1.0e-13 for name in arguments.species
        ):
            fail("selected reactive EB 2D mass fraction is negative")
        if arguments.activity_species and cell_type != 0:
            activity = float(row[f"Y_{arguments.activity_species}"])
            maximum_activity_change = max(
                maximum_activity_change,
                abs(activity - arguments.initial_mass_fraction),
            )

    expected_counts = {
        2: arguments.expected_regular,
        1: arguments.expected_cut,
        0: arguments.expected_covered,
    }
    if type_counts != expected_counts:
        fail(
            "selected reactive EB 2D cell counts mismatch: "
            f"observed={type_counts}, expected={expected_counts}"
        )
    if maximum_closure_error > 5.0e-12:
        fail("selected reactive EB 2D composition closure is too large")
    if arguments.activity_species:
        if maximum_activity_change < arguments.minimum_change:
            fail("selected reactive EB 2D activity is below the requested gate")

    print(f"cells=regular:{type_counts[2]},cut:{type_counts[1]},covered:{type_counts[0]}")
    print(f"minimum_density={minimum_density:.16e}")
    print(f"minimum_pressure={minimum_pressure:.16e}")
    print(f"minimum_temperature={minimum_temperature:.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    print(f"maximum_requested_species_change={maximum_activity_change:.16e}")
    print("selected reactive EB 2D regression: PASS")


if __name__ == "__main__":
    main()
