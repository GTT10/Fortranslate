#!/usr/bin/env python3
"""Validate dynamic-schema output from a selected regular 3D application."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


BASE_COLUMNS = (
    "time",
    "x",
    "y",
    "z",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
    "rhou",
    "rhov",
    "rhow",
)


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(1.0, abs(expected))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--x-lower", type=float, default=0.0)
    parser.add_argument("--x-upper", type=float, default=0.004)
    parser.add_argument("--y-lower", type=float, default=0.0)
    parser.add_argument("--y-upper", type=float, default=0.004)
    parser.add_argument("--z-lower", type=float, default=0.0)
    parser.add_argument("--z-upper", type=float, default=0.004)
    parser.add_argument("--require-uniform", action="store_true")
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    args = parser.parse_args()

    species = tuple(args.species)
    columns = (
        *BASE_COLUMNS,
        *(f"Y_{name}" for name in species),
        *(f"rhoY_{name}" for name in species),
    )
    with args.input.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != columns:
            raise AssertionError(f"unexpected columns: {reader.fieldnames}")
        rows = [{name: float(row[name]) for name in columns} for row in reader]

    extents = (args.nx, args.ny, args.nz)
    bounds = (
        (args.x_lower, args.x_upper),
        (args.y_lower, args.y_upper),
        (args.z_lower, args.z_upper),
    )
    if len(rows) != math.prod(extents):
        raise AssertionError("selected reactive 3D row count is incorrect")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("selected reactive 3D output contains nonfinite data")

    maximum_relationship_error = 0.0
    maximum_closure_error = 0.0
    for index, row in enumerate(rows):
        indices = (
            index % args.nx,
            (index // args.nx) % args.ny,
            index // (args.nx * args.ny),
        )
        expected_coordinates = tuple(
            lower + (cell + 0.5) * (upper - lower) / extent
            for cell, extent, (lower, upper) in zip(indices, extents, bounds)
        )
        if max(
            abs(row[name] - expected)
            for name, expected in zip(("x", "y", "z"), expected_coordinates)
        ) > 5.0e-15:
            raise AssertionError(f"row {index} is not in x-fastest order")
        if abs(row["time"] - args.final_time) > 5.0e-18:
            raise AssertionError(f"row {index} has an unexpected time")
        if min(row["rho"], row["pressure"], row["temperature"], row["rhoE"]) <= 0.0:
            raise AssertionError("selected reactive 3D state is nonphysical")

        maximum_relationship_error = max(
            maximum_relationship_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
        )
        mass_fraction_sum = 0.0
        species_density_sum = 0.0
        for name in species:
            mass_fraction = row[f"Y_{name}"]
            species_density = row[f"rhoY_{name}"]
            if mass_fraction < 0.0 or species_density < 0.0:
                raise AssertionError("selected reactive 3D species is negative")
            mass_fraction_sum += mass_fraction
            species_density_sum += species_density
            maximum_relationship_error = max(
                maximum_relationship_error,
                relative_error(species_density, row["rho"] * mass_fraction),
            )
        maximum_closure_error = max(
            maximum_closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )

    if maximum_relationship_error > 2.0e-11:
        raise AssertionError(
            f"selected reactive 3D relationship error {maximum_relationship_error}"
        )
    if maximum_closure_error > 2.0e-11:
        raise AssertionError(
            f"selected reactive 3D closure error {maximum_closure_error}"
        )

    state_columns = tuple(
        name for name in columns if name not in {"x", "y", "z"}
    )
    reference = rows[0]
    maximum_uniform_error = max(
        abs(row[name] - reference[name])
        for row in rows
        for name in state_columns
    )
    if args.require_uniform and maximum_uniform_error > 2.0e-13:
        raise AssertionError(
            f"selected reactive 3D uniformity error {maximum_uniform_error}"
        )

    maximum_requested_species_change = 0.0
    if args.activity_species is not None:
        if args.initial_mass_fraction is None:
            raise AssertionError("activity check requires an initial mass fraction")
        column = f"Y_{args.activity_species}"
        if column not in columns:
            raise AssertionError("activity species is not in the output schema")
        maximum_requested_species_change = max(
            abs(row[column] - args.initial_mass_fraction) for row in rows
        )
        if maximum_requested_species_change < args.minimum_change:
            raise AssertionError("selected reactive 3D chemistry was not active")

    print(f"minimum_density={min(row['rho'] for row in rows):.16e}")
    print(f"minimum_pressure={min(row['pressure'] for row in rows):.16e}")
    print(f"minimum_temperature={min(row['temperature'] for row in rows):.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    print(f"maximum_uniform_error={maximum_uniform_error:.16e}")
    print(
        "maximum_requested_species_change="
        f"{maximum_requested_species_change:.16e}"
    )
    print("selected reactive 3D regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
