#!/usr/bin/env python3
"""Validate dynamic-schema output from a selected planar reactive EB 3D run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


BASE_COLUMNS = (
    "time",
    "i",
    "j",
    "k",
    "x",
    "y",
    "z",
    "cell_type",
    "volume_fraction",
    "fluid_centroid_x",
    "fluid_centroid_y",
    "fluid_centroid_z",
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


def expected_geometry(
    indices: tuple[int, int, int],
    centers: tuple[float, float, float],
    spacings: tuple[float, float, float],
    plane_axis: int,
    plane_position: float,
) -> tuple[int, float, tuple[float, float, float]]:
    axis_center = centers[plane_axis]
    half_spacing = 0.5 * spacings[plane_axis]
    lower_face = axis_center - half_spacing
    upper_face = axis_center + half_spacing
    centroid = list(centers)
    tolerance = 128.0 * math.ulp(max(1.0, abs(plane_position)))
    if upper_face <= plane_position + tolerance:
        return 0, 0.0, tuple(centroid)
    if lower_face >= plane_position - tolerance:
        return 2, 1.0, tuple(centroid)
    fluid_length = upper_face - plane_position
    volume_fraction = fluid_length / spacings[plane_axis]
    centroid[plane_axis] = plane_position + 0.5 * fluid_length
    return 1, volume_fraction, tuple(centroid)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--x-lower", type=float, default=0.0)
    parser.add_argument("--x-upper", type=float, required=True)
    parser.add_argument("--y-lower", type=float, default=0.0)
    parser.add_argument("--y-upper", type=float, required=True)
    parser.add_argument("--z-lower", type=float, default=0.0)
    parser.add_argument("--z-upper", type=float, required=True)
    parser.add_argument("--plane-axis", choices=("x", "y", "z"), required=True)
    parser.add_argument("--plane-position", type=float, required=True)
    parser.add_argument("--expected-regular", type=int, required=True)
    parser.add_argument("--expected-cut", type=int, required=True)
    parser.add_argument("--expected-covered", type=int, required=True)
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    parser.add_argument("--minimum-active-density-span", type=float, default=0.0)
    parser.add_argument("--expected-sha256")
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    species = tuple(args.species)
    columns = (
        *BASE_COLUMNS,
        *(f"Y_{name}" for name in species),
        *(f"rhoY_{name}" for name in species),
    )
    if args.expected_sha256 is not None:
        actual_hash = hashlib.sha256(args.input.read_bytes()).hexdigest()
        if actual_hash != args.expected_sha256.lower():
            raise AssertionError(
                f"selected reactive EB 3D hash mismatch: {actual_hash}"
            )
    with args.input.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != columns:
            raise AssertionError(
                f"selected reactive EB 3D schema mismatch: {reader.fieldnames}"
            )
        rows = [{name: float(row[name]) for name in columns} for row in reader]

    extents = (args.nx, args.ny, args.nz)
    bounds = (
        (args.x_lower, args.x_upper),
        (args.y_lower, args.y_upper),
        (args.z_lower, args.z_upper),
    )
    spacings = tuple(
        (upper - lower) / extent
        for extent, (lower, upper) in zip(extents, bounds)
    )
    axis = {"x": 0, "y": 1, "z": 2}[args.plane_axis]
    if len(rows) != math.prod(extents):
        raise AssertionError("selected reactive EB 3D row count mismatch")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("selected reactive EB 3D output contains nonfinite data")

    type_counts = {0: 0, 1: 0, 2: 0}
    maximum_relationship_error = 0.0
    maximum_closure_error = 0.0
    maximum_activity_change = 0.0
    active_densities: list[float] = []
    for row_index, row in enumerate(rows):
        zero_based = (
            row_index % args.nx,
            (row_index // args.nx) % args.ny,
            row_index // (args.nx * args.ny),
        )
        one_based = tuple(index + 1 for index in zero_based)
        if tuple(round(row[name]) for name in ("i", "j", "k")) != one_based:
            raise AssertionError(f"row {row_index} is not in x-fastest index order")
        centers = tuple(
            lower + (index + 0.5) * spacing
            for index, spacing, (lower, _) in zip(zero_based, spacings, bounds)
        )
        if max(
            abs(row[name] - expected)
            for name, expected in zip(("x", "y", "z"), centers)
        ) > 5.0e-14:
            raise AssertionError(f"row {row_index} has a bad Cartesian center")
        cell_type, volume_fraction, centroid = expected_geometry(
            zero_based, centers, spacings, axis, args.plane_position
        )
        observed_type = round(row["cell_type"])
        if row["cell_type"] != observed_type or observed_type != cell_type:
            raise AssertionError(f"row {row_index} has a bad cell type")
        if abs(row["volume_fraction"] - volume_fraction) > 5.0e-13:
            raise AssertionError(f"row {row_index} has a bad volume fraction")
        if max(
            abs(row[name] - expected)
            for name, expected in zip(
                ("fluid_centroid_x", "fluid_centroid_y", "fluid_centroid_z"),
                centroid,
            )
        ) > 5.0e-13:
            raise AssertionError(f"row {row_index} has a bad fluid centroid")
        type_counts[cell_type] += 1
        if abs(row["time"] - args.final_time) > 5.0e-18:
            raise AssertionError(f"row {row_index} has an unexpected time")
        if min(row["rho"], row["pressure"], row["temperature"], row["rhoE"]) <= 0.0:
            raise AssertionError("selected reactive EB 3D state is nonphysical")
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
            if mass_fraction < -1.0e-13 or species_density < -1.0e-13:
                raise AssertionError("selected reactive EB 3D species is negative")
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
        if cell_type != 0:
            active_densities.append(row["rho"])
            if args.activity_species is not None:
                if args.initial_mass_fraction is None:
                    raise AssertionError(
                        "activity check requires an initial mass fraction"
                    )
                maximum_activity_change = max(
                    maximum_activity_change,
                    abs(
                        row[f"Y_{args.activity_species}"]
                        - args.initial_mass_fraction
                    ),
                )

    expected_counts = {
        2: args.expected_regular,
        1: args.expected_cut,
        0: args.expected_covered,
    }
    if type_counts != expected_counts:
        raise AssertionError(
            "selected reactive EB 3D cell counts mismatch: "
            f"observed={type_counts}, expected={expected_counts}"
        )
    if maximum_relationship_error > 2.0e-11:
        raise AssertionError(
            f"selected reactive EB 3D relationship error {maximum_relationship_error}"
        )
    if maximum_closure_error > 2.0e-11:
        raise AssertionError(
            f"selected reactive EB 3D closure error {maximum_closure_error}"
        )
    if maximum_activity_change < args.minimum_change:
        raise AssertionError("selected reactive EB 3D activity is below the gate")
    density_span = max(active_densities) - min(active_densities)
    if density_span < args.minimum_active_density_span:
        raise AssertionError("selected reactive EB 3D density response is unresolved")

    print(
        "cells="
        f"regular:{type_counts[2]},cut:{type_counts[1]},covered:{type_counts[0]}"
    )
    print(f"minimum_density={min(active_densities):.16e}")
    print(f"minimum_pressure={min(row['pressure'] for row in rows):.16e}")
    print(f"minimum_temperature={min(row['temperature'] for row in rows):.16e}")
    print(f"active_density_span={density_span:.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    print(f"maximum_requested_species_change={maximum_activity_change:.16e}")
    print("selected reactive EB 3D regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
