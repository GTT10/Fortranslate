#!/usr/bin/env python3
"""Validate the public 3D entropy-wave CSV and its conserved integrals."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


COLUMNS = (
    "x",
    "y",
    "z",
    "rho",
    "u",
    "v",
    "w",
    "p",
    "total_energy_density",
    "momentum_x_density",
    "momentum_y_density",
    "momentum_z_density",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--maximum-l1", type=float, required=True)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise AssertionError(f"unexpected columns: {reader.fieldnames}")
        rows = [{name: float(row[name]) for name in COLUMNS} for row in reader]

    expected_rows = args.nx * args.ny * args.nz
    if len(rows) != expected_rows:
        raise AssertionError(f"expected {expected_rows} rows, found {len(rows)}")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("CSV contains a non-finite value")

    unique_x = sorted({row["x"] for row in rows})
    unique_y = sorted({row["y"] for row in rows})
    unique_z = sorted({row["z"] for row in rows})
    if (len(unique_x), len(unique_y), len(unique_z)) != (
        args.nx,
        args.ny,
        args.nz,
    ):
        raise AssertionError("CSV coordinate cardinality does not match the grid")
    for index, row in enumerate(rows):
        i = index % args.nx
        j = (index // args.nx) % args.ny
        k = index // (args.nx * args.ny)
        expected_coordinates = (
            (i + 0.5) / args.nx,
            (j + 0.5) / args.ny,
            (k + 0.5) / args.nz,
        )
        actual_coordinates = (row["x"], row["y"], row["z"])
        if max(
            abs(actual - expected)
            for actual, expected in zip(actual_coordinates, expected_coordinates)
        ) > 5.0e-15:
            raise AssertionError(f"row {index} is not in x-fastest cell order")

    gamma = 1.4
    rho0 = 1.0
    pressure = 1.0
    velocity = (0.7, 0.2, -0.1)
    kinetic_speed_squared = sum(component * component for component in velocity)
    phase_speed = sum(velocity)
    density_error = 0.0
    relation_error = 0.0
    for row in rows:
        exact_density = rho0 + 0.1 * math.sin(
            2.0
            * math.pi
            * (row["x"] + row["y"] + row["z"] - args.time * phase_speed)
        )
        density_error += abs(row["rho"] - exact_density)
        if row["rho"] <= 0.0 or row["p"] <= 0.0:
            raise AssertionError("CSV contains a non-physical state")
        relation_error = max(
            relation_error,
            abs(row["u"] - velocity[0]),
            abs(row["v"] - velocity[1]),
            abs(row["w"] - velocity[2]),
            abs(row["p"] - pressure),
            abs(row["momentum_x_density"] - row["rho"] * velocity[0]),
            abs(row["momentum_y_density"] - row["rho"] * velocity[1]),
            abs(row["momentum_z_density"] - row["rho"] * velocity[2]),
            abs(
                row["total_energy_density"]
                - pressure / (gamma - 1.0)
                - 0.5 * row["rho"] * kinetic_speed_squared
            ),
        )
    density_error /= len(rows)
    if density_error > args.maximum_l1:
        raise AssertionError(
            f"density L1 {density_error:.16e} exceeds {args.maximum_l1:.16e}"
        )
    if relation_error > 5.0e-12:
        raise AssertionError(f"primitive/conserved relation error {relation_error}")

    dx = 1.0 / args.nx
    dy = 1.0 / args.ny
    dz = 1.0 / args.nz
    cell_volume = dx * dy * dz
    actual_totals = (
        sum(row["rho"] for row in rows) * cell_volume,
        sum(row["momentum_x_density"] for row in rows) * cell_volume,
        sum(row["momentum_y_density"] for row in rows) * cell_volume,
        sum(row["momentum_z_density"] for row in rows) * cell_volume,
        sum(row["total_energy_density"] for row in rows) * cell_volume,
    )
    expected_totals = (
        rho0,
        rho0 * velocity[0],
        rho0 * velocity[1],
        rho0 * velocity[2],
        pressure / (gamma - 1.0) + 0.5 * rho0 * kinetic_speed_squared,
    )
    conservation_error = max(
        abs(actual - expected) / max(1.0, abs(expected))
        for actual, expected in zip(actual_totals, expected_totals)
    )
    if conservation_error > 2.0e-11:
        raise AssertionError(f"conservation error {conservation_error}")

    print(
        "3D entropy wave: PASS "
        f"(rows={len(rows)}, L1={density_error:.8e}, "
        f"conservation={conservation_error:.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
