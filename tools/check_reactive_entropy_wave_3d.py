#!/usr/bin/env python3
"""Validate the public general-EOS multispecies 3D entropy-wave CSV."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SPECIES_BY_MODEL = {
    "elementary": ("H2", "H", "O", "O2", "OH", "H2O", "N2"),
    "full_h2o2": (
        "H2",
        "H",
        "O",
        "O2",
        "OH",
        "H2O",
        "HO2",
        "H2O2",
        "AR",
        "N2",
    ),
}

MOLE_FRACTIONS = {
    "H2": 0.29570,
    "H": 1.0e-5,
    "O": 1.0e-5,
    "O2": 0.14784,
    "OH": 1.0e-5,
    "H2O": 0.0,
    "HO2": 0.0,
    "H2O2": 0.0,
    "AR": 0.0,
    "N2": 0.55643,
}

MOLECULAR_WEIGHTS = {
    "H2": 2.016,
    "H": 1.008,
    "O": 15.999,
    "O2": 31.998,
    "OH": 17.007,
    "H2O": 18.015,
    "HO2": 33.006,
    "H2O2": 34.014,
    "AR": 39.950,
    "N2": 28.014,
}

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
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--thermo-model", choices=SPECIES_BY_MODEL, required=True)
    parser.add_argument("--maximum-l1", type=float, required=True)
    args = parser.parse_args()

    species = SPECIES_BY_MODEL[args.thermo_model]
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

    expected_rows = args.nx * args.ny * args.nz
    if len(rows) != expected_rows:
        raise AssertionError(f"expected {expected_rows} rows, found {len(rows)}")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("CSV contains a non-finite value")

    lower = (0.0, 0.0, 0.0)
    lengths = (0.01, 0.01, 0.01)
    extents = (args.nx, args.ny, args.nz)
    unique_coordinates = tuple(
        sorted({row[name] for row in rows}) for name in ("x", "y", "z")
    )
    if tuple(len(values) for values in unique_coordinates) != extents:
        raise AssertionError("CSV coordinate cardinality does not match the grid")
    for index, row in enumerate(rows):
        indices = (
            index % args.nx,
            (index // args.nx) % args.ny,
            index // (args.nx * args.ny),
        )
        expected_coordinates = tuple(
            lo + (cell + 0.5) * length / extent
            for lo, cell, length, extent in zip(lower, indices, lengths, extents)
        )
        actual_coordinates = (row["x"], row["y"], row["z"])
        if max(
            abs(actual - expected)
            for actual, expected in zip(actual_coordinates, expected_coordinates)
        ) > 5.0e-15:
            raise AssertionError(f"row {index} is not in x-fastest cell order")
        if abs(row["time"] - args.time) > 5.0e-18:
            raise AssertionError(f"row {index} has an unexpected output time")

    mean_molecular_weight = sum(
        MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name] for name in species
    )
    denominator = mean_molecular_weight
    expected_y = {
        name: MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name] / denominator
        for name in species
    }
    gas_constant = 8.31446261815324e3 / mean_molecular_weight
    pressure0 = 101325.0
    temperature0 = 1000.0
    base_density = pressure0 / (gas_constant * temperature0)
    velocity = (300.0, 200.0, -100.0)
    phase_speed = sum(component / length for component, length in zip(velocity, lengths))

    density_error = 0.0
    relationship_error = 0.0
    composition_error = 0.0
    closure_error = 0.0
    pressure_error = 0.0
    velocity_error = 0.0
    for row in rows:
        phase = 2.0 * math.pi * (
            row["x"] / lengths[0]
            + row["y"] / lengths[1]
            + row["z"] / lengths[2]
            - args.time * phase_speed
        )
        exact_density = base_density * (1.0 + 0.08 * math.sin(phase))
        density_error += abs(row["rho"] - exact_density)
        if row["rho"] <= 0.0 or row["pressure"] <= 0.0 or row["temperature"] <= 0.0:
            raise AssertionError("CSV contains a non-physical state")
        relationship_error = max(
            relationship_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
            relative_error(
                row["pressure"], row["rho"] * gas_constant * row["temperature"]
            ),
        )
        mass_fraction_sum = 0.0
        species_density_sum = 0.0
        for name in species:
            mass_fraction = row[f"Y_{name}"]
            species_density = row[f"rhoY_{name}"]
            mass_fraction_sum += mass_fraction
            species_density_sum += species_density
            composition_error = max(
                composition_error,
                abs(mass_fraction - expected_y[name]),
                relative_error(species_density, row["rho"] * mass_fraction),
            )
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )
        pressure_error = max(pressure_error, relative_error(row["pressure"], pressure0))
        velocity_error = max(
            velocity_error,
            *(abs(row[name] - expected) for name, expected in zip(("u", "v", "w"), velocity)),
        )

    density_error /= len(rows)
    if density_error > args.maximum_l1:
        raise AssertionError(
            f"density L1 {density_error:.16e} exceeds {args.maximum_l1:.16e}"
        )
    if relationship_error > 2.0e-12:
        raise AssertionError(f"primitive/conserved relationship error {relationship_error}")
    if composition_error > 2.0e-12 or closure_error > 2.0e-12:
        raise AssertionError(
            f"species contract error: composition={composition_error}, closure={closure_error}"
        )
    if pressure_error > 5.0e-4 or velocity_error > 5.0e-2:
        raise AssertionError(
            f"manufactured-state drift: pressure={pressure_error}, velocity={velocity_error}"
        )

    cell_volume = math.prod(length / extent for length, extent in zip(lengths, extents))
    expected_integrals = {
        "rho": base_density,
        "rhou": base_density * velocity[0],
        "rhov": base_density * velocity[1],
        "rhow": base_density * velocity[2],
    }
    integral_error = 0.0
    volume = math.prod(lengths)
    for name, mean_value in expected_integrals.items():
        actual = sum(row[name] for row in rows) * cell_volume
        expected = mean_value * volume
        integral_error = max(integral_error, relative_error(actual, expected))
    for name in species:
        actual = sum(row[f"rhoY_{name}"] for row in rows) * cell_volume
        expected = base_density * expected_y[name] * volume
        integral_error = max(integral_error, relative_error(actual, expected))
    if integral_error > 2.0e-11:
        raise AssertionError(f"periodic conserved-integral error {integral_error}")

    print(
        "Reactive 3D entropy wave: PASS "
        f"(rows={len(rows)}, L1={density_error:.8e}, "
        f"integrals={integral_error:.3e}, closure={closure_error:.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
