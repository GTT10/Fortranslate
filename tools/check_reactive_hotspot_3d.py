#!/usr/bin/env python3
"""Validate the public elementary-chemistry 3D hotspot CSV."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")
MOLE_FRACTIONS = {
    "H2": 0.29570,
    "H": 1.0e-5,
    "O": 1.0e-5,
    "O2": 0.14784,
    "OH": 1.0e-5,
    "H2O": 0.0,
    "N2": 0.55643,
}
MOLECULAR_WEIGHTS = {
    "H2": 2.016,
    "H": 1.008,
    "O": 15.999,
    "O2": 31.998,
    "OH": 17.007,
    "H2O": 18.015,
    "N2": 28.014,
}
ATOMS = {
    "H2": (2, 0, 0),
    "H": (1, 0, 0),
    "O": (0, 1, 0),
    "O2": (0, 2, 0),
    "OH": (1, 1, 0),
    "H2O": (2, 1, 0),
    "N2": (0, 0, 2),
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
    return abs(actual - expected) / max(abs(expected), 1.0e-30)


def periodic_displacement(value: float, center: float, length: float) -> float:
    delta = value - center
    return delta - length * round(delta / length)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--time", type=float, required=True)
    args = parser.parse_args()

    columns = (
        *BASE_COLUMNS,
        *(f"Y_{name}" for name in SPECIES),
        *(f"rhoY_{name}" for name in SPECIES),
    )
    with args.input.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != columns:
            raise AssertionError(f"unexpected columns: {reader.fieldnames}")
        rows = [{name: float(row[name]) for name in columns} for row in reader]

    extents = (args.nx, args.ny, args.nz)
    if len(rows) != math.prod(extents):
        raise AssertionError("reactive hotspot row count does not match the grid")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("reactive hotspot CSV contains a non-finite value")

    length = 0.01
    for index, row in enumerate(rows):
        indices = (
            index % args.nx,
            (index // args.nx) % args.ny,
            index // (args.nx * args.ny),
        )
        expected_coordinates = tuple(
            (cell + 0.5) * length / extent
            for cell, extent in zip(indices, extents)
        )
        if max(
            abs(row[name] - expected)
            for name, expected in zip(("x", "y", "z"), expected_coordinates)
        ) > 5.0e-15:
            raise AssertionError(f"row {index} is not in x-fastest cell order")
        if abs(row["time"] - args.time) > 5.0e-18:
            raise AssertionError(f"row {index} has an unexpected output time")

    mean_molecular_weight = sum(
        MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name] for name in SPECIES
    )
    mass_fractions = {
        name: MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name]
        / mean_molecular_weight
        for name in SPECIES
    }
    gas_constant = 8.31446261815324e3 / mean_molecular_weight
    cell_volume = math.prod(length / extent for extent in extents)
    initial_mass = 0.0
    initial_elements = [0.0, 0.0, 0.0]
    for row in rows:
        radius_squared = sum(
            periodic_displacement(row[name], 0.005, length) ** 2
            for name in ("x", "y", "z")
        )
        initial_temperature = 1200.0 + 250.0 * math.exp(
            -0.5 * radius_squared / 0.0012**2
        )
        initial_density = 101325.0 / (gas_constant * initial_temperature)
        initial_mass += initial_density * cell_volume
        for species_name in SPECIES:
            species_amount = (
                initial_density
                * mass_fractions[species_name]
                / MOLECULAR_WEIGHTS[species_name]
                * cell_volume
            )
            for element in range(3):
                initial_elements[element] += (
                    ATOMS[species_name][element] * species_amount
                )

    relationship_error = 0.0
    closure_error = 0.0
    composition_change = 0.0
    final_elements = [0.0, 0.0, 0.0]
    for row in rows:
        if (
            row["rho"] <= 0.0
            or row["pressure"] <= 0.0
            or row["temperature"] <= 0.0
            or row["rhoE"] <= 0.0
        ):
            raise AssertionError("reactive hotspot contains a non-physical state")
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
        for species_name in SPECIES:
            mass_fraction = row[f"Y_{species_name}"]
            species_density = row[f"rhoY_{species_name}"]
            if mass_fraction < 0.0 or species_density < 0.0:
                raise AssertionError("reactive hotspot has a negative species")
            mass_fraction_sum += mass_fraction
            species_density_sum += species_density
            relationship_error = max(
                relationship_error,
                relative_error(species_density, row["rho"] * mass_fraction),
            )
            composition_change = max(
                composition_change,
                abs(mass_fraction - mass_fractions[species_name]),
            )
            species_amount = (
                species_density / MOLECULAR_WEIGHTS[species_name] * cell_volume
            )
            for element in range(3):
                final_elements[element] += ATOMS[species_name][element] * species_amount
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )

    if relationship_error > 2.0e-11 or closure_error > 2.0e-11:
        raise AssertionError(
            f"state contract error: relation={relationship_error}, closure={closure_error}"
        )
    if composition_change < 1.0e-6:
        raise AssertionError("reactive hotspot did not produce measurable chemistry")

    final_mass = sum(row["rho"] for row in rows) * cell_volume
    mass_error = relative_error(final_mass, initial_mass)
    element_error = max(
        relative_error(actual, expected)
        for actual, expected in zip(final_elements, initial_elements)
    )
    momentum_integrals = tuple(
        sum(row[name] for row in rows) * cell_volume
        for name in ("rhou", "rhov", "rhow")
    )
    if mass_error > 2.0e-11 or element_error > 2.0e-10:
        raise AssertionError(
            f"reactive hotspot conservation error: mass={mass_error}, elements={element_error}"
        )
    if max(abs(value) for value in momentum_integrals) > 2.0e-12:
        raise AssertionError(f"reactive hotspot net momentum {momentum_integrals}")
    temperature_range = max(row["temperature"] for row in rows) - min(
        row["temperature"] for row in rows
    )
    if temperature_range < 20.0:
        raise AssertionError("reactive hotspot lost its resolved temperature structure")

    print(
        "Reactive 3D hotspot: PASS "
        f"(rows={len(rows)}, chemistry={composition_change:.3e}, "
        f"mass={mass_error:.3e}, elements={element_error:.3e}, "
        f"closure={closure_error:.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
