#!/usr/bin/env python3
"""Check the periodic molecular-transport pulse regression."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")
MOLECULAR_WEIGHTS = {
    "H2": 2.016,
    "H": 1.008,
    "O": 15.999,
    "O2": 31.998,
    "OH": 17.007,
    "H2O": 18.015,
    "N2": 28.014,
}
BASE_MOLE_FRACTIONS = {
    "H2": 0.29570,
    "H": 1.0e-5,
    "O": 1.0e-5,
    "O2": 0.14784,
    "OH": 1.0e-5,
    "H2O": 0.0,
    "N2": 0.55643,
}
UNIVERSAL_GAS_CONSTANT = 8314.46261815324


def initial_state(x: float, length: float) -> tuple[float, dict[str, float]]:
    phase = math.sin(2.0 * math.pi * x / length)
    mole = dict(BASE_MOLE_FRACTIONS)
    mole["H2"] += 0.02 * phase
    mole["N2"] -= 0.02 * phase
    molecular_weight = sum(mole[name] * MOLECULAR_WEIGHTS[name] for name in SPECIES)
    mass = {
        name: mole[name] * MOLECULAR_WEIGHTS[name] / molecular_weight
        for name in SPECIES
    }
    density = 101325.0 * molecular_weight / (
        UNIVERSAL_GAS_CONSTANT * 1000.0
    )
    return density, mass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--nx", type=int, default=96)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if len(rows) != args.nx:
        raise AssertionError("transport output has an unexpected number of cells")

    length = 0.012
    dx = length / args.nx
    initial_species_mass = {name: 0.0 for name in SPECIES}
    initial_h2 = []
    for index in range(args.nx):
        x = (index + 0.5) * dx
        density, mass = initial_state(x, length)
        initial_h2.append(mass["H2"])
        for name in SPECIES:
            initial_species_mass[name] += dx * density * mass[name]

    final_h2 = [row["Y_H2"] for row in rows]
    initial_h2_range = max(initial_h2) - min(initial_h2)
    final_h2_range = max(final_h2) - min(final_h2)
    maximum_closure_error = max(
        abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0) for row in rows
    )
    maximum_species_mass_error = 0.0
    for name in SPECIES:
        final_mass = dx * sum(row["rho"] * row[f"Y_{name}"] for row in rows)
        initial_mass = initial_species_mass[name]
        if initial_mass == 0.0:
            error = abs(final_mass)
        else:
            error = abs(final_mass - initial_mass) / abs(initial_mass)
        maximum_species_mass_error = max(maximum_species_mass_error, error)

    temperature_span = max(row["temperature"] for row in rows) - min(
        row["temperature"] for row in rows
    )
    pressure_span = max(row["pressure"] for row in rows) - min(
        row["pressure"] for row in rows
    )
    maximum_speed = max(abs(row["u"]) for row in rows)
    minimum_density = min(row["rho"] for row in rows)
    minimum_pressure = min(row["pressure"] for row in rows)
    minimum_temperature = min(row["temperature"] for row in rows)
    minimum_mass_fraction = min(
        row[f"Y_{name}"] for row in rows for name in SPECIES
    )

    print(f"initial_h2_range={initial_h2_range:.16e}")
    print(f"final_h2_range={final_h2_range:.16e}")
    print(f"maximum_species_mass_error={maximum_species_mass_error:.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    print(f"temperature_span={temperature_span:.16e}")
    print(f"pressure_span={pressure_span:.16e}")
    print(f"maximum_speed={maximum_speed:.16e}")

    if not (0.90 * initial_h2_range < final_h2_range < 0.99 * initial_h2_range):
        raise AssertionError("H2 wave did not undergo the qualified diffusion")
    if maximum_species_mass_error > 5.0e-12:
        raise AssertionError("periodic species masses are not conserved")
    if maximum_closure_error > 5.0e-12:
        raise AssertionError("mass fractions lost closure")
    if min(minimum_density, minimum_pressure, minimum_temperature) <= 0.0:
        raise AssertionError("transport pulse lost thermodynamic positivity")
    if minimum_mass_fraction < -1.0e-13:
        raise AssertionError("transport pulse produced a negative mass fraction")
    if not (0.01 < temperature_span < 0.20):
        raise AssertionError("transport temperature signature is outside bounds")
    if not (5.0 < pressure_span < 60.0):
        raise AssertionError("transport pressure signature is outside bounds")
    if not (1.0e-3 < maximum_speed < 0.10):
        raise AssertionError("transport-induced velocity signature is outside bounds")
    print("Reactive molecular-transport pulse regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
