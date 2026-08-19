#!/usr/bin/env python3
"""Check the periodic general-EOS H2/N2 composition-wave regression."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")
MOLECULAR_WEIGHTS = {
    "H2": 2.01588,
    "H": 1.00794,
    "O": 15.9994,
    "O2": 31.9988,
    "OH": 17.00734,
    "H2O": 18.01528,
    "N2": 28.0134,
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


def exact_state(x: float, time: float) -> tuple[float, dict[str, float]]:
    length = 0.012
    shifted = (x - 200.0 * time) % length
    phase = math.sin(2.0 * math.pi * shifted / length)
    mole = dict(BASE_MOLE_FRACTIONS)
    mole["H2"] += 0.04 * phase
    mole["N2"] -= 0.04 * phase
    denominator = sum(mole[name] * MOLECULAR_WEIGHTS[name] for name in SPECIES)
    mass = {
        name: mole[name] * MOLECULAR_WEIGHTS[name] / denominator
        for name in SPECIES
    }
    density = 101325.0 * denominator / (UNIVERSAL_GAS_CONSTANT * 1000.0)
    return density, mass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if len(rows) != 160:
        raise AssertionError("composition-wave output does not contain 160 cells")

    density_relative_l1 = 0.0
    h2_l1 = 0.0
    maximum_pressure_error = 0.0
    maximum_temperature_error = 0.0
    maximum_closure_error = 0.0
    for row in rows:
        exact_density, exact_mass = exact_state(row["x"], row["time"])
        density_relative_l1 += abs(row["rho"] - exact_density) / exact_density
        h2_l1 += abs(row["Y_H2"] - exact_mass["H2"])
        maximum_pressure_error = max(
            maximum_pressure_error, abs(row["pressure"] - 101325.0)
        )
        maximum_temperature_error = max(
            maximum_temperature_error, abs(row["temperature"] - 1000.0)
        )
        maximum_closure_error = max(
            maximum_closure_error,
            abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0),
        )

    density_relative_l1 /= len(rows)
    h2_l1 /= len(rows)
    print(f"density_relative_l1={density_relative_l1:.16e}")
    print(f"h2_mass_fraction_l1={h2_l1:.16e}")
    print(f"maximum_pressure_error={maximum_pressure_error:.16e}")
    print(f"maximum_temperature_error={maximum_temperature_error:.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")

    if density_relative_l1 > 2.0e-5:
        raise AssertionError("composition-wave density error is too large")
    if h2_l1 > 2.5e-6:
        raise AssertionError("composition-wave H2 error is too large")
    if maximum_pressure_error > 8.0:
        raise AssertionError("composition wave generated excessive pressure error")
    if maximum_temperature_error > 0.08:
        raise AssertionError("composition wave generated excessive temperature error")
    if maximum_closure_error > 5.0e-12:
        raise AssertionError("composition-wave mass fractions lost closure")
    print("Reactive composition-wave regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
