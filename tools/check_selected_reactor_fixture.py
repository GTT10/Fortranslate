#!/usr/bin/env python3
"""Independent structure and conservation checks for the selected 0D fixture."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SPECIES = ("H2", "H")
MOLECULAR_WEIGHT = {"H2": 2.016, "H": 1.008}
H_ATOMS = {"H2": 2.0, "H": 1.0}
BASE_COLUMNS = (
    "time",
    "temperature",
    "pressure",
    "density",
    "specific_internal_energy",
    "relative_energy_error",
    "closure_error",
)
EXPECTED_COLUMNS = (
    *BASE_COLUMNS,
    *(f"Y_{name}" for name in SPECIES),
    *(f"wdot_{name}" for name in SPECIES),
)


def inventory(row: dict[str, float]) -> float:
    return sum(
        H_ATOMS[name] * row[f"Y_{name}"] / MOLECULAR_WEIGHT[name]
        for name in SPECIES
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--final-time", type=float, default=1.0e-7)
    parser.add_argument("--expected-rows", type=int, default=5)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != EXPECTED_COLUMNS:
            raise AssertionError("selected fixture CSV header is incorrect")
        raw_rows = list(reader)

    if len(raw_rows) != args.expected_rows:
        raise AssertionError(
            "selected fixture output row count is incorrect: "
            f"expected {args.expected_rows}, got {len(raw_rows)}"
        )
    if any(None in row for row in raw_rows):
        raise AssertionError("selected fixture CSV contains extra fields")
    rows = [{name: float(value) for name, value in row.items()} for row in raw_rows]
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("selected fixture output contains non-finite values")

    times = [row["time"] for row in rows]
    if abs(times[0]) > 1.0e-15 or abs(times[-1] - args.final_time) > 1.0e-14:
        raise AssertionError("selected fixture time interval is incorrect")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise AssertionError("selected fixture times are not strictly increasing")

    initial_expected = {"H2": 8.0 / 9.0, "H": 1.0 / 9.0}
    mapping_error = max(
        abs(rows[0][f"Y_{name}"] - initial_expected[name]) for name in SPECIES
    )
    densities = [row["density"] for row in rows]
    density_drift = max(abs(value - densities[0]) for value in densities)
    closure_error = max(
        max(abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0),
            row["closure_error"])
        for row in rows
    )
    energy_error = max(row["relative_energy_error"] for row in rows)
    minimum_mass_fraction = min(
        row[f"Y_{name}"] for row in rows for name in SPECIES
    )
    initial_inventory = inventory(rows[0])
    inventory_error = max(abs(inventory(row) - initial_inventory) for row in rows)
    maximum_temperature_change = max(
        abs(row["temperature"] - rows[0]["temperature"]) for row in rows
    )
    maximum_species_change = max(
        abs(row[f"Y_{name}"] - rows[0][f"Y_{name}"])
        for row in rows
        for name in SPECIES
    )

    maximum_rate = 0.0
    source_mass_error = 0.0
    source_atom_error = 0.0
    for row in rows:
        rates = {name: row[f"wdot_{name}"] for name in SPECIES}
        maximum_rate = max(maximum_rate, *(abs(value) for value in rates.values()))
        source_mass_error = max(
            source_mass_error,
            abs(sum(MOLECULAR_WEIGHT[name] * rates[name] for name in SPECIES)),
        )
        source_atom_error = max(
            source_atom_error,
            abs(sum(H_ATOMS[name] * rates[name] for name in SPECIES)),
        )

    metrics = {
        "rows": float(len(rows)),
        "mapping_error": mapping_error,
        "density_drift": density_drift,
        "closure_error": closure_error,
        "energy_error": energy_error,
        "hydrogen_inventory_error": inventory_error,
        "minimum_mass_fraction": minimum_mass_fraction,
        "maximum_temperature_change": maximum_temperature_change,
        "maximum_species_change": maximum_species_change,
        "maximum_molar_rate": maximum_rate,
        "source_mass_error": source_mass_error,
        "source_atom_error": source_atom_error,
    }
    for name, value in metrics.items():
        print(f"{name}={value:.16e}")

    if mapping_error > 5.0e-14:
        raise AssertionError("selected fixture composition was not mapped by name")
    if density_drift > 5.0e-14 * max(1.0, abs(densities[0])):
        raise AssertionError("selected fixture density drifted")
    if closure_error > 5.0e-12:
        raise AssertionError("selected fixture mass-fraction closure failed")
    if energy_error > 5.0e-9:
        raise AssertionError("selected fixture internal energy drifted")
    if inventory_error > 5.0e-10:
        raise AssertionError("selected fixture hydrogen inventory drifted")
    if minimum_mass_fraction < -1.0e-13:
        raise AssertionError("selected fixture generated a negative mass fraction")
    if maximum_temperature_change < 1.0e-8 or maximum_species_change < 1.0e-12:
        raise AssertionError("selected fixture reactor did not evolve")
    if maximum_rate <= 1.0e-12:
        raise AssertionError("selected fixture production rates are inactive")
    source_scale = max(1.0, maximum_rate)
    if source_mass_error > 2.0e-12 * source_scale:
        raise AssertionError("selected fixture rates violate mass conservation")
    if source_atom_error > 2.0e-12 * source_scale:
        raise AssertionError("selected fixture rates violate atom conservation")

    print("Selected reactor fixture regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
