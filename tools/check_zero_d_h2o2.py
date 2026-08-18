#!/usr/bin/env python3
"""Independent structural and conservation checks for the full H2/O2 reactor."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2")
MW = {
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
H_ATOMS = {
    "H2": 2.0, "H": 1.0, "O": 0.0, "O2": 0.0, "OH": 1.0,
    "H2O": 2.0, "HO2": 1.0, "H2O2": 2.0, "AR": 0.0, "N2": 0.0,
}
O_ATOMS = {
    "H2": 0.0, "H": 0.0, "O": 1.0, "O2": 2.0, "OH": 1.0,
    "H2O": 1.0, "HO2": 2.0, "H2O2": 2.0, "AR": 0.0, "N2": 0.0,
}


def atom_inventory(row: dict[str, float], atoms: dict[str, float]) -> float:
    return sum(atoms[name] * row[f"Y_{name}"] / MW[name] for name in SPECIES)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--final-time", type=float, default=2.0e-3)
    parser.add_argument("--energy-error-max", type=float, default=5.0e-9)
    parser.add_argument("--closure-error-max", type=float, default=5.0e-12)
    parser.add_argument("--atom-error-max", type=float, default=5.0e-10)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        raw_rows = list(csv.DictReader(handle))
    if len(raw_rows) < 3:
        raise AssertionError("H2/O2 output has too few rows")
    rows = [{key: float(value) for key, value in row.items()} for row in raw_rows]
    values = [value for row in rows for value in row.values()]
    if not all(math.isfinite(value) for value in values):
        raise AssertionError("H2/O2 output contains non-finite values")

    times = [row["time"] for row in rows]
    if abs(times[0]) > 1.0e-15 or abs(times[-1] - args.final_time) > 1.0e-12:
        raise AssertionError("H2/O2 time interval is incorrect")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise AssertionError("H2/O2 output times are not strictly increasing")

    densities = [row["density"] for row in rows]
    density_drift = max(abs(value - densities[0]) for value in densities)
    closure_error = max(
        abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0) for row in rows
    )
    reported_closure = max(row["closure_error"] for row in rows)
    energy_error = max(row["relative_energy_error"] for row in rows)
    minimum_mass_fraction = min(
        row[f"Y_{name}"] for row in rows for name in SPECIES
    )

    h0 = atom_inventory(rows[0], H_ATOMS)
    o0 = atom_inventory(rows[0], O_ATOMS)
    h_error = max(abs(atom_inventory(row, H_ATOMS) - h0) for row in rows)
    o_error = max(abs(atom_inventory(row, O_ATOMS) - o0) for row in rows)
    n2_drift = max(abs(row["Y_N2"] - rows[0]["Y_N2"]) for row in rows)
    ar_drift = max(abs(row["Y_AR"] - rows[0]["Y_AR"]) for row in rows)

    maximum_temperature_change = max(
        abs(row["temperature"] - rows[0]["temperature"]) for row in rows
    )
    maximum_species_change = max(
        abs(row[f"Y_{name}"] - rows[0][f"Y_{name}"])
        for row in rows
        for name in SPECIES
    )

    source_mass_error = 0.0
    source_h_error = 0.0
    source_o_error = 0.0
    maximum_rate = 0.0
    for row in rows:
        wdot = {name: row[f"wdot_{name}"] for name in SPECIES}
        maximum_rate = max(maximum_rate, *(abs(value) for value in wdot.values()))
        source_mass_error = max(
            source_mass_error,
            abs(sum(MW[name] * wdot[name] for name in SPECIES)),
        )
        source_h_error = max(
            source_h_error,
            abs(sum(H_ATOMS[name] * wdot[name] for name in SPECIES)),
        )
        source_o_error = max(
            source_o_error,
            abs(sum(O_ATOMS[name] * wdot[name] for name in SPECIES)),
        )

    metrics = {
        "rows": float(len(rows)),
        "density_drift": density_drift,
        "maximum_closure_error": max(closure_error, reported_closure),
        "maximum_energy_error": energy_error,
        "hydrogen_inventory_error": h_error,
        "oxygen_inventory_error": o_error,
        "n2_drift": n2_drift,
        "ar_drift": ar_drift,
        "minimum_mass_fraction": minimum_mass_fraction,
        "maximum_temperature_change": maximum_temperature_change,
        "maximum_species_change": maximum_species_change,
        "maximum_molar_rate": maximum_rate,
        "source_mass_error": source_mass_error,
        "source_h_error": source_h_error,
        "source_o_error": source_o_error,
        "final_temperature": rows[-1]["temperature"],
        "final_pressure": rows[-1]["pressure"],
    }
    for name, value in metrics.items():
        print(f"{name}={value:.16e}")

    if density_drift > 5.0e-14 * max(1.0, abs(densities[0])):
        raise AssertionError("constant-volume density drifted")
    if max(closure_error, reported_closure) > args.closure_error_max:
        raise AssertionError("mass-fraction closure exceeded tolerance")
    if energy_error > args.energy_error_max:
        raise AssertionError("internal-energy conservation exceeded tolerance")
    if max(h_error, o_error) > args.atom_error_max:
        raise AssertionError("element inventory exceeded tolerance")
    if max(n2_drift, ar_drift) > 5.0e-12:
        raise AssertionError("inert species changed")
    if minimum_mass_fraction < -1.0e-13:
        raise AssertionError("negative mass fraction generated")
    if maximum_temperature_change < 1.0 or maximum_species_change < 1.0e-7:
        raise AssertionError("full H2/O2 reactor did not evolve")
    source_scale = max(1.0, maximum_rate)
    if source_mass_error > 2.0e-12 * source_scale:
        raise AssertionError("instantaneous production rates violate mass conservation")
    if max(source_h_error, source_o_error) > 2.0e-12 * source_scale:
        raise AssertionError("instantaneous production rates violate atom conservation")

    print("Full H2/O2 structural regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
