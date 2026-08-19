#!/usr/bin/env python3
"""Compare PeleF full H2/O2 rates and reactor history with Cantera."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import cantera as ct

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2")


@dataclass
class WorstCase:
    scaled: float = 0.0
    absolute: float = 0.0
    time: float = 0.0
    name: str = ""
    value: float = 0.0
    reference: float = 0.0


def scaled_error(value: float, reference: float, absolute: float, relative: float) -> float:
    return abs(value - reference) / (absolute + relative * max(abs(value), abs(reference)))


def update_worst(
    worst: WorstCase,
    *,
    value: float,
    reference: float,
    time: float,
    name: str,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> None:
    absolute_error = abs(value - reference)
    scaled = scaled_error(value, reference, absolute_tolerance, relative_tolerance)
    if scaled > worst.scaled:
        worst.scaled = scaled
        worst.absolute = absolute_error
        worst.time = time
        worst.name = name
        worst.value = value
        worst.reference = reference


def print_worst(prefix: str, worst: WorstCase) -> None:
    print(f"{prefix}_name={worst.name}")
    print(f"{prefix}_time={worst.time:.16e}")
    print(f"{prefix}_pelef={worst.value:.16e}")
    print(f"{prefix}_cantera={worst.reference:.16e}")
    print(f"{prefix}_absolute_error={worst.absolute:.16e}")
    print(f"{prefix}_scaled_error={worst.scaled:.16e}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--mechanism", default="h2o2.yaml")
    parser.add_argument("--phase", default="ohmech")
    parser.add_argument("--thermo-relative-tolerance", type=float, default=5.0e-5)
    parser.add_argument("--species-relative-tolerance", type=float, default=2.0e-4)
    parser.add_argument("--species-absolute-tolerance", type=float, default=2.0e-10)
    parser.add_argument("--rate-relative-tolerance", type=float, default=2.0e-8)
    parser.add_argument("--rate-absolute-tolerance", type=float, default=1.0e-10)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise AssertionError("PeleF full H2/O2 CSV is empty")

    gas = ct.Solution(args.mechanism, args.phase)
    missing = [name for name in SPECIES if name not in gas.species_names]
    if missing:
        raise AssertionError(f"Cantera mechanism is missing species: {missing}")
    initial_y = {name: rows[0][f"Y_{name}"] for name in SPECIES}
    gas.TPY = rows[0]["temperature"], rows[0]["pressure"], initial_y
    reactor = ct.IdealGasReactor(gas, energy="on", clone=True)
    reactor.volume = 1.0
    network = ct.ReactorNet([reactor])
    network.rtol = 1.0e-11
    network.atol = 1.0e-18

    rate_gas = ct.Solution(args.mechanism, args.phase)
    species_indices = {name: rate_gas.species_index(name) for name in SPECIES}
    worst_temperature = WorstCase(name="temperature")
    worst_pressure = WorstCase(name="pressure")
    worst_species = WorstCase()
    worst_rate = WorstCase()
    maximum_reference_rate = 0.0
    maximum_absolute_rate_error = 0.0

    for row in rows:
        if row["time"] > 0.0:
            network.advance(row["time"])
        phase = reactor.phase
        update_worst(
            worst_temperature,
            value=row["temperature"],
            reference=phase.T,
            time=row["time"],
            name="temperature",
            absolute_tolerance=1.0e-6,
            relative_tolerance=args.thermo_relative_tolerance,
        )
        update_worst(
            worst_pressure,
            value=row["pressure"],
            reference=phase.P,
            time=row["time"],
            name="pressure",
            absolute_tolerance=1.0e-3,
            relative_tolerance=args.thermo_relative_tolerance,
        )
        for name in SPECIES:
            update_worst(
                worst_species,
                value=row[f"Y_{name}"],
                reference=phase[name].Y[0],
                time=row["time"],
                name=name,
                absolute_tolerance=args.species_absolute_tolerance,
                relative_tolerance=args.species_relative_tolerance,
            )

        # Evaluate rates at the exact PeleF state. This isolates the generated
        # production-rate kernel from trajectory divergence between integrators.
        rate_gas.TDY = (
            row["temperature"],
            row["density"],
            {name: row[f"Y_{name}"] for name in SPECIES},
        )
        reference_rates = rate_gas.net_production_rates
        for name in SPECIES:
            value = row[f"wdot_{name}"]
            reference = reference_rates[species_indices[name]]
            maximum_reference_rate = max(maximum_reference_rate, abs(reference))
            maximum_absolute_rate_error = max(
                maximum_absolute_rate_error, abs(value - reference)
            )
            update_worst(
                worst_rate,
                value=value,
                reference=reference,
                time=row["time"],
                name=name,
                absolute_tolerance=args.rate_absolute_tolerance,
                relative_tolerance=args.rate_relative_tolerance,
            )

    print(f"cantera_version={ct.__version__}")
    print_worst("worst_temperature", worst_temperature)
    print_worst("worst_pressure", worst_pressure)
    print_worst("worst_species", worst_species)
    print_worst("worst_production_rate", worst_rate)
    print(f"maximum_reference_rate={maximum_reference_rate:.16e}")
    print(f"maximum_absolute_rate_error={maximum_absolute_rate_error:.16e}")
    print(f"cantera_final_temperature={reactor.phase.T:.16e}")
    print(f"pelef_final_temperature={rows[-1]['temperature']:.16e}")

    if worst_temperature.scaled > 1.0:
        raise AssertionError("temperature history disagrees with Cantera")
    if worst_pressure.scaled > 1.0:
        raise AssertionError("pressure history disagrees with Cantera")
    if worst_species.scaled > 1.0:
        raise AssertionError("species history disagrees with Cantera")
    if worst_rate.scaled > 1.0:
        raise AssertionError("production rates disagree with Cantera")

    print("Cantera full H2/O2 parity: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
