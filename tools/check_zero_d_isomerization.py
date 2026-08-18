#!/usr/bin/env python3
"""Verify the synthetic adiabatic zero-dimensional isomerization case."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--initial-temperature", type=float, default=700.0)
    parser.add_argument("--final-temperature-reference", type=float, default=1840.0)
    parser.add_argument("--temperature-relative-tolerance", type=float, default=2.0e-10)
    parser.add_argument("--energy-error-max", type=float, default=2.0e-11)
    parser.add_argument("--closure-error-max", type=float, default=2.0e-12)
    parser.add_argument("--final-reactant-max", type=float, default=1.0e-8)
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) < 3:
        raise AssertionError("Expected at least three reactor output rows")

    required = {
        "time",
        "temperature",
        "pressure",
        "Y_A",
        "Y_B",
        "specific_internal_energy",
        "relative_energy_error",
        "rate_constant",
    }
    if set(rows[0]) != required:
        raise AssertionError(f"Unexpected CSV columns: {set(rows[0])}")

    values = [{key: float(value) for key, value in row.items()} for row in rows]
    if not all(math.isfinite(value) for row in values for value in row.values()):
        raise AssertionError("Reactor output contains a non-finite value")

    times = [row["time"] for row in values]
    temperatures = [row["temperature"] for row in values]
    reactant = [row["Y_A"] for row in values]
    product = [row["Y_B"] for row in values]
    pressures = [row["pressure"] for row in values]
    energy_errors = [row["relative_energy_error"] for row in values]

    if abs(times[0]) > 1.0e-15 or abs(times[-1] - args.final_time) > 1.0e-12:
        raise AssertionError("Unexpected reactor time interval")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise AssertionError("Reactor output times are not strictly increasing")
    if any(value <= 0.0 for value in temperatures + pressures):
        raise AssertionError("Reactor output contains non-positive T or p")
    if any(value < -1.0e-13 or value > 1.0 + 1.0e-13 for value in reactant + product):
        raise AssertionError("Reactor output contains an invalid mass fraction")

    closure_error = max(abs(a + b - 1.0) for a, b in zip(reactant, product))
    if closure_error > args.closure_error_max:
        raise AssertionError(f"Mass-fraction closure error {closure_error:.6e}")
    if any(right > left + 5.0e-14 for left, right in zip(reactant, reactant[1:])):
        raise AssertionError("Reactant mass fraction is not monotone decreasing")
    if any(right < left - 5.0e-11 for left, right in zip(temperatures, temperatures[1:])):
        raise AssertionError("Adiabatic temperature is not monotone increasing")

    final_temperature_error = abs(
        temperatures[-1] - args.final_temperature_reference
    ) / args.final_temperature_reference
    if final_temperature_error > args.temperature_relative_tolerance:
        raise AssertionError(
            f"Final-temperature relative error {final_temperature_error:.6e}"
        )
    if reactant[-1] > args.final_reactant_max:
        raise AssertionError(f"Final reactant fraction {reactant[-1]:.6e}")
    if max(energy_errors) > args.energy_error_max:
        raise AssertionError(f"Energy error {max(energy_errors):.6e}")
    if temperatures[-1] <= args.initial_temperature:
        raise AssertionError("Adiabatic reaction did not heat the mixture")

    metrics = {
        "rows": float(len(rows)),
        "final_temperature": temperatures[-1],
        "final_reactant_fraction": reactant[-1],
        "maximum_energy_error": max(energy_errors),
        "maximum_closure_error": closure_error,
        "final_pressure": pressures[-1],
    }
    for name, value in metrics.items():
        print(f"{name}={value:.16e}")
    print("Zero-dimensional isomerization regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
