#!/usr/bin/env python3
"""Verify passive-species transport in the constant-gamma MultiSpecSod case."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from compare_sod import (
    PrimitiveState,
    exact_state_at,
    l1_error,
    solve_star_region,
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--nx", type=int, default=400)
    parser.add_argument("--time", type=float, default=0.2)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=1.0)
    parser.add_argument("--discontinuity", type=float, default=0.5)
    parser.add_argument("--density-l1-max", type=float, default=1.7e-3)
    parser.add_argument("--pressure-l1-max", type=float, default=1.0e-3)
    parser.add_argument("--mass-fraction-l1-max", type=float, default=7.0e-3)
    parser.add_argument("--species1-density-l1-max", type=float, default=3.0e-3)
    parser.add_argument("--species2-density-l1-max", type=float, default=2.5e-3)
    parser.add_argument("--conservation-error-max", type=float, default=2.0e-12)
    parser.add_argument("--closure-error-max", type=float, default=2.0e-12)
    parser.add_argument("--transition-cells-max", type=int, default=20)
    parser.add_argument("--far-leakage-max", type=float, default=1.0e-8)
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != args.nx:
        raise AssertionError(f"Expected {args.nx} rows, found {len(rows)}")

    required = {
        "x",
        "rho",
        "u",
        "p",
        "internal_energy_density",
        "temperature",
        "total_energy_density",
        "rhoY1",
        "Y1",
        "rhoY2",
        "Y2",
    }
    if not required.issubset(rows[0]):
        missing = sorted(required.difference(rows[0]))
        raise AssertionError(f"Missing output columns: {missing}")

    x = [float(row["x"]) for row in rows]
    rho = [float(row["rho"]) for row in rows]
    velocity = [float(row["u"]) for row in rows]
    pressure = [float(row["p"]) for row in rows]
    internal_energy_density = [
        float(row["internal_energy_density"]) for row in rows
    ]
    temperature = [float(row["temperature"]) for row in rows]
    total_energy_density = [
        float(row["total_energy_density"]) for row in rows
    ]
    rho_y1 = [float(row["rhoY1"]) for row in rows]
    rho_y2 = [float(row["rhoY2"]) for row in rows]
    y1 = [float(row["Y1"]) for row in rows]
    y2 = [float(row["Y2"]) for row in rows]

    if not all(
        math.isfinite(value)
        for row in rows
        for value in map(float, row.values())
    ):
        raise AssertionError("Output contains a non-finite value")
    if min(rho) <= 0.0 or min(pressure) <= 0.0:
        raise AssertionError("Output contains non-positive density or pressure")

    left = PrimitiveState(1.0, 0.0, 1.0)
    right = PrimitiveState(0.125, 0.0, 0.1)
    star_pressure, star_velocity = solve_star_region(left, right, args.gamma)
    contact_position = args.discontinuity + star_velocity * args.time
    exact = [
        exact_state_at(
            (position - args.discontinuity) / args.time,
            left,
            right,
            args.gamma,
            star_pressure,
            star_velocity,
        )
        for position in x
    ]
    exact_y1 = [1.0 if position <= contact_position else 0.0 for position in x]
    exact_rho_y1 = [state.rho * value for state, value in zip(exact, exact_y1)]
    exact_rho_y2 = [
        state.rho * (1.0 - value) for state, value in zip(exact, exact_y1)
    ]

    density_l1 = l1_error(rho, (state.rho for state in exact))
    pressure_l1 = l1_error(pressure, (state.pressure for state in exact))
    y1_l1 = l1_error(y1, exact_y1)
    rho_y1_l1 = l1_error(rho_y1, exact_rho_y1)
    rho_y2_l1 = l1_error(rho_y2, exact_rho_y2)

    closure_error = max(
        abs(first + second - total)
        for first, second, total in zip(rho_y1, rho_y2, rho)
    )
    mass_fraction_closure_error = max(
        abs(first + second - 1.0) for first, second in zip(y1, y2)
    )
    internal_energy_layout_error = max(
        abs(value - local_pressure / (args.gamma - 1.0))
        for value, local_pressure in zip(internal_energy_density, pressure)
    )
    temperature_layout_error = max(
        abs(value - local_pressure / local_density)
        for value, local_pressure, local_density in zip(
            temperature, pressure, rho
        )
    )
    total_energy_layout_error = max(
        abs(
            value
            - local_pressure / (args.gamma - 1.0)
            - 0.5 * local_density * local_velocity * local_velocity
        )
        for value, local_pressure, local_density, local_velocity in zip(
            total_energy_density, pressure, rho, velocity
        )
    )

    dx = (args.x_max - args.x_min) / args.nx
    species1_mass_error = abs(sum(rho_y1) * dx - 0.5)
    species2_mass_error = abs(sum(rho_y2) * dx - 0.0625)
    transition_cells = sum(1 for value in y1 if 1.0e-3 < value < 1.0 - 1.0e-3)
    left_leakage = max(
        y2[index]
        for index, position in enumerate(x)
        if position < contact_position - 0.05
    )
    right_leakage = max(
        y1[index]
        for index, position in enumerate(x)
        if position > contact_position + 0.05
    )

    metrics = {
        "density_l1": density_l1,
        "pressure_l1": pressure_l1,
        "mass_fraction_l1": y1_l1,
        "species1_density_l1": rho_y1_l1,
        "species2_density_l1": rho_y2_l1,
        "species_density_closure_error": closure_error,
        "mass_fraction_closure_error": mass_fraction_closure_error,
        "internal_energy_layout_error": internal_energy_layout_error,
        "temperature_layout_error": temperature_layout_error,
        "total_energy_layout_error": total_energy_layout_error,
        "species1_mass_error": species1_mass_error,
        "species2_mass_error": species2_mass_error,
        "minimum_mass_fraction": min(min(y1), min(y2)),
        "maximum_mass_fraction": max(max(y1), max(y2)),
        "transition_cells": float(transition_cells),
        "left_far_leakage": left_leakage,
        "right_far_leakage": right_leakage,
        "contact_position": contact_position,
    }
    for name, value in metrics.items():
        print(f"{name}={value:.16e}")

    failures: list[str] = []
    limits = {
        "density_l1": args.density_l1_max,
        "pressure_l1": args.pressure_l1_max,
        "mass_fraction_l1": args.mass_fraction_l1_max,
        "species1_density_l1": args.species1_density_l1_max,
        "species2_density_l1": args.species2_density_l1_max,
        "species_density_closure_error": args.closure_error_max,
        "mass_fraction_closure_error": args.closure_error_max,
        "internal_energy_layout_error": args.closure_error_max,
        "temperature_layout_error": args.closure_error_max,
        "total_energy_layout_error": args.closure_error_max,
        "species1_mass_error": args.conservation_error_max,
        "species2_mass_error": args.conservation_error_max,
        "left_far_leakage": args.far_leakage_max,
        "right_far_leakage": args.far_leakage_max,
    }
    for name, limit in limits.items():
        if metrics[name] > limit:
            failures.append(f"{name}={metrics[name]:.6e} exceeds {limit:.6e}")
    if metrics["minimum_mass_fraction"] < -1.0e-12:
        failures.append("negative mass fraction")
    if metrics["maximum_mass_fraction"] > 1.0 + 1.0e-12:
        failures.append("mass fraction exceeds one")
    if transition_cells > args.transition_cells_max:
        failures.append(
            f"transition_cells={transition_cells} exceeds {args.transition_cells_max}"
        )
    if failures:
        raise AssertionError("; ".join(failures))

    print("MultiSpecSod regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
