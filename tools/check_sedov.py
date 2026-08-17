#!/usr/bin/env python3
"""Verify the symmetric one-dimensional Sedov-type blast regression."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--nx", type=int, default=800)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=1.0)
    parser.add_argument("--blast-center", type=float, default=0.5)
    parser.add_argument("--blast-radius", type=float, default=0.0125)
    parser.add_argument("--ambient-density", type=float, default=1.0)
    parser.add_argument("--ambient-pressure", type=float, default=1.0e-5)
    parser.add_argument("--blast-pressure", type=float, default=100.0)
    parser.add_argument("--initial-velocity", type=float, default=0.0)
    parser.add_argument("--balance-error-max", type=float, default=2.0e-9)
    parser.add_argument("--symmetry-error-max", type=float, default=2.0e-10)
    parser.add_argument("--maximum-density-min", type=float, default=2.0)
    parser.add_argument("--maximum-pressure-min", type=float, default=0.1)
    parser.add_argument("--shock-radius-min", type=float, default=0.05)
    parser.add_argument("--shock-radius-max", type=float, default=0.45)
    parser.add_argument("--density-squared-reference", type=float)
    parser.add_argument("--pressure-integral-reference", type=float)
    parser.add_argument("--shock-radius-reference", type=float)
    parser.add_argument("--signature-relative-tolerance", type=float, default=2.0e-6)
    return parser.parse_args()


def relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0)


def main() -> int:
    args = parse_arguments()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != args.nx:
        raise AssertionError(f"Expected {args.nx} rows, found {len(rows)}")

    x = [float(row["x"]) for row in rows]
    rho = [float(row["rho"]) for row in rows]
    velocity = [float(row["u"]) for row in rows]
    pressure = [float(row["p"]) for row in rows]
    momentum = [float(row["momentum_density"]) for row in rows]
    energy = [float(row["total_energy_density"]) for row in rows]

    values = [value for row in rows for value in map(float, row.values())]
    if not all(math.isfinite(value) for value in values):
        raise AssertionError("Output contains a non-finite value")
    if min(rho) <= 0.0 or min(pressure) <= 0.0:
        raise AssertionError("Output contains non-positive density or pressure")

    dx = (args.x_max - args.x_min) / args.nx
    initial_mass = 0.0
    initial_momentum = 0.0
    initial_energy = 0.0
    for position in x:
        initial_pressure = (
            args.blast_pressure
            if abs(position - args.blast_center) <= args.blast_radius
            else args.ambient_pressure
        )
        initial_mass += args.ambient_density * dx
        initial_momentum += (
            args.ambient_density * args.initial_velocity * dx
        )
        initial_energy += (
            initial_pressure / (args.gamma - 1.0)
            + 0.5
            * args.ambient_density
            * args.initial_velocity
            * args.initial_velocity
        ) * dx

    final_mass = sum(rho) * dx
    final_momentum = sum(momentum) * dx
    final_energy = sum(energy) * dx

    mass_error = abs(final_mass - initial_mass)
    momentum_error = abs(final_momentum - initial_momentum)
    energy_error = abs(final_energy - initial_energy)

    density_symmetry_error = max(
        abs(left - right) for left, right in zip(rho, reversed(rho))
    )
    pressure_symmetry_error = max(
        abs(left - right) for left, right in zip(pressure, reversed(pressure))
    )
    velocity_antisymmetry_error = max(
        abs(left + right) for left, right in zip(velocity, reversed(velocity))
    )

    pressure_threshold = max(
        100.0 * args.ambient_pressure, 1.0e-4 * max(pressure)
    )
    active_radii = [
        abs(position - args.blast_center)
        for position, value in zip(x, pressure)
        if value > pressure_threshold
    ]
    if not active_radii:
        raise AssertionError("No blast region detected")
    shock_radius = max(active_radii)

    density_squared_integral = sum(value * value for value in rho) * dx
    pressure_integral = sum(pressure) * dx
    kinetic_energy_integral = sum(
        0.5 * density * speed * speed
        for density, speed in zip(rho, velocity)
    ) * dx

    metrics = {
        "minimum_density": min(rho),
        "maximum_density": max(rho),
        "minimum_pressure": min(pressure),
        "maximum_pressure": max(pressure),
        "shock_radius": shock_radius,
        "mass_balance_error": mass_error,
        "momentum_balance_error": momentum_error,
        "energy_balance_error": energy_error,
        "density_symmetry_error": density_symmetry_error,
        "pressure_symmetry_error": pressure_symmetry_error,
        "velocity_antisymmetry_error": velocity_antisymmetry_error,
        "density_squared_integral": density_squared_integral,
        "pressure_integral": pressure_integral,
        "kinetic_energy_integral": kinetic_energy_integral,
    }
    for name, value in metrics.items():
        print(f"{name}={value:.16e}")

    failures: list[str] = []
    for name, error in (
        ("mass_balance_error", mass_error),
        ("momentum_balance_error", momentum_error),
        ("energy_balance_error", energy_error),
    ):
        if error > args.balance_error_max:
            failures.append(
                f"{name}={error:.6e} exceeds {args.balance_error_max:.6e}"
            )

    for name, error in (
        ("density_symmetry_error", density_symmetry_error),
        ("pressure_symmetry_error", pressure_symmetry_error),
        ("velocity_antisymmetry_error", velocity_antisymmetry_error),
    ):
        if error > args.symmetry_error_max:
            failures.append(
                f"{name}={error:.6e} exceeds {args.symmetry_error_max:.6e}"
            )

    if max(rho) < args.maximum_density_min:
        failures.append(
            f"maximum_density={max(rho):.6e} is below "
            f"{args.maximum_density_min:.6e}"
        )
    if max(pressure) < args.maximum_pressure_min:
        failures.append(
            f"maximum_pressure={max(pressure):.6e} is below "
            f"{args.maximum_pressure_min:.6e}"
        )
    if not args.shock_radius_min <= shock_radius <= args.shock_radius_max:
        failures.append(
            f"shock_radius={shock_radius:.6e} is outside "
            f"[{args.shock_radius_min:.6e}, {args.shock_radius_max:.6e}]"
        )

    signatures = (
        ("density_squared_integral", density_squared_integral, args.density_squared_reference),
        ("pressure_integral", pressure_integral, args.pressure_integral_reference),
        ("shock_radius", shock_radius, args.shock_radius_reference),
    )
    for name, value, reference in signatures:
        if reference is None:
            continue
        error = relative_error(value, reference)
        if error > args.signature_relative_tolerance:
            failures.append(
                f"{name} relative error {error:.6e} exceeds "
                f"{args.signature_relative_tolerance:.6e}"
            )

    if failures:
        raise AssertionError("; ".join(failures))

    print("Sedov regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
