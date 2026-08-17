#!/usr/bin/env python3
"""Verify the canonical one-dimensional Shu-Osher regression."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--nx", type=int, default=800)
    parser.add_argument("--time", type=float, default=1.8)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--x-min", type=float, default=-5.0)
    parser.add_argument("--x-max", type=float, default=5.0)
    parser.add_argument("--shock-location", type=float, default=-4.0)
    parser.add_argument("--left-density", type=float, default=3.857143)
    parser.add_argument("--left-velocity", type=float, default=2.629369)
    parser.add_argument("--left-pressure", type=float, default=10.33333)
    parser.add_argument("--density-base", type=float, default=1.0)
    parser.add_argument("--density-amplitude", type=float, default=0.2)
    parser.add_argument("--density-wavenumber", type=float, default=5.0)
    parser.add_argument("--right-velocity", type=float, default=0.0)
    parser.add_argument("--right-pressure", type=float, default=1.0)
    parser.add_argument("--balance-error-max", type=float, default=2.0e-10)
    parser.add_argument("--maximum-density-min", type=float, default=4.4)
    parser.add_argument("--window-range-min", type=float, default=3.0)
    parser.add_argument("--window-extrema-min", type=int, default=15)
    parser.add_argument(
        "--density-squared-reference", type=float, default=112.75925276923157
    )
    parser.add_argument(
        "--density-moment-reference", type=float, default=-7.070225172008236
    )
    parser.add_argument("--signature-relative-tolerance", type=float, default=2.0e-6)
    return parser.parse_args()


def physical_flux(
    rho: float,
    velocity: float,
    pressure: float,
    gamma: float,
) -> tuple[float, float, float]:
    energy = pressure / (gamma - 1.0) + 0.5 * rho * velocity * velocity
    return (
        rho * velocity,
        rho * velocity * velocity + pressure,
        (energy + pressure) * velocity,
    )


def relative_error(value: float, reference: float) -> float:
    scale = max(abs(reference), 1.0)
    return abs(value - reference) / scale


def main() -> int:
    args = parse_arguments()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != args.nx:
        raise AssertionError(f"Expected {args.nx} rows, found {len(rows)}")

    x = [float(row["x"]) for row in rows]
    rho = [float(row["rho"]) for row in rows]
    pressure = [float(row["p"]) for row in rows]
    momentum = [float(row["momentum_density"]) for row in rows]
    energy = [float(row["total_energy_density"]) for row in rows]

    values = [
        value
        for row in rows
        for value in map(float, row.values())
    ]
    if not all(math.isfinite(value) for value in values):
        raise AssertionError("Output contains a non-finite value")
    if min(rho) <= 0.0 or min(pressure) <= 0.0:
        raise AssertionError("Output contains non-positive density or pressure")

    dx = (args.x_max - args.x_min) / args.nx
    initial_mass = 0.0
    initial_momentum = 0.0
    initial_energy = 0.0
    for position in x:
        if position < args.shock_location:
            initial_rho = args.left_density
            initial_velocity = args.left_velocity
            initial_pressure = args.left_pressure
        else:
            initial_rho = args.density_base + args.density_amplitude * math.sin(
                args.density_wavenumber * position
            )
            initial_velocity = args.right_velocity
            initial_pressure = args.right_pressure
        initial_mass += initial_rho * dx
        initial_momentum += initial_rho * initial_velocity * dx
        initial_energy += (
            initial_pressure / (args.gamma - 1.0)
            + 0.5 * initial_rho * initial_velocity * initial_velocity
        ) * dx

    left_mass_flux, left_momentum_flux, left_energy_flux = physical_flux(
        args.left_density,
        args.left_velocity,
        args.left_pressure,
        args.gamma,
    )
    right_mass_flux, right_momentum_flux, right_energy_flux = physical_flux(
        args.density_base
        + args.density_amplitude
        * math.sin(args.density_wavenumber * (args.x_max - 0.5 * dx)),
        args.right_velocity,
        args.right_pressure,
        args.gamma,
    )

    expected_mass = initial_mass + args.time * (left_mass_flux - right_mass_flux)
    expected_momentum = initial_momentum + args.time * (
        left_momentum_flux - right_momentum_flux
    )
    expected_energy = initial_energy + args.time * (
        left_energy_flux - right_energy_flux
    )

    final_mass = sum(rho) * dx
    final_momentum = sum(momentum) * dx
    final_energy = sum(energy) * dx
    mass_error = abs(final_mass - expected_mass)
    momentum_error = abs(final_momentum - expected_momentum)
    energy_error = abs(final_energy - expected_energy)

    window_indices = [
        index
        for index, position in enumerate(x)
        if -2.0 <= position <= 2.5
    ]
    window_density = [rho[index] for index in window_indices]
    density_range = max(window_density) - min(window_density)

    extrema_count = 0
    for index in window_indices:
        if index == 0 or index == len(rho) - 1:
            continue
        left_difference = rho[index] - rho[index - 1]
        right_difference = rho[index + 1] - rho[index]
        if left_difference * right_difference < 0.0:
            extrema_count += 1

    density_squared_integral = sum(value * value for value in rho) * dx
    density_moment = sum(
        value * math.sin(0.37 * position)
        for value, position in zip(rho, x)
    ) * dx

    metrics = {
        "minimum_density": min(rho),
        "maximum_density": max(rho),
        "minimum_pressure": min(pressure),
        "window_density_range": density_range,
        "window_extrema": float(extrema_count),
        "mass_balance_error": mass_error,
        "momentum_balance_error": momentum_error,
        "energy_balance_error": energy_error,
        "density_squared_integral": density_squared_integral,
        "density_weighted_moment": density_moment,
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

    if max(rho) < args.maximum_density_min:
        failures.append(
            f"maximum_density={max(rho):.6e} is below "
            f"{args.maximum_density_min:.6e}"
        )
    if density_range < args.window_range_min:
        failures.append(
            f"window_density_range={density_range:.6e} is below "
            f"{args.window_range_min:.6e}"
        )
    if extrema_count < args.window_extrema_min:
        failures.append(
            f"window_extrema={extrema_count} is below "
            f"{args.window_extrema_min}"
        )

    squared_error = relative_error(
        density_squared_integral, args.density_squared_reference
    )
    moment_error = relative_error(
        density_moment, args.density_moment_reference
    )
    if squared_error > args.signature_relative_tolerance:
        failures.append(
            "density_squared_integral relative error "
            f"{squared_error:.6e} exceeds "
            f"{args.signature_relative_tolerance:.6e}"
        )
    if moment_error > args.signature_relative_tolerance:
        failures.append(
            "density_weighted_moment relative error "
            f"{moment_error:.6e} exceeds "
            f"{args.signature_relative_tolerance:.6e}"
        )

    if failures:
        raise AssertionError("; ".join(failures))

    print("Shu-Osher regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
