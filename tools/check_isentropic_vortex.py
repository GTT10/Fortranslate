#!/usr/bin/env python3
"""Verify the periodic two-dimensional isentropic-vortex regression."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--ny", type=int, default=64)
    parser.add_argument("--time", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=10.0)
    parser.add_argument("--y-min", type=float, default=0.0)
    parser.add_argument("--y-max", type=float, default=10.0)
    parser.add_argument("--center-x", type=float, default=5.0)
    parser.add_argument("--center-y", type=float, default=5.0)
    parser.add_argument("--strength", type=float, default=5.0)
    parser.add_argument("--base-density", type=float, default=1.0)
    parser.add_argument("--base-pressure", type=float, default=1.0)
    parser.add_argument("--base-velocity-x", type=float, default=1.0)
    parser.add_argument("--base-velocity-y", type=float, default=1.0)
    parser.add_argument("--density-l1-max", type=float, default=5.0e-4)
    parser.add_argument("--pressure-l1-max", type=float, default=7.0e-4)
    parser.add_argument("--velocity-l1-max", type=float, default=1.2e-3)
    parser.add_argument("--conservation-error-max", type=float, default=2.0e-10)
    return parser.parse_args()


def periodic_displacement(delta: float, period: float) -> float:
    return (delta + 0.5 * period) % period - 0.5 * period


def exact_primitive(
    x: float,
    y: float,
    time: float,
    args: argparse.Namespace,
) -> tuple[float, float, float, float]:
    length_x = args.x_max - args.x_min
    length_y = args.y_max - args.y_min
    center_x = args.x_min + (
        args.center_x + args.base_velocity_x * time - args.x_min
    ) % length_x
    center_y = args.y_min + (
        args.center_y + args.base_velocity_y * time - args.y_min
    ) % length_y

    displacement_x = periodic_displacement(x - center_x, length_x)
    displacement_y = periodic_displacement(y - center_y, length_y)
    radius_squared = displacement_x**2 + displacement_y**2

    velocity_factor = args.strength / (2.0 * math.pi) * math.exp(
        0.5 * (1.0 - radius_squared)
    )
    temperature_perturbation = -(
        (args.gamma - 1.0)
        * args.strength**2
        / (8.0 * args.gamma * math.pi**2)
        * math.exp(1.0 - radius_squared)
    )
    base_temperature = args.base_pressure / args.base_density
    temperature_ratio = (base_temperature + temperature_perturbation) / base_temperature

    density = args.base_density * temperature_ratio ** (1.0 / (args.gamma - 1.0))
    velocity_x = args.base_velocity_x - velocity_factor * displacement_y
    velocity_y = args.base_velocity_y + velocity_factor * displacement_x
    pressure = args.base_pressure * temperature_ratio ** (
        args.gamma / (args.gamma - 1.0)
    )
    return density, velocity_x, velocity_y, pressure


def conserved_from_primitive(
    density: float,
    velocity_x: float,
    velocity_y: float,
    pressure: float,
    gamma: float,
) -> tuple[float, float, float, float, float]:
    energy = pressure / (gamma - 1.0) + 0.5 * density * (
        velocity_x**2 + velocity_y**2
    )
    return (
        density,
        density * velocity_x,
        density * velocity_y,
        0.0,
        energy,
    )


def main() -> int:
    args = parse_arguments()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    expected_rows = args.nx * args.ny
    if len(rows) != expected_rows:
        raise AssertionError(f"Expected {expected_rows} rows, found {len(rows)}")

    finite_values = [
        value
        for row in rows
        for value in map(float, row.values())
    ]
    if not all(math.isfinite(value) for value in finite_values):
        raise AssertionError("Output contains a non-finite value")

    dx = (args.x_max - args.x_min) / args.nx
    dy = (args.y_max - args.y_min) / args.ny
    cell_area = dx * dy

    errors = {"rho": 0.0, "u": 0.0, "v": 0.0, "p": 0.0}
    initial_totals = [0.0] * 5
    final_totals = [0.0] * 5
    minimum_density = math.inf
    minimum_pressure = math.inf

    for row in rows:
        x = float(row["x"])
        y = float(row["y"])
        numerical = {
            "rho": float(row["rho"]),
            "u": float(row["u"]),
            "v": float(row["v"]),
            "p": float(row["p"]),
        }
        exact = exact_primitive(x, y, args.time, args)
        exact_mapping = dict(zip(("rho", "u", "v", "p"), exact))
        for name in errors:
            errors[name] += abs(numerical[name] - exact_mapping[name])

        minimum_density = min(minimum_density, numerical["rho"])
        minimum_pressure = min(minimum_pressure, numerical["p"])

        initial_primitive = exact_primitive(x, y, 0.0, args)
        initial_conserved = conserved_from_primitive(*initial_primitive, args.gamma)
        final_conserved = (
            numerical["rho"],
            float(row["momentum_x_density"]),
            float(row["momentum_y_density"]),
            0.0,
            float(row["total_energy_density"]),
        )
        for component in range(5):
            initial_totals[component] += initial_conserved[component] * cell_area
            final_totals[component] += final_conserved[component] * cell_area

    for name in errors:
        errors[name] /= expected_rows

    conservation_error = max(
        abs(final - initial)
        for final, initial in zip(final_totals, initial_totals)
    )

    metrics = {
        "density_l1": errors["rho"],
        "velocity_x_l1": errors["u"],
        "velocity_y_l1": errors["v"],
        "pressure_l1": errors["p"],
        "minimum_density": minimum_density,
        "minimum_pressure": minimum_pressure,
        "conservation_error": conservation_error,
    }
    for name, value in metrics.items():
        print(f"{name}={value:.16e}")

    failures: list[str] = []
    if errors["rho"] > args.density_l1_max:
        failures.append("density L1 error exceeds threshold")
    if errors["p"] > args.pressure_l1_max:
        failures.append("pressure L1 error exceeds threshold")
    if max(errors["u"], errors["v"]) > args.velocity_l1_max:
        failures.append("velocity L1 error exceeds threshold")
    if conservation_error > args.conservation_error_max:
        failures.append("periodic conservation error exceeds threshold")
    if minimum_density <= 0.0 or minimum_pressure <= 0.0:
        failures.append("non-positive density or pressure")

    if failures:
        raise AssertionError("; ".join(failures))

    print("Isentropic-vortex regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
