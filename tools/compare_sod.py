#!/usr/bin/env python3
"""Compare a PeleF Sod result against the exact ideal-gas Riemann solution."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PrimitiveState:
    rho: float
    velocity: float
    pressure: float


def pressure_function(
    pressure: float, state: PrimitiveState, gamma: float
) -> tuple[float, float]:
    sound_speed = math.sqrt(gamma * state.pressure / state.rho)
    if pressure > state.pressure:
        coefficient_a = 2.0 / ((gamma + 1.0) * state.rho)
        coefficient_b = (gamma - 1.0) / (gamma + 1.0) * state.pressure
        square_root = math.sqrt(coefficient_a / (pressure + coefficient_b))
        value = (pressure - state.pressure) * square_root
        derivative = square_root * (
            1.0 - 0.5 * (pressure - state.pressure) / (pressure + coefficient_b)
        )
        return value, derivative

    exponent = (gamma - 1.0) / (2.0 * gamma)
    pressure_ratio = pressure / state.pressure
    value = 2.0 * sound_speed / (gamma - 1.0) * (
        pressure_ratio**exponent - 1.0
    )
    derivative = (
        1.0
        / (state.rho * sound_speed)
        * pressure_ratio ** (-(gamma + 1.0) / (2.0 * gamma))
    )
    return value, derivative


def solve_star_region(
    left: PrimitiveState, right: PrimitiveState, gamma: float
) -> tuple[float, float]:
    left_sound_speed = math.sqrt(gamma * left.pressure / left.rho)
    right_sound_speed = math.sqrt(gamma * right.pressure / right.rho)
    pressure = max(
        1.0e-12,
        0.5 * (left.pressure + right.pressure)
        - 0.125
        * (right.velocity - left.velocity)
        * (left.rho + right.rho)
        * (left_sound_speed + right_sound_speed),
    )

    for _ in range(100):
        left_value, left_derivative = pressure_function(pressure, left, gamma)
        right_value, right_derivative = pressure_function(pressure, right, gamma)
        residual = (
            left_value
            + right_value
            + right.velocity
            - left.velocity
        )
        updated_pressure = pressure - residual / (left_derivative + right_derivative)
        updated_pressure = max(updated_pressure, 1.0e-12)
        if abs(updated_pressure - pressure) <= 1.0e-12 * (
            updated_pressure + pressure
        ):
            pressure = updated_pressure
            break
        pressure = updated_pressure
    else:
        raise RuntimeError("Exact Riemann pressure solve did not converge")

    left_value, _ = pressure_function(pressure, left, gamma)
    right_value, _ = pressure_function(pressure, right, gamma)
    velocity = 0.5 * (
        left.velocity + right.velocity + right_value - left_value
    )
    return pressure, velocity


def shock_star_density(
    state: PrimitiveState, star_pressure: float, gamma: float
) -> float:
    ratio = star_pressure / state.pressure
    gamma_ratio = (gamma - 1.0) / (gamma + 1.0)
    return state.rho * (ratio + gamma_ratio) / (gamma_ratio * ratio + 1.0)


def rarefaction_star_density(
    state: PrimitiveState, star_pressure: float, gamma: float
) -> float:
    return state.rho * (star_pressure / state.pressure) ** (1.0 / gamma)


def exact_state_at(
    similarity_coordinate: float,
    left: PrimitiveState,
    right: PrimitiveState,
    gamma: float,
    star_pressure: float,
    star_velocity: float,
) -> PrimitiveState:
    left_sound_speed = math.sqrt(gamma * left.pressure / left.rho)
    right_sound_speed = math.sqrt(gamma * right.pressure / right.rho)

    if similarity_coordinate <= star_velocity:
        if star_pressure > left.pressure:
            shock_speed = left.velocity - left_sound_speed * math.sqrt(
                (gamma + 1.0) / (2.0 * gamma) * star_pressure / left.pressure
                + (gamma - 1.0) / (2.0 * gamma)
            )
            if similarity_coordinate <= shock_speed:
                return left
            return PrimitiveState(
                shock_star_density(left, star_pressure, gamma),
                star_velocity,
                star_pressure,
            )

        star_sound_speed = left_sound_speed * (
            star_pressure / left.pressure
        ) ** ((gamma - 1.0) / (2.0 * gamma))
        head_speed = left.velocity - left_sound_speed
        tail_speed = star_velocity - star_sound_speed
        if similarity_coordinate <= head_speed:
            return left
        if similarity_coordinate >= tail_speed:
            return PrimitiveState(
                rarefaction_star_density(left, star_pressure, gamma),
                star_velocity,
                star_pressure,
            )

        velocity = 2.0 / (gamma + 1.0) * (
            left_sound_speed
            + 0.5 * (gamma - 1.0) * left.velocity
            + similarity_coordinate
        )
        sound_speed = 2.0 / (gamma + 1.0) * (
            left_sound_speed
            + 0.5 * (gamma - 1.0) * (left.velocity - similarity_coordinate)
        )
        ratio = sound_speed / left_sound_speed
        return PrimitiveState(
            left.rho * ratio ** (2.0 / (gamma - 1.0)),
            velocity,
            left.pressure * ratio ** (2.0 * gamma / (gamma - 1.0)),
        )

    if star_pressure > right.pressure:
        shock_speed = right.velocity + right_sound_speed * math.sqrt(
            (gamma + 1.0) / (2.0 * gamma) * star_pressure / right.pressure
            + (gamma - 1.0) / (2.0 * gamma)
        )
        if similarity_coordinate >= shock_speed:
            return right
        return PrimitiveState(
            shock_star_density(right, star_pressure, gamma),
            star_velocity,
            star_pressure,
        )

    star_sound_speed = right_sound_speed * (
        star_pressure / right.pressure
    ) ** ((gamma - 1.0) / (2.0 * gamma))
    head_speed = right.velocity + right_sound_speed
    tail_speed = star_velocity + star_sound_speed
    if similarity_coordinate >= head_speed:
        return right
    if similarity_coordinate <= tail_speed:
        return PrimitiveState(
            rarefaction_star_density(right, star_pressure, gamma),
            star_velocity,
            star_pressure,
        )

    velocity = 2.0 / (gamma + 1.0) * (
        -right_sound_speed
        + 0.5 * (gamma - 1.0) * right.velocity
        + similarity_coordinate
    )
    sound_speed = 2.0 / (gamma + 1.0) * (
        right_sound_speed
        + 0.5 * (gamma - 1.0) * (similarity_coordinate - right.velocity)
    )
    ratio = sound_speed / right_sound_speed
    return PrimitiveState(
        right.rho * ratio ** (2.0 / (gamma - 1.0)),
        velocity,
        right.pressure * ratio ** (2.0 * gamma / (gamma - 1.0)),
    )


def l1_error(numerical: Iterable[float], exact: Iterable[float]) -> float:
    differences = [
        abs(numerical_value - exact_value)
        for numerical_value, exact_value in zip(numerical, exact)
    ]
    if not differences:
        raise ValueError("No values supplied to error calculation")
    return sum(differences) / len(differences)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("sod.csv"))
    parser.add_argument("--nx", type=int, default=400)
    parser.add_argument("--time", type=float, default=0.2)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=1.0)
    parser.add_argument("--discontinuity", type=float, default=0.5)
    parser.add_argument("--density-l1-max", type=float, default=3.0e-2)
    parser.add_argument("--pressure-l1-max", type=float, default=3.0e-2)
    parser.add_argument("--mass-error-max", type=float, default=2.0e-12)
    parser.add_argument("--energy-error-max", type=float, default=2.0e-12)
    parser.add_argument("--momentum-error-max", type=float, default=2.0e-12)
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    if args.time <= 0.0:
        raise ValueError("Comparison time must be positive")

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != args.nx:
        raise AssertionError(f"Expected {args.nx} rows, found {len(rows)}")

    x = [float(row["x"]) for row in rows]
    rho = [float(row["rho"]) for row in rows]
    velocity = [float(row["u"]) for row in rows]
    pressure = [float(row["p"]) for row in rows]
    total_energy_density = [float(row["total_energy_density"]) for row in rows]
    momentum_density = [float(row["momentum_density"]) for row in rows]

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

    density_error = l1_error(rho, (state.rho for state in exact))
    pressure_error = l1_error(pressure, (state.pressure for state in exact))

    dx = (args.x_max - args.x_min) / args.nx
    mass = sum(rho) * dx
    momentum = sum(momentum_density) * dx
    energy = sum(total_energy_density) * dx

    initial_mass = (
        left.rho * (args.discontinuity - args.x_min)
        + right.rho * (args.x_max - args.discontinuity)
    )
    initial_energy = (
        left.pressure
        / (args.gamma - 1.0)
        * (args.discontinuity - args.x_min)
        + right.pressure
        / (args.gamma - 1.0)
        * (args.x_max - args.discontinuity)
    )
    expected_momentum = (left.pressure - right.pressure) * args.time

    mass_error = abs(mass - initial_mass)
    energy_error = abs(energy - initial_energy)
    momentum_error = abs(momentum - expected_momentum)

    metrics = {
        "density_l1": density_error,
        "pressure_l1": pressure_error,
        "mass_error": mass_error,
        "energy_error": energy_error,
        "momentum_error": momentum_error,
        "minimum_density": min(rho),
        "minimum_pressure": min(pressure),
        "star_pressure": star_pressure,
        "star_velocity": star_velocity,
    }
    for name, value in metrics.items():
        print(f"{name}={value:.16e}")

    limits = {
        "density_l1": args.density_l1_max,
        "pressure_l1": args.pressure_l1_max,
        "mass_error": args.mass_error_max,
        "energy_error": args.energy_error_max,
        "momentum_error": args.momentum_error_max,
    }
    failures = [
        f"{name}={metrics[name]:.6e} exceeds {limit:.6e}"
        for name, limit in limits.items()
        if metrics[name] > limit
    ]
    if failures:
        raise AssertionError("; ".join(failures))

    print("Sod regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
