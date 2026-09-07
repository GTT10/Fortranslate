#!/usr/bin/env python3
"""Validate selected-mechanism composite AMR reactive 1D CSV output."""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


BASE_COLUMNS = (
    "level",
    "cell_dx",
    "time",
    "x",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
)


def close(left: float, right: float, scale: float) -> bool:
    return abs(left - right) <= 5.0e-12 * max(scale, abs(left), abs(right), 1.0e-30)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--expected-levels", nargs="+", type=int, required=True)
    parser.add_argument("--x-lower", type=float, default=0.0)
    parser.add_argument("--x-upper", type=float, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--refinement-ratio", type=int, required=True)
    parser.add_argument("--minimum-finest-regions", type=int, default=1)
    parser.add_argument("--require-finest-touches-lower", action="store_true")
    parser.add_argument("--require-finest-touches-upper", action="store_true")
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    args = parser.parse_args()

    expected_levels = sorted(set(args.expected_levels))
    if expected_levels != list(range(len(expected_levels))):
        raise AssertionError("expected AMR levels must be consecutive from zero")
    if args.x_upper <= args.x_lower:
        raise AssertionError("invalid AMR domain")
    if args.refinement_ratio < 2:
        raise AssertionError("AMR refinement ratio must be at least two")
    if args.minimum_finest_regions < 1:
        raise AssertionError("minimum finest-region count must be positive")

    with args.input.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected_columns = [*BASE_COLUMNS, *(f"Y_{name}" for name in args.species)]
        if reader.fieldnames != expected_columns:
            raise AssertionError(
                f"selected AMR reactive 1D columns differ: {reader.fieldnames!r}"
            )
        rows = [{name: float(value) for name, value in row.items()} for row in reader]

    if not rows:
        raise AssertionError("selected AMR reactive 1D output is empty")
    domain_length = args.x_upper - args.x_lower
    levels: dict[int, list[dict[str, float]]] = defaultdict(list)
    minimum_density = math.inf
    minimum_pressure = math.inf
    minimum_temperature = math.inf
    minimum_fraction = math.inf
    closure_error = 0.0
    time_error = 0.0
    for row in rows:
        if not all(math.isfinite(value) for value in row.values()):
            raise AssertionError("selected AMR reactive 1D output contains nonfinite data")
        level = round(row["level"])
        if row["level"] != level or level < 0:
            raise AssertionError(f"invalid AMR level value: {row['level']}")
        if row["cell_dx"] <= 0.0:
            raise AssertionError("selected AMR reactive 1D output has nonpositive cell_dx")
        levels[level].append(row)
        time_error = max(time_error, abs(row["time"] - args.final_time))
        minimum_density = min(minimum_density, row["rho"])
        minimum_pressure = min(minimum_pressure, row["pressure"])
        minimum_temperature = min(minimum_temperature, row["temperature"])
        fractions = [row[f"Y_{name}"] for name in args.species]
        minimum_fraction = min(minimum_fraction, *fractions)
        closure_error = max(closure_error, abs(sum(fractions) - 1.0))

    actual_levels = sorted(levels)
    if actual_levels != expected_levels:
        raise AssertionError(
            f"expected AMR levels {expected_levels}, found {actual_levels}"
        )
    if time_error > 5.0e-13 * max(abs(args.final_time), 1.0e-30):
        raise AssertionError(f"selected AMR reactive 1D final-time error is {time_error}")
    if min(minimum_density, minimum_pressure, minimum_temperature) <= 0.0:
        raise AssertionError("selected AMR reactive 1D output lost thermodynamic positivity")
    if minimum_fraction < -2.0e-12 or closure_error > 2.0e-10:
        raise AssertionError("selected AMR reactive 1D composition is invalid")

    dx_by_level: dict[int, float] = {}
    for level in actual_levels:
        level_rows = levels[level]
        coordinates = [row["x"] for row in level_rows]
        if coordinates != sorted(coordinates) or len(set(coordinates)) != len(coordinates):
            raise AssertionError(f"AMR level {level} coordinates are not strictly ordered")
        dx = level_rows[0]["cell_dx"]
        if any(not close(row["cell_dx"], dx, domain_length) for row in level_rows):
            raise AssertionError(f"AMR level {level} has inconsistent cell widths")
        if any(
            right - left < dx - 5.0e-12 * domain_length
            for left, right in zip(coordinates, coordinates[1:])
        ):
            raise AssertionError(f"AMR level {level} contains overlapping cells")
        dx_by_level[level] = dx
    for level in actual_levels[1:]:
        measured_ratio = dx_by_level[level - 1] / dx_by_level[level]
        if not close(measured_ratio, float(args.refinement_ratio), 1.0):
            raise AssertionError(
                f"AMR level {level} refinement ratio is {measured_ratio}"
            )

    intervals = sorted(
        (
            row["x"] - 0.5 * row["cell_dx"],
            row["x"] + 0.5 * row["cell_dx"],
            round(row["level"]),
        )
        for row in rows
    )
    coverage_error = max(
        abs(intervals[0][0] - args.x_lower),
        abs(intervals[-1][1] - args.x_upper),
        abs(sum(right - left for left, right, _ in intervals) - domain_length),
    )
    for (_, left_right, _), (right_left, _, _) in zip(intervals, intervals[1:]):
        coverage_error = max(coverage_error, abs(right_left - left_right))
    if coverage_error > 5.0e-12 * domain_length:
        raise AssertionError(f"composite AMR coverage error is {coverage_error}")

    finest_level = actual_levels[-1]
    finest_dx = dx_by_level[finest_level]
    finest_coordinates = [row["x"] for row in levels[finest_level]]
    finest_lower = min(finest_coordinates) - 0.5 * finest_dx
    finest_upper = max(finest_coordinates) + 0.5 * finest_dx
    if args.require_finest_touches_lower and not close(
        finest_lower, args.x_lower, domain_length
    ):
        raise AssertionError(
            f"finest AMR level does not touch lower boundary: {finest_lower}"
        )
    if args.require_finest_touches_upper and not close(
        finest_upper, args.x_upper, domain_length
    ):
        raise AssertionError(
            f"finest AMR level does not touch upper boundary: {finest_upper}"
        )
    finest_regions = 1 + sum(
        right - left > 1.5 * finest_dx
        for left, right in zip(finest_coordinates, finest_coordinates[1:])
    )
    if finest_regions < args.minimum_finest_regions:
        raise AssertionError(
            f"expected at least {args.minimum_finest_regions} finest regions, "
            f"found {finest_regions}"
        )

    activity = 0.0
    activity_options = (args.activity_species, args.initial_mass_fraction)
    if any(value is not None for value in activity_options):
        if any(value is None for value in activity_options):
            raise AssertionError("both activity options are required")
        column = f"Y_{args.activity_species}"
        if column not in rows[0]:
            raise AssertionError(f"activity species is absent: {args.activity_species}")
        activity = max(abs(row[column] - args.initial_mass_fraction) for row in rows)
        if activity < args.minimum_change:
            raise AssertionError(
                f"selected AMR reactive 1D chemistry activity is only {activity}"
            )

    level_summary = ",".join(f"{level}:{len(levels[level])}" for level in actual_levels)
    print(f"levels={level_summary}")
    print(f"finest_regions={finest_regions}")
    print(f"finest_bounds={finest_lower:.16e}:{finest_upper:.16e}")
    print(f"minimum_density={minimum_density:.16e}")
    print(f"minimum_pressure={minimum_pressure:.16e}")
    print(f"minimum_temperature={minimum_temperature:.16e}")
    print(f"minimum_mass_fraction={minimum_fraction:.16e}")
    print(f"maximum_closure_error={closure_error:.16e}")
    print(f"maximum_coverage_error={coverage_error:.16e}")
    print(f"maximum_requested_species_change={activity:.16e}")
    print("selected AMR reactive 1D regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
