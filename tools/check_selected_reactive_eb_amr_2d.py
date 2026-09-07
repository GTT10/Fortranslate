#!/usr/bin/env python3
"""Validate selected static two-level reactive EB AMR CSV output."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


BASE_COLUMNS = [
    "time",
    "x",
    "y",
    "volume_fraction",
    "cell_type",
    "boundary_length",
    "boundary_normal_x",
    "boundary_normal_y",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine", required=True, type=Path)
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--coarse-nx", required=True, type=int)
    parser.add_argument("--coarse-ny", required=True, type=int)
    parser.add_argument("--x-lower", required=True, type=float)
    parser.add_argument("--x-upper", required=True, type=float)
    parser.add_argument("--y-lower", required=True, type=float)
    parser.add_argument("--y-upper", required=True, type=float)
    parser.add_argument("--coarse-i-lower", required=True, type=int)
    parser.add_argument("--coarse-i-upper", required=True, type=int)
    parser.add_argument("--coarse-j-lower", required=True, type=int)
    parser.add_argument("--coarse-j-upper", required=True, type=int)
    parser.add_argument("--refinement-ratio", required=True, type=int)
    parser.add_argument("--final-time", required=True, type=float)
    parser.add_argument("--coarse-counts", required=True, nargs=3, type=int)
    parser.add_argument("--fine-counts", required=True, nargs=3, type=int)
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    parser.add_argument("--expected-coarse-sha256")
    parser.add_argument("--expected-fine-sha256")
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(message)


def read_rows(path: Path, expected_header: list[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_header:
            fail(f"{path.name}: selected EB AMR CSV schema mismatch")
        return list(reader)


def check_level(
    rows: list[dict[str, str]],
    expected_header: list[str],
    species: list[str],
    nx: int,
    ny: int,
    x_lower: float,
    x_upper: float,
    y_lower: float,
    y_upper: float,
    final_time: float,
    expected_counts: list[int],
    activity_species: str | None,
    initial_mass_fraction: float | None,
) -> tuple[dict[int, int], float, float]:
    if len(rows) != nx * ny:
        fail("selected EB AMR CSV row count mismatch")
    if nx < 1 or ny < 1 or x_upper <= x_lower or y_upper <= y_lower:
        fail("selected EB AMR expected level geometry is invalid")

    dx = (x_upper - x_lower) / nx
    dy = (y_upper - y_lower) / ny
    coordinate_tolerance = 5.0e-15
    type_counts = {0: 0, 1: 0, 2: 0}
    maximum_closure_error = 0.0
    maximum_activity_change = 0.0
    for index, row in enumerate(rows):
        try:
            values = [float(row[column]) for column in expected_header]
        except ValueError as error:
            fail(f"selected EB AMR CSV contains nonnumeric data: {error}")
        if not all(math.isfinite(value) for value in values):
            fail("selected EB AMR CSV contains nonfinite data")
        if not math.isclose(
            float(row["time"]), final_time, rel_tol=0.0, abs_tol=1.0e-14
        ):
            fail("selected EB AMR final time mismatch")

        expected_x = x_lower + (index % nx + 0.5) * dx
        expected_y = y_lower + (index // nx + 0.5) * dy
        if abs(float(row["x"]) - expected_x) > coordinate_tolerance:
            fail("selected EB AMR x-coordinate order mismatch")
        if abs(float(row["y"]) - expected_y) > coordinate_tolerance:
            fail("selected EB AMR y-coordinate order mismatch")

        cell_type_value = float(row["cell_type"])
        cell_type = int(cell_type_value)
        if cell_type_value != cell_type or cell_type not in type_counts:
            fail("selected EB AMR cell type is invalid")
        type_counts[cell_type] += 1
        volume_fraction = float(row["volume_fraction"])
        if not 0.0 <= volume_fraction <= 1.0:
            fail("selected EB AMR volume fraction is outside [0,1]")
        if cell_type == 0 and volume_fraction != 0.0:
            fail("selected EB AMR covered cell has nonzero volume")
        if cell_type == 1 and not 0.0 < volume_fraction < 1.0:
            fail("selected EB AMR cut cell lacks a partial volume")
        if cell_type == 2 and volume_fraction != 1.0:
            fail("selected EB AMR regular cell lacks unit volume")
        if float(row["boundary_length"]) < 0.0:
            fail("selected EB AMR boundary length is negative")
        if (
            float(row["rho"]) <= 0.0
            or float(row["pressure"]) <= 0.0
            or float(row["temperature"]) <= 0.0
        ):
            fail("selected EB AMR state is not positive")

        closure = sum(float(row[f"Y_{name}"]) for name in species)
        maximum_closure_error = max(
            maximum_closure_error, abs(closure - 1.0)
        )
        if any(float(row[f"Y_{name}"]) < -1.0e-13 for name in species):
            fail("selected EB AMR mass fraction is negative")
        if activity_species and cell_type != 0:
            assert initial_mass_fraction is not None
            maximum_activity_change = max(
                maximum_activity_change,
                abs(
                    float(row[f"Y_{activity_species}"])
                    - initial_mass_fraction
                ),
            )

    expected_by_type = {
        2: expected_counts[0],
        1: expected_counts[1],
        0: expected_counts[2],
    }
    if type_counts != expected_by_type:
        fail(
            "selected EB AMR cell counts mismatch: "
            f"observed={type_counts}, expected={expected_by_type}"
        )
    if maximum_closure_error > 5.0e-12:
        fail("selected EB AMR composition closure is too large")
    return type_counts, maximum_closure_error, maximum_activity_change


def conserved_values(
    row: dict[str, str], species: list[str]
) -> list[float]:
    density = float(row["rho"])
    return [
        density,
        density * float(row["u"]),
        density * float(row["v"]),
        density * float(row["w"]),
        float(row["rhoE"]),
        *(density * float(row[f"Y_{name}"]) for name in species),
    ]


def check_average_down(
    coarse_rows: list[dict[str, str]],
    fine_rows: list[dict[str, str]],
    species: list[str],
    coarse_nx: int,
    fine_nx: int,
    i_lower: int,
    i_upper: int,
    j_lower: int,
    j_upper: int,
    ratio: int,
) -> tuple[float, float]:
    maximum_geometry_error = 0.0
    maximum_state_error = 0.0
    scale = ratio * ratio
    for parent_j, coarse_j in enumerate(range(j_lower - 1, j_upper)):
        for parent_i, coarse_i in enumerate(range(i_lower - 1, i_upper)):
            coarse = coarse_rows[coarse_j * coarse_nx + coarse_i]
            coarse_fraction = float(coarse["volume_fraction"])
            coarse_values = conserved_values(coarse, species)
            fine_fraction_sum = 0.0
            fine_sums = [0.0] * len(coarse_values)
            for fine_j in range(parent_j * ratio, (parent_j + 1) * ratio):
                for fine_i in range(
                    parent_i * ratio, (parent_i + 1) * ratio
                ):
                    fine = fine_rows[fine_j * fine_nx + fine_i]
                    fraction = float(fine["volume_fraction"])
                    fine_fraction_sum += fraction
                    for component, value in enumerate(
                        conserved_values(fine, species)
                    ):
                        fine_sums[component] += fraction * value
            averaged_fraction = fine_fraction_sum / scale
            maximum_geometry_error = max(
                maximum_geometry_error,
                abs(coarse_fraction - averaged_fraction),
            )
            for coarse_value, fine_sum in zip(
                coarse_values, fine_sums, strict=True
            ):
                coarse_weighted = coarse_fraction * coarse_value
                fine_weighted = fine_sum / scale
                relative_error = abs(coarse_weighted - fine_weighted) / max(
                    1.0, abs(coarse_weighted), abs(fine_weighted)
                )
                maximum_state_error = max(
                    maximum_state_error, relative_error
                )
    if maximum_geometry_error > 5.0e-13:
        fail("selected EB AMR coarse/fine geometry is inconsistent")
    if maximum_state_error > 5.0e-12:
        fail("selected EB AMR average-down state is inconsistent")
    return maximum_geometry_error, maximum_state_error


def main() -> None:
    arguments = parse_arguments()
    if arguments.refinement_ratio < 2:
        fail("selected EB AMR refinement ratio must be at least two")
    if arguments.activity_species:
        if arguments.activity_species not in arguments.species:
            fail("--activity-species must be present in --species")
        if arguments.initial_mass_fraction is None:
            fail("--activity-species requires --initial-mass-fraction")
    if bool(arguments.expected_coarse_sha256) != bool(
        arguments.expected_fine_sha256
    ):
        fail("coarse and fine SHA-256 expectations must be supplied together")

    patch_nx = arguments.coarse_i_upper - arguments.coarse_i_lower + 1
    patch_ny = arguments.coarse_j_upper - arguments.coarse_j_lower + 1
    if (
        patch_nx < 1
        or patch_ny < 1
        or arguments.coarse_i_lower < 1
        or arguments.coarse_j_lower < 1
        or arguments.coarse_i_upper > arguments.coarse_nx
        or arguments.coarse_j_upper > arguments.coarse_ny
    ):
        fail("selected EB AMR coarse patch is invalid")
    fine_nx = patch_nx * arguments.refinement_ratio
    fine_ny = patch_ny * arguments.refinement_ratio
    coarse_dx = (
        arguments.x_upper - arguments.x_lower
    ) / arguments.coarse_nx
    coarse_dy = (
        arguments.y_upper - arguments.y_lower
    ) / arguments.coarse_ny
    fine_x_lower = arguments.x_lower + (
        arguments.coarse_i_lower - 1
    ) * coarse_dx
    fine_x_upper = arguments.x_lower + arguments.coarse_i_upper * coarse_dx
    fine_y_lower = arguments.y_lower + (
        arguments.coarse_j_lower - 1
    ) * coarse_dy
    fine_y_upper = arguments.y_lower + arguments.coarse_j_upper * coarse_dy

    expected_header = BASE_COLUMNS + [
        f"Y_{name}" for name in arguments.species
    ]
    coarse_rows = read_rows(arguments.coarse, expected_header)
    fine_rows = read_rows(arguments.fine, expected_header)
    coarse_result = check_level(
        coarse_rows,
        expected_header,
        arguments.species,
        arguments.coarse_nx,
        arguments.coarse_ny,
        arguments.x_lower,
        arguments.x_upper,
        arguments.y_lower,
        arguments.y_upper,
        arguments.final_time,
        arguments.coarse_counts,
        arguments.activity_species,
        arguments.initial_mass_fraction,
    )
    fine_result = check_level(
        fine_rows,
        expected_header,
        arguments.species,
        fine_nx,
        fine_ny,
        fine_x_lower,
        fine_x_upper,
        fine_y_lower,
        fine_y_upper,
        arguments.final_time,
        arguments.fine_counts,
        arguments.activity_species,
        arguments.initial_mass_fraction,
    )
    geometry_error, state_error = check_average_down(
        coarse_rows,
        fine_rows,
        arguments.species,
        arguments.coarse_nx,
        fine_nx,
        arguments.coarse_i_lower,
        arguments.coarse_i_upper,
        arguments.coarse_j_lower,
        arguments.coarse_j_upper,
        arguments.refinement_ratio,
    )
    maximum_activity_change = max(coarse_result[2], fine_result[2])
    if (
        arguments.activity_species
        and maximum_activity_change < arguments.minimum_change
    ):
        fail("selected EB AMR activity is below the requested gate")
    if arguments.expected_coarse_sha256:
        coarse_digest = hashlib.sha256(arguments.coarse.read_bytes()).hexdigest()
        fine_digest = hashlib.sha256(arguments.fine.read_bytes()).hexdigest()
        if coarse_digest != arguments.expected_coarse_sha256:
            fail("selected EB AMR coarse SHA-256 mismatch")
        if fine_digest != arguments.expected_fine_sha256:
            fail("selected EB AMR fine SHA-256 mismatch")

    print(f"coarse_cell_counts={coarse_result[0]}")
    print(f"fine_cell_counts={fine_result[0]}")
    print(
        "maximum_closure_error="
        f"{max(coarse_result[1], fine_result[1]):.16e}"
    )
    print(f"maximum_activity_change={maximum_activity_change:.16e}")
    print(f"maximum_geometry_average_error={geometry_error:.16e}")
    print(f"maximum_state_average_error={state_error:.16e}")
    print("selected reactive EB AMR 2D regression: PASS")


if __name__ == "__main__":
    main()
