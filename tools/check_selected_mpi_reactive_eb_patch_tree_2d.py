#!/usr/bin/env python3
"""Validate exact fixed/selected rank parity for sparse MPI reactive EB 2D."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


BASE_COLUMNS = [
    "level",
    "patch",
    "i",
    "j",
    "cell_dx",
    "cell_dy",
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


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_output(
    path: Path, species: list[str], expected_levels: list[int], final_time: float
) -> bytes:
    data = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected_header = [*BASE_COLUMNS, *(f"Y_{name}" for name in species)]
        if reader.fieldnames != expected_header:
            raise AssertionError(f"{path}: sparse MPI EB CSV schema mismatch")
        rows = list(reader)
    if not rows:
        raise AssertionError(f"{path}: sparse MPI EB CSV is empty")

    observed_levels: set[int] = set()
    observed_cell_types: set[int] = set()
    patches_by_level: dict[int, set[int]] = {}
    maximum_closure_error = 0.0
    time_tolerance = max(5.0e-18, 5.0e-8 * abs(final_time))
    for row in rows:
        try:
            values = [float(row[name]) for name in expected_header]
        except ValueError as error:
            raise AssertionError(f"{path}: nonnumeric output: {error}") from error
        if not all(map(math.isfinite, values)):
            raise AssertionError(f"{path}: output contains nonfinite values")
        level = int(row["level"])
        patch = int(row["patch"])
        cell_type = int(row["cell_type"])
        if float(row["level"]) != level or float(row["patch"]) != patch:
            raise AssertionError(f"{path}: nonintegral topology index")
        if float(row["cell_type"]) != cell_type or cell_type not in {0, 1, 2}:
            raise AssertionError(f"{path}: invalid EB cell type")
        observed_levels.add(level)
        observed_cell_types.add(cell_type)
        patches_by_level.setdefault(level, set()).add(patch)
        if abs(float(row["time"]) - final_time) > time_tolerance:
            raise AssertionError(f"{path}: final time mismatch")
        if (
            float(row["rho"]) <= 0.0
            or float(row["pressure"]) <= 0.0
            or float(row["temperature"]) <= 0.0
        ):
            raise AssertionError(f"{path}: nonpositive thermodynamic state")
        volume_fraction = float(row["volume_fraction"])
        if not 0.0 <= volume_fraction <= 1.0:
            raise AssertionError(f"{path}: EB volume fraction outside [0,1]")
        if cell_type == 0 and volume_fraction != 0.0:
            raise AssertionError(f"{path}: covered cell has nonzero volume")
        if cell_type == 1 and not 0.0 < volume_fraction < 1.0:
            raise AssertionError(f"{path}: cut cell lacks partial volume")
        if cell_type == 2 and volume_fraction != 1.0:
            raise AssertionError(f"{path}: regular cell lacks unit volume")
        mass_fractions = [float(row[f"Y_{name}"]) for name in species]
        if min(mass_fractions) < -2.0e-12:
            raise AssertionError(f"{path}: negative species mass fraction")
        maximum_closure_error = max(
            maximum_closure_error, abs(sum(mass_fractions) - 1.0)
        )

    if sorted(observed_levels) != expected_levels:
        raise AssertionError(
            f"{path}: level set changed: {sorted(observed_levels)}"
        )
    if observed_cell_types != {0, 1, 2}:
        raise AssertionError(f"{path}: plane EB cell classes are incomplete")
    if any(patches_by_level[level] != {1} for level in expected_levels):
        raise AssertionError(f"{path}: patch-tree topology changed")
    if maximum_closure_error > 5.0e-10:
        raise AssertionError(f"{path}: species closure failure")
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, action="append", default=[])
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--expected-levels", nargs="+", type=int, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    args = parser.parse_args()

    reference = validate_output(
        args.reference, args.species, args.expected_levels, args.final_time
    )
    for candidate_path in args.candidate:
        candidate = validate_output(
            candidate_path, args.species, args.expected_levels, args.final_time
        )
        if candidate != reference:
            raise AssertionError(
                f"{candidate_path}: exact fixed/selected rank parity failed"
            )
    print(f"reference_sha256={digest(reference)}")
    if args.candidate:
        print("selected_mpi_reactive_eb_patch_tree_2d_exact_match=PASS")
    else:
        print("selected_mpi_reactive_eb_patch_tree_2d_validation=PASS")


if __name__ == "__main__":
    main()
