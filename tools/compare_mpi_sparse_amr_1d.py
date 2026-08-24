#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

PARITY_TOLERANCE = 5.0e-13


def read_csv(path: str) -> tuple[list[str], list[list[float]]]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames
        if names is None:
            raise AssertionError(f"{path}: missing CSV header")
        rows = [[float(row[name]) for name in names] for row in reader]
    return names, rows


def validate_reference(
    names: list[str], rows: list[list[float]], expected_time: float = 1.0e-10
) -> None:
    required = {
        "level",
        "cell_dx",
        "time",
        "x",
        "rho",
        "pressure",
        "temperature",
        "rhoE",
    }
    missing = required.difference(names)
    if missing:
        raise AssertionError(f"missing columns: {sorted(missing)}")
    species = [name for name in names if name.startswith("Y_")]
    if not species:
        raise AssertionError("missing species mass-fraction columns")
    if len(rows) <= 24:
        raise AssertionError("output did not contain refined composite cells")
    column = {name: names.index(name) for name in names}
    if not all(math.isfinite(value) for row in rows for value in row):
        raise AssertionError("sparse MPI AMR output contains non-finite values")
    levels = {int(row[column["level"]]) for row in rows}
    if max(levels) < 2:
        raise AssertionError("sparse MPI AMR case did not reach three levels")
    if min(row[column["cell_dx"]] for row in rows) <= 0.0:
        raise AssertionError("sparse MPI AMR output has invalid cell spacing")
    if min(row[column["rho"]] for row in rows) <= 0.0:
        raise AssertionError("sparse MPI AMR output has non-positive density")
    if min(row[column["pressure"]] for row in rows) <= 0.0:
        raise AssertionError("sparse MPI AMR output has non-positive pressure")
    if min(row[column["temperature"]] for row in rows) <= 0.0:
        raise AssertionError("sparse MPI AMR output has non-positive temperature")
    time_tolerance = max(5.0e-18, 5.0e-8 * abs(expected_time))
    if any(
        abs(row[column["time"]] - expected_time) > time_tolerance for row in rows
    ):
        raise AssertionError("sparse MPI AMR output has the wrong final time")
    coordinates = [row[column["x"]] for row in rows]
    if any(right <= left for left, right in zip(coordinates, coordinates[1:])):
        raise AssertionError("composite AMR cells are not in physical order")
    minimum_species = min(row[column[name]] for row in rows for name in species)
    if minimum_species < -2.0e-12:
        raise AssertionError("sparse MPI AMR output violated species positivity")
    maximum_closure = max(
        abs(sum(row[column[name]] for name in species) - 1.0) for row in rows
    )
    if maximum_closure > 5.0e-10:
        raise AssertionError("sparse MPI AMR output violated species closure")


def compare_rows(
    names: list[str], reference: list[list[float]], filename: str
) -> float:
    candidate_names, candidate = read_csv(filename)
    if candidate_names != names or len(candidate) != len(reference):
        raise AssertionError(f"{filename}: sparse MPI AMR output shape mismatch")
    maximum = 0.0
    for left, right in zip(reference, candidate):
        for a, b in zip(left, right):
            maximum = max(maximum, abs(a - b) / max(1.0, abs(a), abs(b)))
    return maximum


def main(argv: list[str]) -> None:
    if len(argv) < 3:
        raise SystemExit(
            "usage: compare_mpi_sparse_amr_1d.py reference.csv candidate.csv [...]"
        )
    names, reference = read_csv(argv[1])
    validate_reference(names, reference)
    for filename in argv[2:]:
        maximum = compare_rows(names, reference, filename)
        print(f"{filename}: maximum_relative_difference={maximum:.16e}")
        if maximum > PARITY_TOLERANCE:
            raise AssertionError(
                f"{filename}: sparse MPI AMR rank-count parity failed"
            )


if __name__ == "__main__":
    main(sys.argv)
