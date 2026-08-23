#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

EXPECTED_CELLS = 19
PARITY_TOLERANCE = 5.0e-13
SPECIES_COLUMNS = [
    "rhoY_H2",
    "rhoY_H",
    "rhoY_O",
    "rhoY_O2",
    "rhoY_OH",
    "rhoY_H2O",
    "rhoY_HO2",
    "rhoY_H2O2",
    "rhoY_AR",
    "rhoY_N2",
]


def read_csv(path: str) -> tuple[list[str], list[list[float]]]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames
        if names is None:
            raise AssertionError(f"{path}: missing CSV header")
        rows = [[float(row[name]) for name in names] for row in reader]
    return names, rows


def validate_reference(names: list[str], rows: list[list[float]]) -> None:
    required = {
        "x",
        "rho",
        "rhou",
        "rhov",
        "rhow",
        "rhoE",
        "temperature",
        *SPECIES_COLUMNS,
    }
    missing = required.difference(names)
    if missing:
        raise AssertionError(f"missing columns: {sorted(missing)}")
    if len(rows) != EXPECTED_CELLS:
        raise AssertionError(f"expected {EXPECTED_CELLS} cells, found {len(rows)}")
    column = {name: names.index(name) for name in required}
    if not all(math.isfinite(value) for row in rows for value in row):
        raise AssertionError("MPI reactive output contains non-finite values")
    if min(row[column["rho"]] for row in rows) <= 0.0:
        raise AssertionError("MPI reactive output has non-positive density")
    if min(row[column["temperature"]] for row in rows) <= 0.0:
        raise AssertionError("MPI reactive output has non-positive temperature")
    minimum_species = min(
        row[column[name]] for row in rows for name in SPECIES_COLUMNS
    )
    if minimum_species < -2.0e-12:
        raise AssertionError("MPI reactive output violated species positivity")
    maximum_closure = max(
        abs(sum(row[column[name]] for name in SPECIES_COLUMNS) - row[column["rho"]])
        / row[column["rho"]]
        for row in rows
    )
    if maximum_closure > 5.0e-10:
        raise AssertionError("MPI reactive output violated species closure")


if len(sys.argv) < 3:
    raise SystemExit(
        "usage: compare_mpi_reactive_1d.py reference.csv candidate.csv [...]"
    )

names, reference = read_csv(sys.argv[1])
validate_reference(names, reference)
for filename in sys.argv[2:]:
    candidate_names, candidate = read_csv(filename)
    if candidate_names != names or len(candidate) != len(reference):
        raise AssertionError(f"{filename}: MPI reactive output shape mismatch")
    maximum = 0.0
    for left, right in zip(reference, candidate):
        for a, b in zip(left, right):
            maximum = max(maximum, abs(a - b) / max(1.0, abs(a), abs(b)))
    print(f"{filename}: maximum_relative_difference={maximum:.16e}")
    if maximum > PARITY_TOLERANCE:
        raise AssertionError(f"{filename}: MPI reactive rank-count parity failed")
