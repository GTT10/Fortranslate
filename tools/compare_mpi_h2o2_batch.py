#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

EXPECTED_REACTORS = 11
PARITY_TOLERANCE = 5.0e-13


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
        "reactor",
        "final_time",
        "relative_energy_error",
        "closure_error",
        "minimum_mass_fraction",
        "maximum_species_change",
        "accepted_steps",
        "newton_iterations",
    }
    missing = required.difference(names)
    if missing:
        raise AssertionError(f"missing columns: {sorted(missing)}")
    if len(rows) != EXPECTED_REACTORS:
        raise AssertionError(
            f"expected {EXPECTED_REACTORS} reactors, found {len(rows)}"
        )
    column = {name: names.index(name) for name in required}
    if not all(math.isfinite(value) for row in rows for value in row):
        raise AssertionError("MPI chemistry output contains non-finite values")
    reactor_ids = [round(row[column["reactor"]]) for row in rows]
    if reactor_ids != list(range(1, EXPECTED_REACTORS + 1)):
        raise AssertionError("MPI chemistry output has incorrect reactor ordering")
    if max(row[column["relative_energy_error"]] for row in rows) > 5.0e-9:
        raise AssertionError("MPI chemistry energy regression failed")
    if max(row[column["closure_error"]] for row in rows) > 5.0e-12:
        raise AssertionError("MPI chemistry closure regression failed")
    if min(row[column["minimum_mass_fraction"]] for row in rows) < -2.0e-13:
        raise AssertionError("MPI chemistry positivity regression failed")
    if max(abs(row[column["final_time"]] - 5.0e-7) for row in rows) > 5.0e-15:
        raise AssertionError("MPI chemistry final-time regression failed")
    if min(row[column["accepted_steps"]] for row in rows) < 1.0:
        raise AssertionError("MPI chemistry did not accept a time step")
    if max(row[column["newton_iterations"]] for row in rows) < 1.0:
        raise AssertionError("MPI chemistry did not execute Newton iterations")
    if max(row[column["maximum_species_change"]] for row in rows) <= 1.0e-16:
        raise AssertionError("MPI chemistry response was trivial")


if len(sys.argv) < 3:
    raise SystemExit(
        "usage: compare_mpi_h2o2_batch.py reference.csv candidate.csv [...]"
    )

names, reference = read_csv(sys.argv[1])
validate_reference(names, reference)
for filename in sys.argv[2:]:
    candidate_names, candidate = read_csv(filename)
    if candidate_names != names or len(candidate) != len(reference):
        raise AssertionError(f"{filename}: MPI chemistry output shape mismatch")
    maximum = 0.0
    for left, right in zip(reference, candidate):
        for a, b in zip(left, right):
            maximum = max(maximum, abs(a - b) / max(1.0, abs(a), abs(b)))
    print(f"{filename}: maximum_relative_difference={maximum:.16e}")
    if maximum > PARITY_TOLERANCE:
        raise AssertionError(f"{filename}: MPI chemistry rank-count parity failed")
