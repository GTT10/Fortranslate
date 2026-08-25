#!/usr/bin/env python3
"""Check public reactive EB multipatch checkpoint/restart parity."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


FINAL_TIME = 3.0e-6


def load_rows(path: Path, expected_rows: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected_rows:
        raise AssertionError(
            f"{path.name}: expected {expected_rows} rows, got {len(rows)}"
        )
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError(f"{path.name}: nonpositive state")
        species = [name for name in row if name.startswith("Y_")]
        if abs(sum(float(row[name]) for name in species) - 1.0) > 8.0e-13:
            raise AssertionError(f"{path.name}: species closure drift")
    return rows


def close(actual: float, expected: float, tolerance: float = 3.0e-10) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def compare_rows(
    reference: list[dict[str, str]], restarted: list[dict[str, str]], label: str
) -> None:
    if reference[0].keys() != restarted[0].keys():
        raise AssertionError(f"{label}: output columns differ")
    for row_index, (expected, actual) in enumerate(zip(reference, restarted)):
        for name in expected:
            if not close(float(actual[name]), float(expected[name])):
                raise AssertionError(
                    f"{label} row {row_index} {name}: "
                    f"{actual[name]} != {expected[name]}"
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference", required=True, nargs=3, type=Path)
    parser.add_argument("--stopped", required=True, nargs=3, type=Path)
    parser.add_argument("--restarted", required=True, nargs=3, type=Path)
    args = parser.parse_args()

    checkpoint_lines = args.checkpoint.read_text(encoding="utf-8").splitlines()
    if checkpoint_lines[0] != "PELEF_REACTIVE_EB_AMR_PATCH_SET_2D_CHECKPOINT":
        raise AssertionError("multipatch checkpoint magic mismatch")
    header = [int(value) for value in checkpoint_lines[1].split()]
    if len(header) != 4 or header[0] != 1 or header[3] != 2:
        raise AssertionError(f"checkpoint did not preserve two patches: {header}")
    if checkpoint_lines[-1] != "END_CHECKPOINT":
        raise AssertionError("multipatch checkpoint end marker mismatch")

    expected_rows = [14 * 14, 10 * 10, 10 * 10]
    references = [
        load_rows(path, count)
        for path, count in zip(args.reference, expected_rows, strict=True)
    ]
    stopped = [
        load_rows(path, count)
        for path, count in zip(args.stopped, expected_rows, strict=True)
    ]
    restarted = [
        load_rows(path, count)
        for path, count in zip(args.restarted, expected_rows, strict=True)
    ]

    stopped_time = float(stopped[0][0]["time"])
    if not 0.0 < stopped_time < FINAL_TIME:
        raise AssertionError("multipatch checkpoint run did not stop early")
    for rows in references + restarted:
        if abs(float(rows[0]["time"]) - FINAL_TIME) > 5.0e-19:
            raise AssertionError("multipatch final time mismatch")
    if {int(row["cell_type"]) for row in references[0]} != {0, 1, 2}:
        raise AssertionError("reference root lacks complete EB coverage")

    for index, (expected, actual) in enumerate(zip(references, restarted)):
        compare_rows(expected, actual, f"level {index}")
    first_x = max(float(row["x"]) for row in restarted[1])
    second_x = min(float(row["x"]) for row in restarted[2])
    if first_x >= second_x:
        raise AssertionError("restarted child ordering or separation changed")
    print("check_reactive_eb_amr_multipatch_restart_2d: PASS")


if __name__ == "__main__":
    main()
