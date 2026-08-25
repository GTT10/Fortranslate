#!/usr/bin/env python3
"""Check dynamic three-level topology and cadence checkpoint/restart parity."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


MAGIC = "PELEF_REACTIVE_EB_AMR_DYNAMIC_THREE_LEVEL_2D_CHECKPOINT"


def load(path: Path, expected: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected:
        raise AssertionError(f"{path.name}: expected {expected} rows, got {len(rows)}")
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
    return rows


def compare(
    reference: list[dict[str, str]], restarted: list[dict[str, str]], label: str
) -> None:
    if len(reference) != len(restarted):
        raise AssertionError(f"{label}: row-count mismatch")
    for index, (expected, actual) in enumerate(zip(reference, restarted)):
        if expected.keys() != actual.keys():
            raise AssertionError(f"{label}: column mismatch")
        for name in expected:
            lhs = float(actual[name])
            rhs = float(expected[name])
            tolerance = 3.0e-10 * max(1.0, abs(rhs))
            if abs(lhs - rhs) > tolerance:
                raise AssertionError(
                    f"{label} row {index} {name}: {lhs} != {rhs}"
                )


def check_checkpoint(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != MAGIC:
        raise AssertionError("incorrect dynamic three-level checkpoint magic")
    if lines[-1] != "END_CHECKPOINT":
        raise AssertionError("incomplete dynamic three-level checkpoint")
    schema, species_count, _ = (int(value) for value in lines[1].split())
    if schema != 1 or species_count < 1:
        raise AssertionError("invalid dynamic checkpoint schema header")
    metadata = 2 + species_count
    controls = [int(value) for value in lines[metadata + 11].split()]
    if controls != [1, 1, 0, 4, 4, 0]:
        raise AssertionError(f"incorrect regrid controls {controls}")
    finest_patch = [int(value) for value in lines[metadata + 14].split()]
    if finest_patch == [6, 9, 6, 9, 2]:
        raise AssertionError("checkpoint retained the configured 8x8 seed")
    if (
        (finest_patch[1] - finest_patch[0] + 1) * finest_patch[4] != 22
        or (finest_patch[3] - finest_patch[2] + 1) * finest_patch[4] != 28
    ):
        raise AssertionError(f"incorrect stored finest topology {finest_patch}")
    time_record = lines[metadata + 15].split()
    if len(time_record) != 5 or int(time_record[3]) < 1:
        raise AssertionError("checkpoint did not preserve the regrid count")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference", required=True, nargs=3, type=Path)
    parser.add_argument("--stopped", required=True, nargs=3, type=Path)
    parser.add_argument("--restarted", required=True, nargs=3, type=Path)
    args = parser.parse_args()

    check_checkpoint(args.checkpoint)
    counts = [12 * 12, 20 * 20, 22 * 28]
    reference = [load(path, count) for path, count in zip(args.reference, counts)]
    stopped = [load(path, count) for path, count in zip(args.stopped, counts)]
    restarted = [load(path, count) for path, count in zip(args.restarted, counts)]
    stopped_times = {float(row["time"]) for level in stopped for row in level}
    if len(stopped_times) != 1 or next(iter(stopped_times)) >= 2.0e-7:
        raise AssertionError("checkpoint run did not stop before final time")
    for index, label in enumerate(("root", "middle", "finest")):
        compare(reference[index], restarted[index], label)
    print("check_reactive_eb_amr_three_level_dynamic_restart_2d: PASS")


if __name__ == "__main__":
    main()
