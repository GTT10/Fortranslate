#!/usr/bin/env python3
"""Check static three-level EB AMR checkpoint/restart parity."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference", required=True, nargs=3, type=Path)
    parser.add_argument("--stopped", required=True, nargs=3, type=Path)
    parser.add_argument("--restarted", required=True, nargs=3, type=Path)
    args = parser.parse_args()

    text = args.checkpoint.read_text(encoding="utf-8")
    if not text.startswith("PELEF_REACTIVE_EB_AMR_THREE_LEVEL_2D_CHECKPOINT"):
        raise AssertionError("incorrect three-level checkpoint magic")
    if not text.rstrip().endswith("END_CHECKPOINT"):
        raise AssertionError("incomplete three-level checkpoint")

    counts = [8 * 8, 12 * 12, 16 * 16]
    reference = [load(path, count) for path, count in zip(args.reference, counts)]
    stopped = [load(path, count) for path, count in zip(args.stopped, counts)]
    restarted = [load(path, count) for path, count in zip(args.restarted, counts)]
    stopped_times = {float(row["time"]) for level in stopped for row in level}
    if len(stopped_times) != 1 or next(iter(stopped_times)) >= 4.0e-7:
        raise AssertionError("checkpoint run did not stop before final time")
    for index, label in enumerate(("root", "middle", "finest")):
        compare(reference[index], restarted[index], label)
    print("check_reactive_eb_amr_three_level_restart_2d: PASS")


if __name__ == "__main__":
    main()
