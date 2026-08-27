#!/usr/bin/env python3
"""Check reacting EB AMR split-run parity across a root-only checkpoint."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 8 * 8:
        raise AssertionError(f"{path.name}: expected 64 rows, got {len(rows)}")
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
    return rows


def close(actual: float, expected: float, tolerance: float = 3.0e-10) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--stopped", required=True, type=Path)
    parser.add_argument("--restarted", required=True, type=Path)
    parser.add_argument("--fine", required=True, nargs="+", type=Path)
    args = parser.parse_args()

    checkpoint_lines = args.checkpoint.read_text(encoding="utf-8").splitlines()
    if checkpoint_lines[0] != "PELEF_REACTIVE_EB_AMR_2D_CHECKPOINT":
        raise AssertionError("checkpoint magic mismatch")
    header = [int(value) for value in checkpoint_lines[1].split()]
    if len(header) != 4 or header[0] != 3 or header[3] != 0:
        raise AssertionError(f"checkpoint did not preserve root-only state: {header}")
    metadata = 2 + header[1]
    if checkpoint_lines[metadata + 12] != "linear":
        raise AssertionError("checkpoint did not preserve linear prolongation")
    if checkpoint_lines[-1] != "END_CHECKPOINT":
        raise AssertionError("checkpoint end marker mismatch")

    reference = load_rows(args.reference)
    stopped = load_rows(args.stopped)
    restarted = load_rows(args.restarted)
    if abs(float(reference[0]["time"]) - 2.0e-7) > 4.0e-20:
        raise AssertionError("reference final time mismatch")
    if abs(float(restarted[0]["time"]) - 2.0e-7) > 4.0e-20:
        raise AssertionError("restart final time mismatch")
    stopped_time = float(stopped[0]["time"])
    if not 0.0 < stopped_time < 2.0e-7:
        raise AssertionError("checkpoint run did not stop before final time")
    if {int(row["cell_type"]) for row in reference} != {0, 1, 2}:
        raise AssertionError("reference output lacks complete EB coverage")
    if any(path.exists() for path in args.fine):
        raise AssertionError("root-only checkpoint lifecycle wrote stale fine output")

    if reference[0].keys() != restarted[0].keys():
        raise AssertionError("restart output columns differ from reference")
    for row_index, (expected, actual) in enumerate(zip(reference, restarted)):
        for name in expected:
            if not close(float(actual[name]), float(expected[name])):
                raise AssertionError(
                    f"row {row_index} {name}: {actual[name]} != {expected[name]}"
                )

    active = [row for row in restarted if int(row["cell_type"]) != 0]
    species = [name for name in restarted[0] if name.startswith("Y_")]
    if max(abs(float(row["temperature"]) - 1200.0) for row in active) <= 1.0e-12:
        raise AssertionError("restart case chemistry did not advance")
    for row in active:
        if abs(sum(float(row[name]) for name in species) - 1.0) > 5.0e-13:
            raise AssertionError("restart species closure drift")
    print("check_reactive_eb_amr_restart_2d: PASS")


if __name__ == "__main__":
    main()
