#!/usr/bin/env python3
"""Check public arbitrary-depth EB patch-tree checkpoint/restart parity."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


MAGIC = "PELEF_REACTIVE_AMR_EB_PATCH_TREE_2D"
FINAL_TIME = 3.0e-9
IDENTITY = ("level", "patch", "i", "j")


def load(path: Path) -> dict[tuple[int, int, int, int], dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError(f"{path.name}: empty composite output")
    result: dict[tuple[int, int, int, int], dict[str, str]] = {}
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
        key = tuple(int(row[name]) for name in IDENTITY)
        if key in result:
            raise AssertionError(f"{path.name}: duplicate cell identity {key}")
        result[key] = row
    return result


def compare(
    reference: dict[tuple[int, int, int, int], dict[str, str]],
    restarted: dict[tuple[int, int, int, int], dict[str, str]],
) -> None:
    if reference.keys() != restarted.keys():
        raise AssertionError("reference/restart composite topology mismatch")
    for key in reference:
        expected = reference[key]
        actual = restarted[key]
        if expected.keys() != actual.keys():
            raise AssertionError(f"cell {key}: column mismatch")
        for name in expected:
            lhs = float(actual[name])
            rhs = float(expected[name])
            tolerance = 3.0e-10 * max(1.0, abs(rhs))
            if abs(lhs - rhs) > tolerance:
                raise AssertionError(f"cell {key} {name}: {lhs} != {rhs}")


def check_checkpoint(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != MAGIC:
        raise AssertionError("incorrect patch-tree checkpoint magic")
    if lines[-1] != "END_CHECKPOINT":
        raise AssertionError("incomplete patch-tree checkpoint")
    schema, species, nvar, levels = (int(value) for value in lines[1].split())
    if schema != 1 or species != 7 or nvar <= species or levels != 4:
        raise AssertionError("invalid patch-tree checkpoint header")


def unique_time(rows: dict[tuple[int, int, int, int], dict[str, str]]) -> float:
    times = {float(row["time"]) for row in rows.values()}
    if len(times) != 1:
        raise AssertionError("composite output contains multiple times")
    return next(iter(times))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--stopped", required=True, type=Path)
    parser.add_argument("--restarted", required=True, type=Path)
    args = parser.parse_args()

    check_checkpoint(args.checkpoint)
    reference = load(args.reference)
    stopped = load(args.stopped)
    restarted = load(args.restarted)
    for label, rows in (("reference", reference), ("stopped", stopped),
                        ("restarted", restarted)):
        levels = {key[0] for key in rows}
        if levels != {0, 1, 2, 3}:
            raise AssertionError(f"{label}: expected four levels, found {levels}")
    stopped_time = unique_time(stopped)
    if not 0.0 < stopped_time < FINAL_TIME:
        raise AssertionError("checkpoint run did not stop before final time")
    if abs(unique_time(reference) - FINAL_TIME) > 2.0e-20:
        raise AssertionError("reference did not reach final time")
    if abs(unique_time(restarted) - FINAL_TIME) > 2.0e-20:
        raise AssertionError("restart did not reach final time")
    compare(reference, restarted)
    print("check_reactive_eb_patch_tree_restart_2d: PASS")


if __name__ == "__main__":
    main()
