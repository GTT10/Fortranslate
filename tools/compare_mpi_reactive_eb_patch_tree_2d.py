#!/usr/bin/env python3
"""Check rank-count parity for the public sparse MPI reactive EB tree."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


IDENTITY = ("level", "patch", "i", "j")
FINAL_TIME = 1.0e-9


def load(path: Path) -> dict[tuple[int, int, int, int], dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError(f"{path.name}: empty composite output")
    species = [name for name in rows[0] if name.startswith("Y_")]
    if len(species) != 7:
        raise AssertionError(f"{path.name}: expected seven species")
    result: dict[tuple[int, int, int, int], dict[str, str]] = {}
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite field")
        key = tuple(int(row[name]) for name in IDENTITY)
        if key in result:
            raise AssertionError(f"{path.name}: duplicate identity {key}")
        if abs(float(row["time"]) - FINAL_TIME) > 2.0e-22:
            raise AssertionError(f"{path.name}: incorrect final time")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError(f"{path.name}: nonpositive state")
        if abs(sum(float(row[name]) for name in species) - 1.0) > 8.0e-12:
            raise AssertionError(f"{path.name}: species closure drift")
        result[key] = row
    levels = {key[0] for key in result}
    if levels != {0, 1, 2, 3}:
        raise AssertionError(f"{path.name}: expected four levels, found {levels}")
    return result


def compare(
    reference: dict[tuple[int, int, int, int], dict[str, str]],
    candidate: dict[tuple[int, int, int, int], dict[str, str]],
    label: str,
) -> None:
    if reference.keys() != candidate.keys():
        raise AssertionError(f"{label}: composite topology mismatch")
    for key in reference:
        expected = reference[key]
        actual = candidate[key]
        if expected.keys() != actual.keys():
            raise AssertionError(f"{label} cell {key}: column mismatch")
        for name in expected:
            lhs = float(actual[name])
            rhs = float(expected[name])
            tolerance = 3.0e-10 * max(1.0, abs(rhs))
            if abs(lhs - rhs) > tolerance:
                raise AssertionError(
                    f"{label} cell {key} {name}: {lhs} != {rhs}"
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("outputs", nargs=4, type=Path)
    args = parser.parse_args()

    reference = load(args.outputs[0])
    for ranks, path in zip((2, 4, 8), args.outputs[1:]):
        compare(reference, load(path), f"{ranks} ranks")
    print("compare_mpi_reactive_eb_patch_tree_2d: PASS")


if __name__ == "__main__":
    main()
