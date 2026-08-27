#!/usr/bin/env python3
"""Validate the runnable arbitrary-depth reactive 2D EB patch-tree case."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--require-x-upper-boundary", action="store_true")
    args = parser.parse_args()

    with args.output.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError("empty patch-tree output")

    required = {
        "level",
        "patch",
        "i",
        "j",
        "cell_dx",
        "cell_dy",
        "time",
        "volume_fraction",
        "cell_type",
        "rho",
        "pressure",
        "temperature",
    }
    if not required.issubset(rows[0]):
        raise AssertionError("incomplete patch-tree output header")
    species = [name for name in rows[0] if name.startswith("Y_")]
    if len(species) != 7:
        raise AssertionError(f"expected 7 species columns, found {len(species)}")

    identities: set[tuple[int, int, int, int]] = set()
    levels: set[int] = set()
    x_upper_levels: set[int] = set()
    cell_types: set[int] = set()
    root_dx = 0.012 / 12.0
    for row in rows:
        values = [float(value) for value in row.values()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("nonfinite patch-tree output")
        level = int(row["level"])
        identity = (level, int(row["patch"]), int(row["i"]), int(row["j"]))
        if identity in identities:
            raise AssertionError(f"duplicate patch-tree cell identity {identity}")
        identities.add(identity)
        levels.add(level)
        cell_types.add(int(row["cell_type"]))
        expected_spacing = root_dx / 2**level
        if abs(float(row["cell_dx"]) - expected_spacing) > 2.0e-14 * root_dx:
            raise AssertionError("incorrect level x spacing")
        if abs(float(row["cell_dy"]) - expected_spacing) > 2.0e-14 * root_dx:
            raise AssertionError("incorrect level y spacing")
        if abs(float(row["x"]) + 0.5 * float(row["cell_dx"]) - 0.012) <= (
            2.0e-14 * root_dx
        ):
            x_upper_levels.add(level)
        if abs(float(row["time"]) - 1.0e-9) > 2.0e-22:
            raise AssertionError("incorrect patch-tree final time")
        if not 0.0 <= float(row["volume_fraction"]) <= 1.0:
            raise AssertionError("invalid EB volume fraction")
        if float(row["rho"]) <= 0.0 or float(row["pressure"]) <= 0.0:
            raise AssertionError("nonpositive reactive state")
        if float(row["temperature"]) <= 0.0:
            raise AssertionError("nonpositive reactive temperature")
        closure = sum(float(row[name]) for name in species)
        if abs(closure - 1.0) > 8.0e-12:
            raise AssertionError("species mass-fraction closure drift")

    if levels != {0, 1, 2, 3}:
        raise AssertionError(f"expected four populated levels, found {levels}")
    if cell_types != {0, 1, 2}:
        raise AssertionError(f"incomplete EB cell classes {cell_types}")
    if args.require_x_upper_boundary and x_upper_levels != levels:
        raise AssertionError(
            "not every populated AMR level reaches the x-upper physical boundary: "
            f"{x_upper_levels}"
        )
    print("check_reactive_eb_patch_tree_2d: PASS")


if __name__ == "__main__":
    main()
