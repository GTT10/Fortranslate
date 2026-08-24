#!/usr/bin/env python3
"""Check active-cell chemistry parity and covered-cell inertia for EB 2D."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 16:
        raise AssertionError(f"{path}: expected 16 cells, found {len(rows)}")
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path}: nonfinite value")
    return rows


def close(actual: float, expected: float, tolerance: float = 5.0e-10) -> bool:
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reactive", required=True, type=Path)
    parser.add_argument("--inert", required=True, type=Path)
    args = parser.parse_args()

    reference = load_rows(args.reference)
    reactive = load_rows(args.reactive)
    inert = load_rows(args.inert)
    species_columns = [name for name in reference[0] if name.startswith("Y_")]
    state_columns = [
        "rho",
        "u",
        "v",
        "w",
        "pressure",
        "temperature",
        "rhoE",
        *species_columns,
    ]
    geometry_columns = [
        "x",
        "y",
        "volume_fraction",
        "cell_type",
        "boundary_length",
        "boundary_normal_x",
        "boundary_normal_y",
    ]

    cell_types: set[int] = set()
    reaction_change = 0.0
    initial_mass = 0.0
    final_mass = 0.0
    initial_energy = 0.0
    final_energy = 0.0
    for index, (reference_row, reactive_row, inert_row) in enumerate(
        zip(reference, reactive, inert, strict=True)
    ):
        expected_x = (index % 4 + 0.5) * 0.001
        expected_y = (index // 4 + 0.5) * 0.001
        if abs(float(reactive_row["x"]) - expected_x) > 2.0e-16 or abs(
            float(reactive_row["y"]) - expected_y
        ) > 2.0e-16:
            raise AssertionError("incorrect reactive EB coordinate order")
        if abs(float(reference_row["x"]) - expected_x) > 2.0e-16 or abs(
            float(reference_row["y"]) - expected_y
        ) > 2.0e-16:
            raise AssertionError("incorrect regular reference coordinate order")
        for name in geometry_columns:
            if reactive_row[name] != inert_row[name]:
                raise AssertionError(f"geometry changed during chemistry: {name}")

        cell_type = int(reactive_row["cell_type"])
        cell_types.add(cell_type)
        if cell_type == 0:
            for name in state_columns:
                if reactive_row[name] != inert_row[name]:
                    raise AssertionError(f"covered cell changed: {name}")
        else:
            for name in state_columns:
                actual = float(reactive_row[name])
                expected = float(reference_row[name])
                if not close(actual, expected):
                    raise AssertionError(
                        f"active-cell chemistry mismatch for {name}: "
                        f"{actual} != {expected}"
                    )
            reaction_change = max(
                reaction_change,
                abs(
                    float(reactive_row["temperature"])
                    - float(inert_row["temperature"])
                ),
                *( 
                    abs(float(reactive_row[name]) - float(inert_row[name]))
                    for name in species_columns
                ),
            )

        volume = float(reactive_row["volume_fraction"])
        initial_mass += volume * float(inert_row["rho"])
        final_mass += volume * float(reactive_row["rho"])
        initial_energy += volume * float(inert_row["rhoE"])
        final_energy += volume * float(reactive_row["rhoE"])

    if cell_types != {0, 1, 2}:
        raise AssertionError(f"missing EB cell class: {sorted(cell_types)}")
    if reaction_change <= 1.0e-14:
        raise AssertionError("reactive EB run did not change chemistry")
    if not close(final_mass, initial_mass, 2.0e-12):
        raise AssertionError("reactive EB mass conservation failure")
    if not close(final_energy, initial_energy, 2.0e-12):
        raise AssertionError("reactive EB energy conservation failure")

    print("check_reactive_eb_chemistry_2d: PASS")


if __name__ == "__main__":
    main()
