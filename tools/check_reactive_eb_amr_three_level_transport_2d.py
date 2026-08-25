#!/usr/bin/env python3
"""Check public three-level reactive EB AMR thermal conduction."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load(path: Path, expected_rows: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected_rows:
        raise AssertionError(
            f"{path.name}: expected {expected_rows} rows, got {len(rows)}"
        )
    species = [name for name in rows[0] if name.startswith("Y_")]
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
        if abs(float(row["time"]) - 2.0e-7) > 4.0e-20:
            raise AssertionError(f"{path.name}: incorrect final time")
        if int(row["cell_type"]) != 0:
            if float(row["rho"]) <= 0.0 or float(row["temperature"]) <= 0.0:
                raise AssertionError(f"{path.name}: invalid active thermodynamics")
            closure = sum(float(row[name]) for name in species)
            if abs(closure - 1.0) > 8.0e-12:
                raise AssertionError(f"{path.name}: species closure drift")
    return rows


def active_temperatures(rows: list[dict[str, str]]) -> list[float]:
    return [
        float(row["temperature"])
        for row in rows
        if int(row["cell_type"]) != 0
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    for level in ("root", "middle", "finest"):
        parser.add_argument(f"--reference-{level}", required=True, type=Path)
        parser.add_argument(f"--transport-{level}", required=True, type=Path)
    args = parser.parse_args()

    expected_rows = {"root": 8 * 8, "middle": 12 * 12, "finest": 16 * 16}
    reference_hierarchy_temperature: list[float] = []
    transport_hierarchy_temperature: list[float] = []
    maximum_change = 0.0
    for level, count in expected_rows.items():
        reference = load(getattr(args, f"reference_{level}"), count)
        transported = load(getattr(args, f"transport_{level}"), count)
        if [row["cell_type"] for row in reference] != [
            row["cell_type"] for row in transported
        ]:
            raise AssertionError(f"{level}: EB classification changed")
        reference_temperature = active_temperatures(reference)
        transport_temperature = active_temperatures(transported)
        reference_hierarchy_temperature.extend(reference_temperature)
        transport_hierarchy_temperature.extend(transport_temperature)
        maximum_change = max(
            maximum_change,
            max(
                abs(lhs - rhs)
                for lhs, rhs in zip(reference_temperature, transport_temperature)
            ),
        )

    reference_span = max(reference_hierarchy_temperature) - min(
        reference_hierarchy_temperature
    )
    transport_span = max(transport_hierarchy_temperature) - min(
        transport_hierarchy_temperature
    )
    if reference_span - transport_span <= 1.0e-8:
        raise AssertionError(
            "three-level conduction did not reduce the hierarchy temperature span"
        )
    if maximum_change <= 1.0e-8:
        raise AssertionError("three-level thermal transport produced no change")
    print("check_reactive_eb_amr_three_level_transport_2d: PASS")


if __name__ == "__main__":
    main()
