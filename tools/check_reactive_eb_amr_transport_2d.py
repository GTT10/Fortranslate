#!/usr/bin/env python3
"""Check public two-level reactive EB AMR configured-wall transport."""

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
    parser.add_argument("--reference-coarse", required=True, type=Path)
    parser.add_argument("--reference-fine", required=True, type=Path)
    parser.add_argument("--transport-coarse", required=True, type=Path)
    parser.add_argument("--transport-fine", required=True, type=Path)
    args = parser.parse_args()

    reference_coarse = load(args.reference_coarse, 12 * 12)
    reference_fine = load(args.reference_fine, 20 * 20)
    transport_coarse = load(args.transport_coarse, 12 * 12)
    transport_fine = load(args.transport_fine, 20 * 20)
    for reference, transported in (
        (reference_coarse, transport_coarse),
        (reference_fine, transport_fine),
    ):
        if [row["cell_type"] for row in reference] != [
            row["cell_type"] for row in transported
        ]:
            raise AssertionError("EB classification changed")

    reference_temperature = active_temperatures(reference_fine)
    transport_temperature = active_temperatures(transport_fine)
    reference_span = max(reference_temperature) - min(reference_temperature)
    transport_span = max(transport_temperature) - min(transport_temperature)
    if not transport_span < reference_span - 1.0e-8:
        raise AssertionError(
            "AMR conduction did not reduce fine-level span: "
            f"{transport_span} >= {reference_span}"
        )
    changed = max(
        abs(lhs - rhs)
        for lhs, rhs in zip(reference_temperature, transport_temperature)
    )
    if changed <= 1.0e-8:
        raise AssertionError("AMR thermal transport produced no measurable change")
    cut_pairs = [
        (reference_row, transport_row)
        for reference_rows, transport_rows in (
            (reference_coarse, transport_coarse),
            (reference_fine, transport_fine),
        )
        for reference_row, transport_row in zip(reference_rows, transport_rows)
        if int(reference_row["cell_type"]) == 1
    ]
    if not cut_pairs:
        raise AssertionError("AMR wall regression contains no cut cells")
    wall_heating = max(
        float(transport_row["temperature"])
        - float(reference_row["temperature"])
        for reference_row, transport_row in cut_pairs
    )
    if wall_heating <= 1.0e-8:
        raise AssertionError("configured AMR isothermal wall did not heat cut cells")
    wall_velocity_response = max(
        abs(float(transport_row["v"]) - float(reference_row["v"]))
        for reference_row, transport_row in cut_pairs
    )
    if wall_velocity_response <= 1.0e-12:
        raise AssertionError(
            "configured moving no-slip AMR wall transferred no momentum"
        )
    print("check_reactive_eb_amr_transport_2d: PASS")


if __name__ == "__main__":
    main()
