#!/usr/bin/env python3
"""Check public adiabatic-slip EB thermal conduction."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 16 * 16:
        raise AssertionError(f"{path.name}: expected 256 rows, got {len(rows)}")
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--transport", required=True, type=Path)
    args = parser.parse_args()

    reference = load(args.reference)
    transported = load(args.transport)
    if [row["cell_type"] for row in reference] != [
        row["cell_type"] for row in transported
    ]:
        raise AssertionError("EB classification changed")
    active_reference = [
        float(row["temperature"])
        for row in reference
        if int(row["cell_type"]) != 0
    ]
    active_transport = [
        float(row["temperature"])
        for row in transported
        if int(row["cell_type"]) != 0
    ]
    reference_span = max(active_reference) - min(active_reference)
    transport_span = max(active_transport) - min(active_transport)
    if not transport_span < reference_span - 1.0e-8:
        raise AssertionError(
            f"thermal conduction did not reduce span: {transport_span} >= "
            f"{reference_span}"
        )
    changed = max(
        abs(lhs - rhs) for lhs, rhs in zip(active_reference, active_transport)
    )
    if changed <= 1.0e-8:
        raise AssertionError("thermal transport produced no measurable change")
    cut_pairs = [
        (reference_row, transport_row)
        for reference_row, transport_row in zip(reference, transported)
        if int(reference_row["cell_type"]) == 1
    ]
    wall_heating = max(
        float(transport_row["temperature"])
        - float(reference_row["temperature"])
        for reference_row, transport_row in cut_pairs
    )
    if wall_heating <= 1.0e-8:
        raise AssertionError("configured isothermal EB wall did not heat cut cells")
    wall_speed = max(abs(float(row["v"])) for _, row in cut_pairs)
    if wall_speed <= 1.0e-12:
        raise AssertionError("configured moving no-slip EB wall transferred no momentum")
    print("check_reactive_eb_transport_2d: PASS")


if __name__ == "__main__":
    main()
