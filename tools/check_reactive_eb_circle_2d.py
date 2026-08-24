#!/usr/bin/env python3
"""Validate the input-driven stationary reactive EB circle regression."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 400:
        raise AssertionError(f"expected 400 cells, found {len(rows)}")

    required = {
        "time",
        "x",
        "y",
        "volume_fraction",
        "cell_type",
        "boundary_length",
        "boundary_normal_x",
        "boundary_normal_y",
        "rho",
        "u",
        "v",
        "w",
        "pressure",
        "temperature",
        "rhoE",
        "Y_H2",
        "Y_N2",
    }
    missing = required.difference(rows[0])
    if missing:
        raise AssertionError(f"missing columns: {sorted(missing)}")

    cell_types: set[int] = set()
    maximum_speed = 0.0
    maximum_pressure_error = 0.0
    maximum_temperature_error = 0.0
    maximum_closure_error = 0.0
    cut_boundary_length = 0.0
    species_columns = [name for name in rows[0] if name.startswith("Y_")]
    for index, row in enumerate(rows):
        values = [float(value) for value in row.values()]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("nonfinite EB output")
        if abs(float(row["time"]) - 1.0e-7) > 2.0e-20:
            raise AssertionError("incorrect final time")
        expected_x = (index % 20 + 0.5) * 0.01 / 20.0
        expected_y = (index // 20 + 0.5) * 0.01 / 20.0
        if abs(float(row["x"]) - expected_x) > 2.0e-16 or abs(
            float(row["y"]) - expected_y
        ) > 2.0e-16:
            raise AssertionError("incorrect EB CSV coordinate order")
        volume_fraction = float(row["volume_fraction"])
        if not 0.0 <= volume_fraction <= 1.0:
            raise AssertionError("invalid volume fraction")
        cell_type = int(row["cell_type"])
        cell_types.add(cell_type)
        if cell_type == 1:
            cut_boundary_length += float(row["boundary_length"])
            normal_norm = math.hypot(
                float(row["boundary_normal_x"]),
                float(row["boundary_normal_y"]),
            )
            if abs(normal_norm - 1.0) > 2.0e-12:
                raise AssertionError("cut-cell normal is not unit length")
        elif (
            float(row["boundary_length"]) != 0.0
            or float(row["boundary_normal_x"]) != 0.0
            or float(row["boundary_normal_y"]) != 0.0
        ):
            raise AssertionError("non-cut cell has embedded-boundary metrics")
        rho = float(row["rho"])
        pressure = float(row["pressure"])
        temperature = float(row["temperature"])
        if rho <= 0.0 or pressure <= 0.0 or temperature <= 0.0:
            raise AssertionError("nonpositive thermodynamic state")
        maximum_speed = max(
            maximum_speed,
            math.hypot(float(row["u"]), float(row["v"])),
            abs(float(row["w"])),
        )
        maximum_pressure_error = max(
            maximum_pressure_error, abs(pressure - 101325.0)
        )
        maximum_temperature_error = max(
            maximum_temperature_error, abs(temperature - 1000.0)
        )
        closure = sum(float(row[name]) for name in species_columns)
        maximum_closure_error = max(maximum_closure_error, abs(closure - 1.0))

    if cell_types != {0, 1, 2}:
        raise AssertionError(f"missing EB cell class: {sorted(cell_types)}")
    if cut_boundary_length <= 0.0:
        raise AssertionError("cut cells have no embedded-boundary metric")
    if maximum_speed > 2.0e-8:
        raise AssertionError(f"stationary speed drift: {maximum_speed}")
    if maximum_pressure_error > 2.0e-6:
        raise AssertionError(f"pressure drift: {maximum_pressure_error}")
    if maximum_temperature_error > 2.0e-8:
        raise AssertionError(f"temperature drift: {maximum_temperature_error}")
    if maximum_closure_error > 2.0e-13:
        raise AssertionError(f"species closure drift: {maximum_closure_error}")

    print("check_reactive_eb_circle_2d: PASS")


if __name__ == "__main__":
    main()
