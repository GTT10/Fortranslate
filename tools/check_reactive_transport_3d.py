#!/usr/bin/env python3
"""Validate the paired public full-H2O2 3D transport hotspot."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2")
MOLE_FRACTIONS = {
    "H2": 0.29570,
    "H": 1.0e-5,
    "O": 1.0e-5,
    "O2": 0.14784,
    "OH": 1.0e-5,
    "H2O": 0.0,
    "HO2": 0.0,
    "H2O2": 0.0,
    "AR": 0.0,
    "N2": 0.55643,
}
MOLECULAR_WEIGHTS = {
    "H2": 2.016,
    "H": 1.008,
    "O": 15.999,
    "O2": 31.998,
    "OH": 17.007,
    "H2O": 18.015,
    "HO2": 33.006,
    "H2O2": 34.014,
    "AR": 39.950,
    "N2": 28.014,
}
BASE_COLUMNS = (
    "time",
    "x",
    "y",
    "z",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
    "rhou",
    "rhov",
    "rhow",
)
COLUMNS = (
    *BASE_COLUMNS,
    *(f"Y_{name}" for name in SPECIES),
    *(f"rhoY_{name}" for name in SPECIES),
)


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(1.0, abs(expected))


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise AssertionError(f"unexpected columns in {path}: {reader.fieldnames}")
        return [{name: float(row[name]) for name in COLUMNS} for row in reader]


def validate_rows(
    rows: list[dict[str, float]], nx: int, ny: int, nz: int, time: float
) -> tuple[float, float]:
    extents = (nx, ny, nz)
    if len(rows) != math.prod(extents):
        raise AssertionError("3D transport row count does not match the grid")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("3D transport CSV contains a non-finite value")

    relationship_error = 0.0
    closure_error = 0.0
    for index, row in enumerate(rows):
        indices = (
            index % nx,
            (index // nx) % ny,
            index // (nx * ny),
        )
        expected_coordinates = tuple(
            (cell + 0.5) * 0.01 / extent
            for cell, extent in zip(indices, extents)
        )
        coordinate_error = max(
            abs(row[name] - expected)
            for name, expected in zip(("x", "y", "z"), expected_coordinates)
        )
        if coordinate_error > 5.0e-15:
            raise AssertionError(f"row {index} is not in x-fastest cell order")
        if abs(row["time"] - time) > 5.0e-18:
            raise AssertionError(f"row {index} has an unexpected output time")
        if min(row["rho"], row["pressure"], row["temperature"], row["rhoE"]) <= 0.0:
            raise AssertionError("3D transport CSV contains a non-physical state")

        inverse_molecular_weight = sum(
            row[f"Y_{name}"] / MOLECULAR_WEIGHTS[name] for name in SPECIES
        )
        gas_constant = 8.31446261815324e3 * inverse_molecular_weight
        relationship_error = max(
            relationship_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
            relative_error(
                row["pressure"], row["rho"] * gas_constant * row["temperature"]
            ),
        )
        mass_fraction_sum = 0.0
        species_density_sum = 0.0
        for name in SPECIES:
            mass_fraction = row[f"Y_{name}"]
            species_density = row[f"rhoY_{name}"]
            if mass_fraction < 0.0 or species_density < 0.0:
                raise AssertionError("3D transport CSV contains a negative species")
            mass_fraction_sum += mass_fraction
            species_density_sum += species_density
            relationship_error = max(
                relationship_error,
                relative_error(species_density, row["rho"] * mass_fraction),
            )
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )
    if relationship_error > 2.0e-11 or closure_error > 2.0e-11:
        raise AssertionError(
            "3D transport state contract error: "
            f"relation={relationship_error}, closure={closure_error}"
        )
    return relationship_error, closure_error


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--transport", type=Path, required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--time", type=float, required=True)
    args = parser.parse_args()

    control = read_rows(args.control)
    transported = read_rows(args.transport)
    validate_rows(control, args.nx, args.ny, args.nz, args.time)
    relationship_error, closure_error = validate_rows(
        transported, args.nx, args.ny, args.nz, args.time
    )

    cell_volume = (0.01 / args.nx) * (0.01 / args.ny) * (0.01 / args.nz)
    conserved_columns = (
        "rho",
        "rhou",
        "rhov",
        "rhow",
        "rhoE",
        *(f"rhoY_{name}" for name in SPECIES),
    )
    integral_error = 0.0
    for name in conserved_columns:
        control_integral = sum(row[name] for row in control) * cell_volume
        transport_integral = sum(row[name] for row in transported) * cell_volume
        integral_error = max(
            integral_error, relative_error(transport_integral, control_integral)
        )
    if integral_error > 2.0e-11:
        raise AssertionError(f"periodic transport integral error {integral_error}")

    control_span = max(row["temperature"] for row in control) - min(
        row["temperature"] for row in control
    )
    transport_span = max(row["temperature"] for row in transported) - min(
        row["temperature"] for row in transported
    )
    span_reduction = control_span - transport_span
    temperature_change = max(
        abs(after["temperature"] - before["temperature"])
        for before, after in zip(control, transported)
    )
    composition_change = max(
        abs(after[f"Y_{name}"] - before[f"Y_{name}"])
        for before, after in zip(control, transported)
        for name in SPECIES
    )
    maximum_speed = max(
        math.sqrt(row["u"] ** 2 + row["v"] ** 2 + row["w"] ** 2)
        for row in transported
    )
    if span_reduction < 2.0e-2 or temperature_change < 5.0e-2:
        raise AssertionError(
            "Fourier conduction did not measurably smooth the 3D hotspot: "
            f"span={span_reduction}, change={temperature_change}"
        )
    if not 1.0e-10 < composition_change < 1.0e-6:
        raise AssertionError(
            f"unexpected species/barodiffusion response {composition_change}"
        )
    if maximum_speed < 1.0e-4:
        raise AssertionError("transport coupling produced no resolved flow response")

    print(
        "Reactive 3D transport: PASS "
        f"(rows={len(transported)}, dT={temperature_change:.3e}, "
        f"span={span_reduction:.3e}, dY={composition_change:.3e}, "
        f"integrals={integral_error:.3e}, relation={relationship_error:.3e}, "
        f"closure={closure_error:.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
