#!/usr/bin/env python3
"""Check fixed-H2/O2 chemistry activity and composite AMR invariants."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


MOLECULAR_WEIGHTS = {
    "H2": 2.016,
    "H": 1.008,
    "O": 15.999,
    "O2": 31.998,
    "OH": 17.007,
    "H2O": 18.015,
    "HO2": 33.006,
    "H2O2": 34.014,
    "AR": 39.95,
    "N2": 28.014,
}
ATOMS = {
    "H2": (2.0, 0.0, 0.0),
    "H": (1.0, 0.0, 0.0),
    "O": (0.0, 1.0, 0.0),
    "O2": (0.0, 2.0, 0.0),
    "OH": (1.0, 1.0, 0.0),
    "H2O": (2.0, 1.0, 0.0),
    "HO2": (1.0, 2.0, 0.0),
    "H2O2": (2.0, 2.0, 0.0),
    "AR": (0.0, 0.0, 0.0),
    "N2": (0.0, 0.0, 2.0),
}
EULER_FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoE")
FIXED_SPECIES = tuple(MOLECULAR_WEIGHTS)
BASE_FIELDS = (
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
EXPECTED_FIELDS = (
    *BASE_FIELDS,
    *(f"Y_{name}" for name in FIXED_SPECIES),
    *(f"rhoY_{name}" for name in FIXED_SPECIES),
)


def level_path(prefix: Path, level: str) -> Path:
    return Path(f"{prefix}_{level}.csv")


def read_rows(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise AssertionError(f"missing CSV header: {path}")
        fields = reader.fieldnames
        if tuple(fields) != EXPECTED_FIELDS:
            raise AssertionError(f"unexpected fixed-H2/O2 CSV schema: {path}")
        rows = [
            {name: float(value) for name, value in row.items()}
            for row in reader
        ]
    if not rows or any(
        not math.isfinite(value) for row in rows for value in row.values()
    ):
        raise AssertionError(f"invalid values in {path}")
    return fields, rows


def relative_error(left: float, right: float, floor: float = 1.0) -> float:
    return abs(left - right) / max(floor, abs(left), abs(right))


def validate_level(
    rows: list[dict[str, float]], nx: int, ny: int, nz: int, label: str
) -> None:
    if len(rows) != nx * ny * nz:
        raise AssertionError(
            f"{label}: expected {nx * ny * nz} rows, found {len(rows)}"
        )
    x_values = sorted({row["x"] for row in rows})
    y_values = sorted({row["y"] for row in rows})
    z_values = sorted({row["z"] for row in rows})
    if (len(x_values), len(y_values), len(z_values)) != (nx, ny, nz):
        raise AssertionError(f"{label}: unexpected coordinate topology")
    reference_time = rows[0]["time"]
    if reference_time <= 0.0:
        raise AssertionError(f"{label}: nonpositive output time")
    relation_error = 0.0
    closure_error = 0.0
    for index, row in enumerate(rows):
        i = index % nx
        j = (index // nx) % ny
        k = index // (nx * ny)
        if (row["x"], row["y"], row["z"]) != (
            x_values[i],
            y_values[j],
            z_values[k],
        ):
            raise AssertionError(f"{label}: CSV row order is not k-j-i structured")
        if row["time"] != reference_time:
            raise AssertionError(f"{label}: inconsistent row times")
        if min(row["rho"], row["pressure"], row["temperature"], row["rhoE"]) <= 0:
            raise AssertionError(f"{label}: nonphysical thermodynamic state")
        relation_error = max(
            relation_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
        )
        mass_fraction_sum = sum(row[f"Y_{name}"] for name in FIXED_SPECIES)
        species_density_sum = sum(
            row[f"rhoY_{name}"] for name in FIXED_SPECIES
        )
        if min(row[f"Y_{name}"] for name in FIXED_SPECIES) < 0.0:
            raise AssertionError(f"{label}: negative species mass fraction")
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )
        for name in FIXED_SPECIES:
            relation_error = max(
                relation_error,
                relative_error(
                    row[f"rhoY_{name}"], row["rho"] * row[f"Y_{name}"]
                ),
            )
    if relation_error > 2.0e-11 or closure_error > 2.0e-11:
        raise AssertionError(
            f"{label}: state relation={relation_error}, closure={closure_error}"
        )


def require_matching_layout(
    control: list[dict[str, float]],
    reactive: list[dict[str, float]],
    label: str,
) -> None:
    if len(control) != len(reactive):
        raise AssertionError(f"{label}: control/reactive row counts differ")
    for control_row, reactive_row in zip(control, reactive):
        if tuple(control_row[name] for name in ("time", "x", "y", "z")) != tuple(
            reactive_row[name] for name in ("time", "x", "y", "z")
        ):
            raise AssertionError(f"{label}: control/reactive layouts differ")


def composite_sums(
    coarse: list[dict[str, float]],
    fine: list[dict[str, float]],
    fields: list[str],
    nx: int,
    ny: int,
    nz: int,
    lower: tuple[int, int, int],
    upper: tuple[int, int, int],
    ratio: int,
) -> dict[str, float]:
    if len(coarse) != nx * ny * nz:
        raise AssertionError("unexpected coarse row count")
    covered: set[int] = set()
    for k in range(lower[2] - 1, upper[2]):
        for j in range(lower[1] - 1, upper[1]):
            for i in range(lower[0] - 1, upper[0]):
                covered.add((k * ny + j) * nx + i)
    sums = {
        field: sum(
            row[field] for index, row in enumerate(coarse) if index not in covered
        )
        + sum(row[field] for row in fine) / ratio**3
        for field in fields
    }
    return sums


def element_totals(
    component_sums: dict[str, float], species: list[str]
) -> tuple[float, float, float]:
    totals = [0.0, 0.0, 0.0]
    for name in species:
        amount = component_sums[f"rhoY_{name}"] / MOLECULAR_WEIGHTS[name]
        for element in range(3):
            totals[element] += amount * ATOMS[name][element]
    return tuple(totals)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--control-prefix", type=Path, required=True)
    parser.add_argument("--reactive-prefix", type=Path, required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--i-lower", type=int, required=True)
    parser.add_argument("--i-upper", type=int, required=True)
    parser.add_argument("--j-lower", type=int, required=True)
    parser.add_argument("--j-upper", type=int, required=True)
    parser.add_argument("--k-lower", type=int, required=True)
    parser.add_argument("--k-upper", type=int, required=True)
    parser.add_argument("--ratio", type=int, required=True)
    args = parser.parse_args()

    control_coarse_fields, control_coarse = read_rows(
        level_path(args.control_prefix, "coarse")
    )
    control_fine_fields, control_fine = read_rows(
        level_path(args.control_prefix, "fine")
    )
    reactive_coarse_fields, reactive_coarse = read_rows(
        level_path(args.reactive_prefix, "coarse")
    )
    reactive_fine_fields, reactive_fine = read_rows(
        level_path(args.reactive_prefix, "fine")
    )
    if not (
        control_coarse_fields
        == control_fine_fields
        == reactive_coarse_fields
        == reactive_fine_fields
    ):
        raise AssertionError("control/reactive CSV headers differ")
    fine_nx = (args.i_upper - args.i_lower + 1) * args.ratio
    fine_ny = (args.j_upper - args.j_lower + 1) * args.ratio
    fine_nz = (args.k_upper - args.k_lower + 1) * args.ratio
    validate_level(control_coarse, args.nx, args.ny, args.nz, "control coarse")
    validate_level(control_fine, fine_nx, fine_ny, fine_nz, "control fine")
    validate_level(reactive_coarse, args.nx, args.ny, args.nz, "reactive coarse")
    validate_level(reactive_fine, fine_nx, fine_ny, fine_nz, "reactive fine")
    require_matching_layout(control_coarse, reactive_coarse, "coarse")
    require_matching_layout(control_fine, reactive_fine, "fine")

    species = list(FIXED_SPECIES)
    fields = [*EULER_FIELDS, *(f"rhoY_{name}" for name in species)]
    bounds_lower = (args.i_lower, args.j_lower, args.k_lower)
    bounds_upper = (args.i_upper, args.j_upper, args.k_upper)
    control = composite_sums(
        control_coarse,
        control_fine,
        fields,
        args.nx,
        args.ny,
        args.nz,
        bounds_lower,
        bounds_upper,
        args.ratio,
    )
    reactive = composite_sums(
        reactive_coarse,
        reactive_fine,
        fields,
        args.nx,
        args.ny,
        args.nz,
        bounds_lower,
        bounds_upper,
        args.ratio,
    )
    euler_error = max(
        relative_error(control[name], reactive[name]) for name in EULER_FIELDS
    )
    control_elements = element_totals(control, species)
    reactive_elements = element_totals(reactive, species)
    element_error = max(
        relative_error(left, right, 1.0e-30)
        for left, right in zip(control_elements, reactive_elements)
        if max(abs(left), abs(right)) > 1.0e-30
    )
    species_change = max(
        relative_error(control[f"rhoY_{name}"], reactive[f"rhoY_{name}"])
        for name in species
    )
    temperature_change = max(
        abs(left["temperature"] - right["temperature"])
        for left, right in zip(control_fine, reactive_fine)
    )
    if euler_error > 5.0e-11:
        raise AssertionError(f"reactive Euler conservation error {euler_error}")
    if element_error > 5.0e-10:
        raise AssertionError(f"reactive elemental conservation error {element_error}")
    if species_change <= 1.0e-13 or temperature_change <= 1.0e-8:
        raise AssertionError(
            "chemistry activity is missing: "
            f"species={species_change}, temperature={temperature_change}"
        )
    print(
        "Reactive 3D AMR chemistry: PASS "
        f"(euler={euler_error:.3e}, elements={element_error:.3e}, "
        f"species_change={species_change:.3e}, "
        f"temperature_change={temperature_change:.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
