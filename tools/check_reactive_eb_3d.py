#!/usr/bin/env python3
"""Validate the public planar 3D EB StateRedist and FluxRedist CSV files."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")
MOLE_FRACTIONS = {
    "H2": 0.29570,
    "H": 1.0e-5,
    "O": 1.0e-5,
    "O2": 0.14784,
    "OH": 1.0e-5,
    "H2O": 0.0,
    "N2": 0.55643,
}
MOLECULAR_WEIGHTS = {
    "H2": 2.016,
    "H": 1.008,
    "O": 15.999,
    "O2": 31.998,
    "OH": 17.007,
    "H2O": 18.015,
    "N2": 28.014,
}
BASE_COLUMNS = (
    "time",
    "i",
    "j",
    "k",
    "x",
    "y",
    "z",
    "cell_type",
    "volume_fraction",
    "fluid_centroid_x",
    "fluid_centroid_y",
    "fluid_centroid_z",
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


def load(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise AssertionError(f"{path}: unexpected columns {reader.fieldnames}")
        rows = [{name: float(row[name]) for name in COLUMNS} for row in reader]
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError(f"{path}: non-finite CSV value")
    return rows


def expected_geometry(i: int, j: int, k: int) -> tuple[int, float, tuple[float, ...]]:
    dx, dy, dz = 0.1, 0.125, 1.0 / 6.0
    x = (i - 0.5) * dx
    y = (j - 0.5) * dy
    z = (k - 0.5) * dz
    if i <= 3:
        return 0, 0.0, (x, y, z)
    if i == 4:
        return 1, 0.05, (0.3975, y, z)
    return 2, 1.0, (x, y, z)


def validate_rows(path: Path, rows: list[dict[str, float]]) -> dict[str, float]:
    nx, ny, nz = 10, 8, 6
    if len(rows) != nx * ny * nz:
        raise AssertionError(f"{path}: expected 480 rows, found {len(rows)}")
    counts = {0: 0, 1: 0, 2: 0}
    mean_molecular_weight = sum(
        MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name] for name in SPECIES
    )
    expected_y = {
        name: MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name]
        / mean_molecular_weight
        for name in SPECIES
    }
    gas_constant = 8.31446261815324e3 / mean_molecular_weight
    cell_volume = 0.1 * 0.125 / 6.0
    integrals = {
        name: 0.0
        for name in ("rho", "rhou", "rhov", "rhow", "rhoE")
    }
    integrals.update({f"rhoY_{name}": 0.0 for name in SPECIES})
    closure_error = 0.0
    relationship_error = 0.0
    composition_error = 0.0
    active_density: list[float] = []
    active_pressure: list[float] = []
    for index, row in enumerate(rows):
        i = index % nx + 1
        j = (index // nx) % ny + 1
        k = index // (nx * ny) + 1
        if tuple(round(row[name]) for name in ("i", "j", "k")) != (i, j, k):
            raise AssertionError(f"{path}: row {index} is not x-fastest")
        expected_center = ((i - 0.5) / nx, (j - 0.5) / ny, (k - 0.5) / nz)
        if max(
            abs(row[name] - expected)
            for name, expected in zip(("x", "y", "z"), expected_center)
        ) > 3.0e-15:
            raise AssertionError(f"{path}: bad Cartesian center at row {index}")
        cell_type, kappa, centroid = expected_geometry(i, j, k)
        if round(row["cell_type"]) != cell_type:
            raise AssertionError(f"{path}: bad cell type at row {index}")
        if abs(row["volume_fraction"] - kappa) > 3.0e-14:
            raise AssertionError(f"{path}: bad volume fraction at row {index}")
        if max(
            abs(row[name] - expected)
            for name, expected in zip(
                ("fluid_centroid_x", "fluid_centroid_y", "fluid_centroid_z"),
                centroid,
            )
        ) > 3.0e-14:
            raise AssertionError(f"{path}: bad fluid centroid at row {index}")
        counts[cell_type] += 1
        if abs(row["time"] - 5.0e-5) > 5.0e-18:
            raise AssertionError(f"{path}: unexpected output time")
        if row["rho"] <= 0.0 or row["pressure"] <= 0.0 or row["temperature"] <= 0.0:
            raise AssertionError(f"{path}: nonphysical state at row {index}")
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
            mass_fraction_sum += mass_fraction
            species_density_sum += species_density
            composition_error = max(
                composition_error,
                abs(mass_fraction - expected_y[name]),
                relative_error(species_density, row["rho"] * mass_fraction),
            )
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )
        for name in integrals:
            integrals[name] += kappa * row[name] * cell_volume
        if cell_type != 0:
            active_density.append(row["rho"])
            active_pressure.append(row["pressure"])

    if counts != {0: 144, 1: 48, 2: 288}:
        raise AssertionError(f"{path}: geometry counts {counts}")
    expected_mass = 0.005 * 0.93 + 0.6 * 0.31
    if relative_error(integrals["rho"], expected_mass) > 3.0e-13:
        raise AssertionError(f"{path}: mass integral {integrals['rho']}")
    for name in SPECIES:
        expected_species_mass = expected_mass * expected_y[name]
        if relative_error(integrals[f"rhoY_{name}"], expected_species_mass) > 3.0e-13:
            raise AssertionError(f"{path}: {name} integral mismatch")
    if max(abs(integrals["rhov"]), abs(integrals["rhow"])) > 3.0e-13:
        raise AssertionError(f"{path}: tangential momentum is not conserved")
    if relationship_error > 3.0e-12:
        raise AssertionError(f"{path}: EOS relationship error {relationship_error}")
    if composition_error > 3.0e-12 or closure_error > 3.0e-12:
        raise AssertionError(
            f"{path}: composition={composition_error}, closure={closure_error}"
        )
    if max(active_density) - min(active_density) < 2.0e-2:
        raise AssertionError(f"{path}: density response is unresolved")
    if max(active_pressure) - min(active_pressure) < 5.0e2:
        raise AssertionError(f"{path}: pressure response is unresolved")
    integrals["maximum_density"] = max(active_density)
    integrals["pressure_span"] = max(active_pressure) - min(active_pressure)
    integrals["closure_error"] = closure_error
    return integrals


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--flux", type=Path, required=True)
    args = parser.parse_args()

    state = validate_rows(args.state, load(args.state))
    flux = validate_rows(args.flux, load(args.flux))
    for name in ("rho", "rhov", "rhow", "rhoE", *(f"rhoY_{s}" for s in SPECIES)):
        if relative_error(state[name], flux[name]) > 5.0e-13:
            raise AssertionError(f"cross-method invariant mismatch for {name}")
    if abs(state["maximum_density"] - 0.3514845375860745) > 5.0e-12:
        raise AssertionError("StateRedist density signature drift")
    if abs(flux["maximum_density"] - 0.5712733525097978) > 5.0e-12:
        raise AssertionError("FluxRedist density signature drift")
    if not (state["rhou"] < -6.0e-2 and flux["rhou"] > 6.0e-2):
        raise AssertionError("normal wall-response signatures are unresolved")
    print(
        "Reactive planar EB 3D: PASS "
        f"(rows=480, state_rho_max={state['maximum_density']:.10e}, "
        f"flux_rho_max={flux['maximum_density']:.10e}, "
        f"invariant_delta={relative_error(state['rhoE'], flux['rhoE']):.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
