#!/usr/bin/env python3
"""Validate inert/reacting chemistry splits against the frozen 0.214 StateRedist run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re


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
ELEMENT_ATOMS = {
    "H2": (2.0, 0.0, 0.0),
    "H": (1.0, 0.0, 0.0),
    "O": (0.0, 1.0, 0.0),
    "O2": (0.0, 2.0, 0.0),
    "OH": (1.0, 1.0, 0.0),
    "H2O": (2.0, 1.0, 0.0),
    "N2": (0.0, 0.0, 2.0),
}
EXPECTED_SPECIES_CHANGE = 0.05331900521091094
EXPECTED_REACTIVE_MAXIMUM_DENSITY = 0.3525627952434091
EXPECTED_REACTIVE_SHA256 = (
    "bacc446c59ba455cd88c26e605c469e8fd38e69e1c56e7a9c3ff18849138bd86"
)
EXPECTED_STEPS = 2
EXPECTED_FINAL_TIME = 5.0e-5
EXPECTED_MINIMUM_DT = 2.2731185428659142e-5
NASA7_COEFFICIENTS = {
    "H2": (
        (2.34433112, 7.98052075e-3, -1.94781510e-5, 2.01572094e-8,
         -7.37611761e-12, -917.935173, 0.683010238),
        (3.33727920, -4.94024731e-5, 4.99456778e-7, -1.79566394e-10,
         2.00255376e-14, -950.158922, -3.20502331),
    ),
    "H": (
        (2.5, 7.05332819e-13, -1.99591964e-15, 2.30081632e-18,
         -9.27732332e-22, 2.54736599e4, -0.446682853),
        (2.50000001, -2.30842973e-11, 1.61561948e-14, -4.73515235e-18,
         4.98197357e-22, 2.54736599e4, -0.446682914),
    ),
    "O": (
        (3.16826710, -3.27931884e-3, 6.64306396e-6, -6.12806624e-9,
         2.11265971e-12, 2.91222592e4, 2.05193346),
        (2.56942078, -8.59741137e-5, 4.19484589e-8, -1.00177799e-11,
         1.22833691e-15, 2.92175791e4, 4.78433864),
    ),
    "O2": (
        (3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9,
         3.24372837e-12, -1063.94356, 3.65767573),
        (3.28253784, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10,
         -2.16717794e-14, -1088.45772, 5.45323129),
    ),
    "OH": (
        (3.99201543, -2.40131752e-3, 4.61793841e-6, -3.88113333e-9,
         1.36411470e-12, 3615.08056, -0.103925458),
        (3.09288767, 5.48429716e-4, 1.26505228e-7, -8.79461556e-11,
         1.17412376e-14, 3858.65700, 4.47669610),
    ),
    "H2O": (
        (4.19864056, -2.03643410e-3, 6.52040211e-6, -5.48797062e-9,
         1.77197817e-12, -3.02937267e4, -0.849032208),
        (3.03399249, 2.17691804e-3, -1.64072518e-7, -9.70419870e-11,
         1.68200992e-14, -3.00042971e4, 4.96677010),
    ),
    "N2": (
        (3.29867700, 1.40824040e-3, -3.96322200e-6, 5.64151500e-9,
         -2.44485400e-12, -1020.89990, 3.95037200),
        (2.92664000, 1.48797680e-3, -5.68476000e-7, 1.00970380e-10,
         -6.75335100e-15, -922.797700, 5.98052800),
    ),
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


def load(path: Path) -> tuple[list[dict[str, float]], list[bytes]]:
    raw_lines = path.read_bytes().splitlines()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise AssertionError(f"{path}: unexpected columns {reader.fieldnames}")
        rows = [{name: float(row[name]) for name in COLUMNS} for row in reader]
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError(f"{path}: non-finite CSV value")
    if len(raw_lines) != len(rows) + 1:
        raise AssertionError(f"{path}: unexpected physical line count")
    return rows, raw_lines[1:]


def specific_internal_energy(name: str, temperature: float) -> float:
    low, high = NASA7_COEFFICIENTS[name]
    coefficients = low if temperature <= 1000.0 else high
    t2 = temperature * temperature
    t3 = t2 * temperature
    t4 = t2 * t2
    h_over_rt = (
        coefficients[0]
        + 0.5 * coefficients[1] * temperature
        + coefficients[2] * t2 / 3.0
        + 0.25 * coefficients[3] * t3
        + coefficients[4] * t4 / 5.0
        + coefficients[5] / temperature
    )
    return 8.31446261815324e3 * temperature * (h_over_rt - 1.0) / (
        MOLECULAR_WEIGHTS[name]
    )


def summary_value(text: str, label: str) -> float:
    pattern = rf"^{re.escape(label)}\s*([+\-0-9.eEdD]+)\s*$"
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise AssertionError(f"summary has {len(matches)} values for {label}")
    return float(matches[0].replace("D", "E").replace("d", "e"))


def check_summary(path: Path, chemistry_enabled: bool) -> None:
    text = path.read_text(encoding="utf-8")
    expected_flag = "T" if chemistry_enabled else "F"
    if not re.search(rf"^Chemistry:\s+{expected_flag}\s*$", text, re.MULTILINE):
        raise AssertionError(f"{path}: unexpected chemistry flag")
    if not re.search(r"^Redistribution:\s+state_redist\s*$", text, re.MULTILINE):
        raise AssertionError(f"{path}: unexpected redistribution")
    steps = summary_value(text, "Completed steps:")
    final_time = summary_value(text, "Final time:")
    minimum_dt = summary_value(text, "Minimum accepted dt:")
    if steps != EXPECTED_STEPS:
        raise AssertionError(f"{path}: expected {EXPECTED_STEPS} steps, found {steps}")
    if abs(final_time - EXPECTED_FINAL_TIME) > 5.0e-18:
        raise AssertionError(f"{path}: final time drift {final_time}")
    if abs(minimum_dt - EXPECTED_MINIMUM_DT) > 5.0e-18:
        raise AssertionError(f"{path}: minimum timestep drift {minimum_dt}")


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


def structure(
    path: Path,
    rows: list[dict[str, float]],
    expected_mass_fractions: dict[str, float] | None = None,
) -> tuple[dict[str, float], tuple[float, float, float]]:
    nx, ny, nz = 10, 8, 6
    if len(rows) != nx * ny * nz:
        raise AssertionError(f"{path}: expected 480 rows, found {len(rows)}")
    counts = {0: 0, 1: 0, 2: 0}
    cell_volume = 0.1 * 0.125 / 6.0
    gas_constant = 8.31446261815324e3
    integrals = {name: 0.0 for name in ("rho", "rhov", "rhow", "rhoE")}
    element_integrals = [0.0, 0.0, 0.0]
    closure_error = 0.0
    relationship_error = 0.0
    energy_error = 0.0
    composition_error = 0.0
    active_density: list[float] = []
    active_pressure: list[float] = []
    for index, row in enumerate(rows):
        i = index % nx + 1
        j = (index // nx) % ny + 1
        k = index // (nx * ny) + 1
        if tuple(round(row[name]) for name in ("i", "j", "k")) != (i, j, k):
            raise AssertionError(f"{path}: row {index} is not x-fastest")
        center = ((i - 0.5) / nx, (j - 0.5) / ny, (k - 0.5) / nz)
        if max(
            abs(row[name] - expected)
            for name, expected in zip(("x", "y", "z"), center)
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
        mass_fraction_sum = 0.0
        species_density_sum = 0.0
        inverse_molecular_weight = 0.0
        for name in SPECIES:
            mass_fraction = row[f"Y_{name}"]
            species_density = row[f"rhoY_{name}"]
            if mass_fraction < 0.0 or species_density < 0.0:
                raise AssertionError(f"{path}: negative species state at row {index}")
            mass_fraction_sum += mass_fraction
            species_density_sum += species_density
            inverse_molecular_weight += mass_fraction / MOLECULAR_WEIGHTS[name]
            if expected_mass_fractions is not None:
                composition_error = max(
                    composition_error,
                    abs(mass_fraction - expected_mass_fractions[name]),
                )
            composition_error = max(
                composition_error,
                relative_error(species_density, row["rho"] * mass_fraction),
            )
            atoms = ELEMENT_ATOMS[name]
            for element in range(3):
                element_integrals[element] += (
                    row["volume_fraction"]
                    * cell_volume
                    * species_density
                    * atoms[element]
                    / MOLECULAR_WEIGHTS[name]
                )
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )
        if inverse_molecular_weight <= 0.0:
            raise AssertionError(f"{path}: invalid mixture at row {index}")
        mixture_molecular_weight = 1.0 / inverse_molecular_weight
        mixture_gas_constant = gas_constant / mixture_molecular_weight
        relationship_error = max(
            relationship_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
            relative_error(
                row["pressure"], row["rho"] * mixture_gas_constant * row["temperature"]
            ),
        )
        mixture_internal_energy = sum(
            row[f"Y_{name}"]
            * specific_internal_energy(name, row["temperature"])
            for name in SPECIES
        )
        expected_total_energy = row["rho"] * (
            mixture_internal_energy
            + 0.5 * (row["u"] ** 2 + row["v"] ** 2 + row["w"] ** 2)
        )
        energy_error = max(
            energy_error, relative_error(row["rhoE"], expected_total_energy)
        )
        for name in integrals:
            integrals[name] += row["volume_fraction"] * row[name] * cell_volume
        if cell_type != 0:
            active_density.append(row["rho"])
            active_pressure.append(row["pressure"])

    if counts != {0: 144, 1: 48, 2: 288}:
        raise AssertionError(f"{path}: geometry counts {counts}")
    if relationship_error > 3.0e-12:
        raise AssertionError(f"{path}: EOS relationship error {relationship_error}")
    if energy_error > 3.0e-12:
        raise AssertionError(f"{path}: local total-energy error {energy_error}")
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
    integrals["energy_error"] = energy_error
    return integrals, tuple(element_integrals)


def inert_mass_fractions() -> dict[str, float]:
    mean_molecular_weight = sum(
        MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name] for name in SPECIES
    )
    return {
        name: MOLE_FRACTIONS[name] * MOLECULAR_WEIGHTS[name] / mean_molecular_weight
        for name in SPECIES
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--inert", type=Path, required=True)
    parser.add_argument("--reactive", type=Path, required=True)
    parser.add_argument("--inert-log", type=Path, required=True)
    parser.add_argument("--reactive-log", type=Path, required=True)
    args = parser.parse_args()

    baseline_rows, _ = load(args.baseline)
    inert_rows, inert_raw_rows = load(args.inert)
    reactive_rows, reactive_raw_rows = load(args.reactive)
    check_summary(args.inert_log, False)
    check_summary(args.reactive_log, True)
    if args.baseline.read_bytes() != args.inert.read_bytes():
        raise AssertionError("frozen StateRedist baseline and chemistry-inert output differ")
    reactive_hash = hashlib.sha256(args.reactive.read_bytes()).hexdigest()
    if reactive_hash != EXPECTED_REACTIVE_SHA256:
        raise AssertionError(f"reactive full-field signature drift: {reactive_hash}")
    inert_expected = inert_mass_fractions()
    baseline, baseline_elements = structure(
        args.baseline, baseline_rows, inert_expected
    )
    inert, inert_elements = structure(args.inert, inert_rows, inert_expected)
    reactive, reactive_elements = structure(args.reactive, reactive_rows)

    conservation_error = max(
        relative_error(reactive[name], baseline[name])
        for name in ("rho", "rhov", "rhow", "rhoE")
    )
    if conservation_error > 5.0e-9:
        raise AssertionError(f"reactive mass/energy/tangential drift {conservation_error}")
    element_error = max(
        abs(reactive_elements[index] - inert_elements[index])
        / max(1.0e-30, abs(inert_elements[index]))
        for index in range(3)
    )
    if element_error > 5.0e-9:
        raise AssertionError(f"reactive elemental drift {element_error}")
    if max(
        abs(reactive_elements[index] - baseline_elements[index])
        / max(1.0e-30, abs(baseline_elements[index]))
        for index in range(3)
    ) > 5.0e-9:
        raise AssertionError("reactive elemental totals do not match the frozen baseline")

    covered_count = 0
    species_changes = {name: 0.0 for name in SPECIES}
    for index, (inert_row, reactive_row) in enumerate(zip(inert_rows, reactive_rows)):
        if round(inert_row["cell_type"]) == 0:
            covered_count += 1
            if inert_raw_rows[index] != reactive_raw_rows[index]:
                raise AssertionError(f"covered row {index} changed during chemistry split")
            continue
        for name in SPECIES:
            species_changes[name] = max(
                species_changes[name],
                abs(reactive_row[f"rhoY_{name}"] - inert_row[f"rhoY_{name}"]),
            )
    chemistry_change = max(species_changes.values())
    if chemistry_change <= 1.0e-10:
        raise AssertionError(f"active species chemistry effect is unresolved: {species_changes}")
    if abs(chemistry_change - EXPECTED_SPECIES_CHANGE) > 2.0e-10:
        raise AssertionError(
            f"reactive species-change signature drift: {chemistry_change}"
        )
    if (
        abs(reactive["maximum_density"] - EXPECTED_REACTIVE_MAXIMUM_DENSITY)
        > 5.0e-11
    ):
        raise AssertionError(
            f"reactive maximum-density signature drift: {reactive['maximum_density']}"
        )
    print(
        "Reactive planar EB 3D chemistry: PASS "
        f"(rows={len(reactive_rows)}, covered={covered_count}, "
        f"species_change={chemistry_change:.10e}, "
        f"conservation_error={conservation_error:.3e}, "
        f"element_error={element_error:.3e}, "
        f"energy_error={reactive['energy_error']:.3e}, "
        f"rho_max={reactive['maximum_density']:.10e}, "
        f"sha256={reactive_hash[:12]})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
