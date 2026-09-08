#!/usr/bin/env python3
"""Validate planar 3D EB transport and coupled R-T-H-T-R outputs."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
import re

from check_reactive_eb_chemistry_3d import (
    COLUMNS,
    ELEMENT_ATOMS,
    MOLECULAR_WEIGHTS,
    MOLE_FRACTIONS,
    SPECIES,
    load,
    relative_error,
    specific_internal_energy,
)


NX, NY, NZ = 10, 8, 6
LENGTH = 1.0e-4
FINAL_TIME = 2.0e-9
EXPECTED_STEPS = 1
EXPECTED_MINIMUM_DT = 2.0e-9
EXPECTED_DIFFUSIVITY = 9.7180419287635654e-4
EXPECTED_THETA = 1.0
EXPECTED_TRANSPORT_RESPONSE = 0.10489647406832604
EXPECTED_CHEMISTRY_RESPONSE = 4.043407673691575e-8
EXPECTED_CONTROL_SHA256 = (
    "98618a9270fc561887e4d2c8dd40849f66d9d17915c6caef74c7453b16d94a91",
)
EXPECTED_TRANSPORT_SHA256 = (
    "f525d24be8bb0c59ef25648474cbf23179c6a5e3b0eb27e3ffd2df091dc0de73",
    "b976f8b94f7d2afe2734aefc1ce6858f8f22b1608807b8008339418dfe682e14",
)
EXPECTED_COUPLED_SHA256 = (
    "7f0f83b35e4fe15fc3648f871184f1aaad5f9e6412f779b2232b1cf7f8d274aa",
    "f2c0ac5a14509dcd5b6d40673c8211f97fdcfb0355211f9fa000c218b16c2f76",
)


def summary_value(text: str, label: str) -> float:
    pattern = rf"^{re.escape(label)}\s*([+\-0-9.eEdD]+)\s*$"
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise AssertionError(f"summary has {len(matches)} values for {label}")
    return float(matches[0].replace("D", "E").replace("d", "e"))


def check_summary(
    path: Path, chemistry: bool, transport: bool, sequence: str
) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    expected_chemistry = "T" if chemistry else "F"
    expected_transport = "T" if transport else "F"
    patterns = (
        rf"^Chemistry:\s+{expected_chemistry}\s*$",
        rf"^Molecular transport:\s+{expected_transport}\s*$",
        r"^Redistribution:\s+state_redist\s*$",
        r"^Viscosity:\s+T\s*$",
        r"^Thermal conduction:\s+T\s*$",
        r"^Species diffusion:\s+T\s*$",
        r"^Barodiffusion:\s+T\s*$",
        rf"^Operator sequence:\s+{re.escape(sequence)}\s*$",
    )
    if any(not re.search(pattern, text, re.MULTILINE) for pattern in patterns):
        raise AssertionError(f"{path}: unexpected physics summary")
    steps = summary_value(text, "Completed steps:")
    final_time = summary_value(text, "Final time:")
    minimum_dt = summary_value(text, "Minimum accepted dt:")
    if steps != EXPECTED_STEPS or abs(final_time - FINAL_TIME) > 5.0e-21:
        raise AssertionError(f"{path}: unexpected step schedule")
    if abs(minimum_dt - EXPECTED_MINIMUM_DT) > 5.0e-21:
        raise AssertionError(f"{path}: invalid minimum timestep {minimum_dt}")
    transport_cfl = summary_value(text, "Transport CFL:")
    if abs(transport_cfl - 0.35) > 5.0e-15:
        raise AssertionError(f"{path}: unexpected transport CFL {transport_cfl}")
    maximum_diffusivity = summary_value(text, "Maximum transport diffusivity:")
    minimum_theta = summary_value(text, "Minimum transport flux theta:")
    invariant_error = summary_value(
        text, "Maximum invariant conservation error:"
    )
    element_error = summary_value(
        text, "Maximum elemental conservation error:"
    )
    if not all(math.isfinite(value) and 0.0 <= value <= 5.0e-10 for value in (
        invariant_error, element_error
    )):
        raise AssertionError(f"{path}: invalid conservation diagnostics")
    if not chemistry and element_error != 0.0:
        raise AssertionError(f"{path}: disabled chemistry changed elements")
    if transport:
        if not math.isclose(
            maximum_diffusivity, EXPECTED_DIFFUSIVITY, rel_tol=2.0e-13
        ) or not math.isclose(minimum_theta, EXPECTED_THETA, abs_tol=1.0e-15):
            raise AssertionError(f"{path}: invalid transport diagnostics")
    elif maximum_diffusivity != 0.0 or minimum_theta != 1.0:
        raise AssertionError(f"{path}: disabled transport diagnostics changed")
    return {
        "steps": steps,
        "minimum_dt": minimum_dt,
        "maximum_diffusivity": maximum_diffusivity,
        "minimum_theta": minimum_theta,
        "invariant_error": invariant_error,
        "element_error": element_error,
    }


def expected_geometry(
    i: int, j: int, k: int
) -> tuple[int, float, tuple[float, float, float]]:
    dx, dy, dz = LENGTH / NX, LENGTH / NY, LENGTH / NZ
    x = (i - 0.5) * dx
    y = (j - 0.5) * dy
    z = (k - 0.5) * dz
    if i <= 3:
        return 0, 0.0, (x, y, z)
    if i == 4:
        return 1, 0.05, (3.975e-5, y, z)
    return 2, 1.0, (x, y, z)


def inspect_rows(path: Path, rows: list[dict[str, float]]) -> dict[str, object]:
    if len(rows) != NX * NY * NZ:
        raise AssertionError(f"{path}: expected 480 rows, found {len(rows)}")
    counts = {0: 0, 1: 0, 2: 0}
    cell_volume = (LENGTH / NX) * (LENGTH / NY) * (LENGTH / NZ)
    integral_names = ("rho", "rhov", "rhow", "rhoE", *(
        f"rhoY_{name}" for name in SPECIES
    ))
    integrals = {name: 0.0 for name in integral_names}
    l1_integrals = {name: 0.0 for name in integral_names}
    elements = [0.0, 0.0, 0.0]
    relationship_error = 0.0
    closure_error = 0.0
    energy_error = 0.0
    active_indices: list[int] = []
    for index, row in enumerate(rows):
        i = index % NX + 1
        j = (index // NX) % NY + 1
        k = index // (NX * NY) + 1
        if tuple(round(row[name]) for name in ("i", "j", "k")) != (i, j, k):
            raise AssertionError(f"{path}: row {index} is not x-fastest")
        if abs(row["time"] - FINAL_TIME) > 5.0e-21:
            raise AssertionError(f"{path}: unexpected time at row {index}")
        cell_type, kappa, centroid = expected_geometry(i, j, k)
        center = ((i - 0.5) * LENGTH / NX, (j - 0.5) * LENGTH / NY,
                  (k - 0.5) * LENGTH / NZ)
        if round(row["cell_type"]) != cell_type:
            raise AssertionError(f"{path}: bad cell type at row {index}")
        if abs(row["volume_fraction"] - kappa) > 3.0e-14:
            raise AssertionError(f"{path}: bad volume fraction at row {index}")
        if max(abs(row[name] - value) for name, value in zip(
            ("x", "y", "z"), center
        )) > 3.0e-18:
            raise AssertionError(f"{path}: bad Cartesian center at row {index}")
        if max(abs(row[name] - value) for name, value in zip(
            ("fluid_centroid_x", "fluid_centroid_y", "fluid_centroid_z"),
            centroid,
        )) > 3.0e-18:
            raise AssertionError(f"{path}: bad fluid centroid at row {index}")
        counts[cell_type] += 1
        if min(row["rho"], row["pressure"], row["temperature"]) <= 0.0:
            raise AssertionError(f"{path}: nonphysical state at row {index}")
        mass_sum = 0.0
        density_sum = 0.0
        inverse_mw = 0.0
        for name in SPECIES:
            mass_fraction = row[f"Y_{name}"]
            species_density = row[f"rhoY_{name}"]
            if mass_fraction < 0.0 or species_density < 0.0:
                raise AssertionError(f"{path}: negative species at row {index}")
            mass_sum += mass_fraction
            density_sum += species_density
            inverse_mw += mass_fraction / MOLECULAR_WEIGHTS[name]
            relationship_error = max(
                relationship_error,
                relative_error(species_density, row["rho"] * mass_fraction),
            )
            for element, atoms in enumerate(ELEMENT_ATOMS[name]):
                elements[element] += (
                    kappa * cell_volume * species_density * atoms
                    / MOLECULAR_WEIGHTS[name]
                )
        closure_error = max(
            closure_error,
            abs(mass_sum - 1.0),
            relative_error(density_sum, row["rho"]),
        )
        mixture_mw = 1.0 / inverse_mw
        relationship_error = max(
            relationship_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
            relative_error(
                row["pressure"],
                row["rho"] * 8.31446261815324e3 / mixture_mw
                * row["temperature"],
            ),
        )
        internal_energy = sum(
            row[f"Y_{name}"]
            * specific_internal_energy(name, row["temperature"])
            for name in SPECIES
        )
        expected_energy = row["rho"] * (
            internal_energy
            + 0.5 * (row["u"] ** 2 + row["v"] ** 2 + row["w"] ** 2)
        )
        energy_error = max(
            energy_error, relative_error(row["rhoE"], expected_energy)
        )
        for name in integral_names:
            contribution = kappa * cell_volume * row[name]
            integrals[name] += contribution
            l1_integrals[name] += abs(contribution)
        if cell_type != 0:
            active_indices.append(index)
    if counts != {0: 144, 1: 48, 2: 288}:
        raise AssertionError(f"{path}: geometry counts {counts}")
    if max(relationship_error, closure_error, energy_error) > 3.0e-12:
        raise AssertionError(
            f"{path}: state contract errors relation={relationship_error}, "
            f"closure={closure_error}, energy={energy_error}"
        )
    return {
        "integrals": integrals,
        "l1_integrals": l1_integrals,
        "elements": tuple(elements),
        "active_indices": active_indices,
        "energy_error": energy_error,
    }


def maximum_relative_difference(
    left: list[dict[str, float]],
    right: list[dict[str, float]],
    indices: list[int],
    fields: tuple[str, ...],
) -> float:
    return max(
        relative_error(left[index][field], right[index][field])
        for index in indices
        for field in fields
    )


def compare_integrals(
    actual: dict[str, float], expected: dict[str, float],
    actual_l1: dict[str, float], expected_l1: dict[str, float],
    names: tuple[str, ...],
) -> float:
    return max(
        abs(actual[name] - expected[name])
        / max(1.0e-300, abs(expected[name]), actual_l1[name], expected_l1[name])
        for name in names
    )


def compare_elements(actual: tuple[float, ...], expected: tuple[float, ...]) -> float:
    return max(
        abs(left - right) / max(1.0e-30, abs(right))
        for left, right in zip(actual, expected)
    )


def check_hash(path: Path, expected: tuple[str, ...]) -> str:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest not in expected:
        raise AssertionError(f"{path}: SHA-256 drift {digest}")
    return digest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--transport", type=Path, required=True)
    parser.add_argument("--coupled", type=Path, required=True)
    parser.add_argument("--control-log", type=Path, required=True)
    parser.add_argument("--transport-log", type=Path, required=True)
    parser.add_argument("--coupled-log", type=Path, required=True)
    args = parser.parse_args()

    control, control_raw = load(args.control)
    transport, transport_raw = load(args.transport)
    coupled, coupled_raw = load(args.coupled)
    control_info = inspect_rows(args.control, control)
    transport_info = inspect_rows(args.transport, transport)
    coupled_info = inspect_rows(args.coupled, coupled)
    control_schedule = check_summary(args.control_log, False, False, "H")
    transport_schedule = check_summary(args.transport_log, False, True, "T-H-T")
    coupled_schedule = check_summary(
        args.coupled_log, True, True, "R-T-H-T-R"
    )

    active = control_info["active_indices"]
    assert isinstance(active, list)
    for index, row in enumerate(control):
        if round(row["cell_type"]) != 0:
            continue
        if control_raw[index] != transport_raw[index] or (
            control_raw[index] != coupled_raw[index]
        ):
            raise AssertionError(f"covered row {index} changed")

    conserved = ("rho", "rhov", "rhow", "rhoE")
    species = tuple(f"rhoY_{name}" for name in SPECIES)
    transport_error = compare_integrals(
        transport_info["integrals"], control_info["integrals"],
        transport_info["l1_integrals"], control_info["l1_integrals"],
        (*conserved, *species),
    )
    coupled_error = compare_integrals(
        coupled_info["integrals"], control_info["integrals"],
        coupled_info["l1_integrals"], control_info["l1_integrals"], conserved,
    )
    transport_element_error = compare_elements(
        transport_info["elements"], control_info["elements"]
    )
    coupled_element_error = compare_elements(
        coupled_info["elements"], control_info["elements"]
    )
    if max(
        transport_error,
        coupled_error,
        transport_element_error,
        coupled_element_error,
    ) > 5.0e-10:
        raise AssertionError("transport/coupling conservation gate failed")

    transport_response = maximum_relative_difference(
        transport, control, active,
        ("temperature", "rhoE", *(f"rhoY_{name}" for name in SPECIES)),
    )
    chemistry_response = maximum_relative_difference(
        coupled, transport, active, tuple(f"rhoY_{name}" for name in SPECIES)
    )
    if transport_response <= 1.0e-10:
        raise AssertionError("public transport case has no resolved response")
    if chemistry_response <= 1.0e-12:
        raise AssertionError("public coupled case has no chemistry response")
    if not math.isclose(
        transport_response, EXPECTED_TRANSPORT_RESPONSE, rel_tol=2.0e-12
    ):
        raise AssertionError("public transport response drifted")
    if not math.isclose(
        chemistry_response, EXPECTED_CHEMISTRY_RESPONSE, rel_tol=2.0e-11
    ):
        raise AssertionError("public chemistry response drifted")

    control_hash = check_hash(args.control, EXPECTED_CONTROL_SHA256)
    transport_hash = check_hash(args.transport, EXPECTED_TRANSPORT_SHA256)
    coupled_hash = check_hash(args.coupled, EXPECTED_COUPLED_SHA256)
    maximum_energy_error = max(
        float(control_info["energy_error"]),
        float(transport_info["energy_error"]),
        float(coupled_info["energy_error"]),
    )
    print(
        "Reactive planar EB 3D transport: PASS "
        f"(rows={len(control)}, covered=144, "
        f"transport_response={transport_response:.3e}, "
        f"chemistry_response={chemistry_response:.3e}, "
        f"conservation_error={max(transport_error, coupled_error):.3e}, "
        f"element_error={max(transport_element_error, coupled_element_error):.3e}, "
        f"energy_error={maximum_energy_error:.3e}, "
        f"steps={int(transport_schedule['steps'])}, "
        f"minimum_dt={transport_schedule['minimum_dt']:.3e}, "
        f"theta={transport_schedule['minimum_theta']:.3e}, "
        f"diffusivity={transport_schedule['maximum_diffusivity']:.3e}, "
        f"sha256={control_hash[:12]}/{transport_hash[:12]}/{coupled_hash[:12]}, "
        f"coupled_steps={int(coupled_schedule['steps'])}, "
        f"control_steps={int(control_schedule['steps'])})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
