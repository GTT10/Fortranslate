#!/usr/bin/env python3
"""Compare the qualified PeleF dilute-gas transport subset with Cantera 3.2."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")


def relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-300)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--viscosity-tolerance",
        type=float,
        default=0.05,
        help="relative qualification bound for mixture viscosity (default: 5%%)",
    )
    parser.add_argument(
        "--conductivity-tolerance",
        type=float,
        default=0.10,
        help="relative qualification bound for thermal conductivity (default: 10%%)",
    )
    parser.add_argument(
        "--diffusion-tolerance",
        type=float,
        default=0.40,
        help=(
            "relative qualification bound for mixture-averaged diffusion "
            "without polar/full-fit corrections (default: 40%%)"
        ),
    )
    args = parser.parse_args()

    try:
        import cantera as ct
    except ImportError as exc:  # pragma: no cover - exercised in CI only
        raise SystemExit("Cantera is required for the transport comparison") from exc

    gas = ct.Solution("h2o2.yaml", "ohmech")
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 4:
        raise AssertionError("transport probe must contain four reference states")

    maximum_mu_error = 0.0
    maximum_lambda_error = 0.0
    maximum_diffusion_error = 0.0
    worst_mu = ""
    worst_lambda = ""
    worst_diffusion = ""

    for row in rows:
        label = row["label"].strip()
        mass_fractions = {name: float(row[f"Y_{name}"]) for name in SPECIES}
        gas.TPY = float(row["temperature"]), float(row["pressure"]), mass_fractions

        mu_error = relative_error(float(row["viscosity"]), gas.viscosity)
        lambda_error = relative_error(
            float(row["thermal_conductivity"]), gas.thermal_conductivity
        )
        if mu_error > maximum_mu_error:
            maximum_mu_error = mu_error
            worst_mu = label
        if lambda_error > maximum_lambda_error:
            maximum_lambda_error = lambda_error
            worst_lambda = label

        cantera_diffusion = gas.mix_diff_coeffs_mass
        for name in SPECIES:
            index = gas.species_index(name)
            error = relative_error(float(row[f"D_{name}"]), cantera_diffusion[index])
            if error > maximum_diffusion_error:
                maximum_diffusion_error = error
                worst_diffusion = f"{label}:{name}"

    print(f"maximum_viscosity_relative_error={maximum_mu_error:.16e} ({worst_mu})")
    print(
        "maximum_conductivity_relative_error="
        f"{maximum_lambda_error:.16e} ({worst_lambda})"
    )
    print(
        "maximum_diffusion_relative_error="
        f"{maximum_diffusion_error:.16e} ({worst_diffusion})"
    )
    print(
        "qualification_bounds="
        f"mu:{args.viscosity_tolerance:.6g},"
        f"lambda:{args.conductivity_tolerance:.6g},"
        f"diffusion:{args.diffusion_tolerance:.6g}"
    )

    if maximum_mu_error > args.viscosity_tolerance:
        raise AssertionError("dilute-gas mixture viscosity is outside qualification")
    if maximum_lambda_error > args.conductivity_tolerance:
        raise AssertionError("modified-Eucken conductivity is outside qualification")
    if maximum_diffusion_error > args.diffusion_tolerance:
        raise AssertionError("mixture-averaged diffusion is outside qualification")
    print("Transport comparison against Cantera: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
