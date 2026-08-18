#!/usr/bin/env python3
"""Generate a deterministic Fortran elementary-mechanism data module."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def fortran_real(value: float) -> str:
    text = f"{value:.12e}"
    return f"{text}_dp"


def validate(data: dict[str, Any]) -> None:
    species = data["species"]
    if not species or len(species) != len(set(species)):
        raise ValueError("species must be non-empty and unique")
    for reaction in data["reactions"]:
        for side in ("reactants", "products"):
            if not reaction[side]:
                raise ValueError(f"{reaction['equation']}: empty {side}")
            for name, coefficient in reaction[side].items():
                if name not in species or float(coefficient) <= 0.0:
                    raise ValueError(f"{reaction['equation']}: invalid {side}")
        rate = reaction["arrhenius"]
        if float(rate["A"]) < 0.0:
            raise ValueError(f"{reaction['equation']}: negative A")


def generate(data: dict[str, Any]) -> str:
    species: list[str] = data["species"]
    reactions: list[dict[str, Any]] = data["reactions"]
    module_name = data["module_name"]
    loader_name = data["loader_name"]
    kernel_name = data["kernel_name"]

    lines: list[str] = [
        f"module {module_name}",
        "  use precision_mod, only: dp",
        "  use nasa7_thermo_mod, only: nasa7_species",
        "  use elementary_kinetics_mod, only: &",
        "    elementary_reaction, elementary_production_rates",
        "  implicit none",
        "  private",
        "",
        f"  integer, parameter, public :: h2o2_nspecies = {len(species)}",
        f"  integer, parameter, public :: h2o2_nreactions = {len(reactions)}",
    ]
    for index, name in enumerate(species, start=1):
        lines.append(
            f"  integer, parameter, public :: h2o2_{name.lower()}_index = {index}"
        )
    lines += [
        "",
        f"  public :: {loader_name}",
        f"  public :: {kernel_name}",
        "",
        "contains",
        "",
        f"  subroutine {loader_name}(reactions, ok)",
        "    type(elementary_reaction), allocatable, intent(out) :: reactions(:)",
        "    logical, intent(out) :: ok",
        "",
        f"    allocate(reactions({len(reactions)}))",
    ]

    for r_index, reaction in enumerate(reactions, start=1):
        lines += [
            "",
            f"    reactions({r_index})%equation = \"{reaction['equation']}\"",
            f"    allocate(reactions({r_index})%reactant_stoich(h2o2_nspecies))",
            f"    allocate(reactions({r_index})%product_stoich(h2o2_nspecies))",
            f"    reactions({r_index})%reactant_stoich = 0.0_dp",
            f"    reactions({r_index})%product_stoich = 0.0_dp",
        ]
        for name, coefficient in reaction["reactants"].items():
            idx = species.index(name) + 1
            lines.append(
                f"    reactions({r_index})%reactant_stoich({idx}) = "
                f"{fortran_real(float(coefficient))}"
            )
        for name, coefficient in reaction["products"].items():
            idx = species.index(name) + 1
            lines.append(
                f"    reactions({r_index})%product_stoich({idx}) = "
                f"{fortran_real(float(coefficient))}"
            )
        rate = reaction["arrhenius"]
        lines += [
            f"    reactions({r_index})%forward_rate%pre_exponential = "
            f"{fortran_real(float(rate['A']))}",
            f"    reactions({r_index})%forward_rate%temperature_exponent = "
            f"{fortran_real(float(rate['b']))}",
            f"    reactions({r_index})%forward_rate%activation_energy = "
            f"{fortran_real(float(rate['Ea']))}",
            f"    reactions({r_index})%reversible = "
            f"{'.true.' if reaction['reversible'] else '.false.'}",
        ]

    lines += [
        "",
        "    ok = .true.",
        f"  end subroutine {loader_name}",
        "",
        f"  subroutine {kernel_name}( &",
        "      species, reactions, temperature, density, mass_fractions, &",
        "      molar_production_rates, ok)",
        "    type(nasa7_species), intent(in) :: species(:)",
        "    type(elementary_reaction), intent(in) :: reactions(:)",
        "    real(dp), intent(in) :: temperature, density, mass_fractions(:)",
        "    real(dp), intent(out) :: molar_production_rates(:)",
        "    logical, intent(out) :: ok",
        "",
        "    ok = size(species) == h2o2_nspecies .and. &",
        "      size(reactions) == h2o2_nreactions",
        "    if (.not. ok) then",
        "      molar_production_rates = 0.0_dp",
        "      return",
        "    end if",
        "    call elementary_production_rates( &",
        "      species, reactions, temperature, density, mass_fractions, &",
        "      molar_production_rates, ok)",
        f"  end subroutine {kernel_name}",
        "",
        f"end module {module_name}",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    validate(data)
    output = generate(data)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
