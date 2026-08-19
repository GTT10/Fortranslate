#!/usr/bin/env python3
"""Generate deterministic Fortran data for elementary/pressure-dependent kinetics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

REACTION_TYPES = {"elementary", "three-body", "falloff"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def fortran_real(value: float) -> str:
    return f"{value:.12e}_dp"


def validate_arrhenius(rate: dict[str, Any], equation: str) -> None:
    if float(rate["A"]) < 0.0:
        raise ValueError(f"{equation}: negative Arrhenius A")
    for key in ("A", "b", "Ea"):
        float(rate[key])


def validate(data: dict[str, Any]) -> None:
    species = data["species"]
    if not species or len(species) != len(set(species)):
        raise ValueError("species must be non-empty and unique")
    if not data.get("reactions"):
        raise ValueError("mechanism must contain reactions")
    for required in ("module_name", "loader_name", "kernel_name"):
        if not data.get(required):
            raise ValueError(f"missing {required}")

    for reaction in data["reactions"]:
        equation = reaction["equation"]
        reaction_type = reaction.get("type", "elementary")
        if reaction_type not in REACTION_TYPES:
            raise ValueError(f"{equation}: unsupported reaction type")
        for side in ("reactants", "products"):
            if not reaction[side]:
                raise ValueError(f"{equation}: empty {side}")
            for name, coefficient in reaction[side].items():
                if name not in species or float(coefficient) <= 0.0:
                    raise ValueError(f"{equation}: invalid {side}")
        if reaction_type == "falloff":
            validate_arrhenius(reaction["low_rate"], equation)
            validate_arrhenius(reaction["high_rate"], equation)
            troe = reaction.get("troe")
            if troe is not None:
                alpha = float(troe["A"])
                if not 0.0 < alpha < 1.0:
                    raise ValueError(f"{equation}: invalid Troe A")
                if float(troe["T3"]) <= 0.0 or float(troe["T1"]) <= 0.0:
                    raise ValueError(f"{equation}: invalid Troe temperatures")
                if float(troe.get("T2", 0.0)) < 0.0:
                    raise ValueError(f"{equation}: invalid Troe T2")
        else:
            validate_arrhenius(reaction["arrhenius"], equation)
        if reaction_type in {"three-body", "falloff"}:
            default_efficiency = float(reaction.get("default_efficiency", 1.0))
            if default_efficiency < 0.0:
                raise ValueError(f"{equation}: negative default efficiency")
            for name, efficiency in reaction.get("efficiencies", {}).items():
                if name not in species or float(efficiency) < 0.0:
                    raise ValueError(f"{equation}: invalid third-body efficiency")


def emit_rate(lines: list[str], target: str, rate: dict[str, Any]) -> None:
    lines += [
        f"    {target}%pre_exponential = {fortran_real(float(rate['A']))}",
        f"    {target}%temperature_exponent = {fortran_real(float(rate['b']))}",
        f"    {target}%activation_energy = {fortran_real(float(rate['Ea']))}",
    ]


def generate(data: dict[str, Any]) -> str:
    species: list[str] = data["species"]
    reactions: list[dict[str, Any]] = data["reactions"]
    module_name = data["module_name"]
    loader_name = data["loader_name"]
    kernel_name = data["kernel_name"]
    jacobian_name = data.get("jacobian_name", f"{kernel_name}_jacobian")
    prefix = data.get("symbol_prefix", "h2o2")

    lines: list[str] = [
        f"module {module_name}",
        "  use precision_mod, only: dp",
        "  use nasa7_thermo_mod, only: nasa7_species",
        "  use elementary_kinetics_mod, only: &",
        "    elementary_reaction, elementary_production_rates, &",
        "    elementary_mass_fraction_jacobian, reaction_kind_elementary, &",
        "    reaction_kind_three_body, reaction_kind_falloff",
        "  implicit none",
        "  private",
        "",
        f"  integer, parameter, public :: {prefix}_nspecies = {len(species)}",
        f"  integer, parameter, public :: {prefix}_nreactions = {len(reactions)}",
    ]
    for index, name in enumerate(species, start=1):
        lines.append(
            f"  integer, parameter, public :: {prefix}_{name.lower()}_index = {index}"
        )
    lines += [
        "",
        f"  public :: {loader_name}",
        f"  public :: {kernel_name}",
        f"  public :: {jacobian_name}",
        "",
        "contains",
        "",
        f"  subroutine {loader_name}(reactions, ok)",
        "    type(elementary_reaction), allocatable, intent(out) :: reactions(:)",
        "    logical, intent(out) :: ok",
        "",
        f"    allocate(reactions({len(reactions)}))",
    ]

    kind_symbol = {
        "elementary": "reaction_kind_elementary",
        "three-body": "reaction_kind_three_body",
        "falloff": "reaction_kind_falloff",
    }
    for r_index, reaction in enumerate(reactions, start=1):
        reaction_type = reaction.get("type", "elementary")
        lines += [
            "",
            f"    reactions({r_index})%equation = \"{reaction['equation']}\"",
            f"    reactions({r_index})%kind = {kind_symbol[reaction_type]}",
            f"    allocate(reactions({r_index})%reactant_stoich({prefix}_nspecies))",
            f"    allocate(reactions({r_index})%product_stoich({prefix}_nspecies))",
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

        if reaction_type == "falloff":
            emit_rate(
                lines,
                f"reactions({r_index})%low_pressure_rate",
                reaction["low_rate"],
            )
            emit_rate(
                lines,
                f"reactions({r_index})%high_pressure_rate",
                reaction["high_rate"],
            )
        else:
            emit_rate(
                lines,
                f"reactions({r_index})%forward_rate",
                reaction["arrhenius"],
            )

        if reaction_type in {"three-body", "falloff"}:
            default_efficiency = float(reaction.get("default_efficiency", 1.0))
            lines += [
                f"    allocate(reactions({r_index})%third_body_efficiencies({prefix}_nspecies))",
                f"    reactions({r_index})%third_body_efficiencies = "
                f"{fortran_real(default_efficiency)}",
            ]
            for name, efficiency in reaction.get("efficiencies", {}).items():
                idx = species.index(name) + 1
                lines.append(
                    f"    reactions({r_index})%third_body_efficiencies({idx}) = "
                    f"{fortran_real(float(efficiency))}"
                )

        troe = reaction.get("troe")
        if troe is not None:
            lines += [
                f"    reactions({r_index})%troe%enabled = .true.",
                f"    reactions({r_index})%troe%alpha = "
                f"{fortran_real(float(troe['A']))}",
                f"    reactions({r_index})%troe%temperature_3 = "
                f"{fortran_real(float(troe['T3']))}",
                f"    reactions({r_index})%troe%temperature_1 = "
                f"{fortran_real(float(troe['T1']))}",
                f"    reactions({r_index})%troe%temperature_2 = "
                f"{fortran_real(float(troe.get('T2', 0.0)))}",
            ]
        lines.append(
            f"    reactions({r_index})%reversible = "
            f"{'.true.' if reaction.get('reversible', True) else '.false.'}"
        )

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
        f"    ok = size(species) == {prefix}_nspecies .and. &",
        f"      size(reactions) == {prefix}_nreactions",
        "    if (.not. ok) then",
        "      molar_production_rates = 0.0_dp",
        "      return",
        "    end if",
        "    call elementary_production_rates( &",
        "      species, reactions, temperature, density, mass_fractions, &",
        "      molar_production_rates, ok)",
        f"  end subroutine {kernel_name}",
        "",
        f"  subroutine {jacobian_name}( &",
        "      species, reactions, temperature, density, mass_fractions, &",
        "      jacobian, ok)",
        "    type(nasa7_species), intent(in) :: species(:)",
        "    type(elementary_reaction), intent(in) :: reactions(:)",
        "    real(dp), intent(in) :: temperature, density, mass_fractions(:)",
        "    real(dp), intent(out) :: jacobian(:, :)",
        "    logical, intent(out) :: ok",
        "",
        f"    ok = size(species) == {prefix}_nspecies .and. &",
        f"      size(reactions) == {prefix}_nreactions",
        "    if (.not. ok) then",
        "      jacobian = 0.0_dp",
        "      return",
        "    end if",
        "    call elementary_mass_fraction_jacobian( &",
        "      species, reactions, temperature, density, mass_fractions, &",
        "      jacobian, ok)",
        f"  end subroutine {jacobian_name}",
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
