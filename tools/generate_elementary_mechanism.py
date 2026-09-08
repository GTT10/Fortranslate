#!/usr/bin/env python3
"""Generate deterministic Fortran data for elementary/pressure-dependent kinetics."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any

REACTION_TYPES = {"elementary", "three-body", "falloff"}
CHEMISTRY_INTEGRATORS = {"explicit", "implicit"}
TRANSPORT_GEOMETRIES = {"atom": 0, "linear": 1, "nonlinear": 2}
MAX_SPECIES = 32
MAX_EQUATION_LENGTH = 128
MAX_FORTRAN_IDENTIFIER_LENGTH = 63
MAX_FORTRAN_LINE_LENGTH = 100
FORTRAN_IDENTIFIER = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
# Preserve chemical labels in strings; only generated index identifiers need
# Fortran spelling. Existing alphanumeric symbols remain byte-for-byte stable.
SPECIES_LABEL = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_()+.*-]*$")


def species_symbol_stem(name: str) -> str:
    if FORTRAN_IDENTIFIER.fullmatch(name):
        return name.lower()
    return "label_" + hashlib.sha256(name.lower().encode("ascii")).hexdigest()[:16]


FORTRAN_DEPENDENCY_MODULES = {
    "ieee_arithmetic",
    "precision_mod",
    "constants_mod",
    "nasa7_thermo_mod",
    "mixture_thermo_mod",
    "elementary_kinetics_mod",
    "constant_volume_reactor_mod",
    "selected_composition_mod",
    "selected_mechanism_runtime_mod",
    "simulation_config_selected_reactor_mod",
    "state_indices_mod",
    "gas_transport_mod",
    "mixture_transport_mod",
    "simulation_config_reactive_1d_mod",
    "slope_limiter_mod",
    "reconstruction_weno_mod",
    "reactive_1d_mod",
    "amr_hierarchy_1d_mod",
    "amr_multipatch_1d_mod",
    "amr_regrid_1d_mod",
    "amr_reactive_1d_mod",
    "amr_multilevel_reactive_1d_mod",
    "amr_multipatch_reactive_1d_mod",
    "amr_reactive_1d_application_mod",
    "mpi_f08",
    "mpi_domain_1d_mod",
    "mpi_reactive_transport_1d_mod",
    "mpi_reactive_1d_mod",
    "mpi_reactive_1d_application_mod",
    "mesh_mod",
    "mesh_3d_mod",
    "simulation_config_reactive_2d_mod",
    "simulation_config_reactive_3d_mod",
    "reactive_boundary_2d_mod",
    "reactive_transport_2d_mod",
    "reactive_2d_mod",
    "reactive_2d_application_mod",
    "eb_geometry_2d_mod",
    "simulation_config_reactive_eb_2d_mod",
    "reactive_eb_cfl_2d_mod",
    "eb_reactive_wall_flux_2d_mod",
    "eb_reactive_redistribution_2d_mod",
    "eb_reactive_reconstruction_2d_mod",
    "eb_reactive_hydro_2d_mod",
    "eb_reactive_transport_2d_mod",
    "reactive_eb_2d_driver_mod",
    "reactive_eb_2d_application_mod",
    "simulation_config_reactive_eb_amr_2d_mod",
    "amr_eb_hierarchy_2d_mod",
    "amr_eb_patch_tree_2d_mod",
    "amr_eb_multilevel_2d_mod",
    "amr_eb_flux_register_2d_mod",
    "amr_eb_reactive_2d_mod",
    "amr_eb_patch_tree_reactive_2d_mod",
    "amr_eb_multilevel_reactive_2d_mod",
    "amr_eb_regrid_2d_mod",
    "amr_eb_transport_2d_mod",
    "amr_eb_multilevel_transport_2d_mod",
    "amr_eb_multipatch_transport_2d_mod",
    "reactive_eb_amr_2d_driver_mod",
    "reactive_eb_amr_2d_application_mod",
    "reactive_directional_flux_3d_mod",
    "reactive_transport_3d_mod",
    "reactive_3d_mod",
    "reactive_entropy_wave_3d_problem_mod",
    "reactive_csv_io_3d_mod",
    "reactive_3d_application_mod",
    "eb_geometry_3d_mod",
    "simulation_config_reactive_eb_3d_mod",
    "reactive_eb_cfl_3d_mod",
    "eb_reactive_wall_flux_3d_mod",
    "eb_reactive_redistribution_3d_mod",
    "eb_reactive_hydro_3d_mod",
    "eb_reactive_transport_3d_mod",
    "reactive_eb_3d_checkpoint_mod",
    "reactive_eb_3d_driver_mod",
    "reactive_eb_3d_application_mod",
}
FORTRAN_IMPORTED_IDENTIFIERS = {
    "dp",
    "nasa7_species",
    "valid_nasa7_species",
    "elementary_reaction",
    "elementary_production_rates",
    "elementary_mass_fraction_jacobian",
    "reaction_kind_elementary",
    "reaction_kind_three_body",
    "reaction_kind_falloff",
}
# Names in the selected-mechanism probe's host scope and in the generated
# module's procedure scopes.  Every public name emitted by a bundle is
# imported into the probe, and each generated procedure has these dummy/local
# names.  Keep this list in sync with app/pelef_mechanism_probe.F90.in and
# generate() below.
FORTRAN_PROBE_RESERVED_IDENTIFIERS = {
    "pelef_mechanism_probe",
    "pelef0d_selected",
    "pelef_reactive_1d_selected",
    "pelef_amr_reactive_1d_selected",
    "pelef_mpi_reactive_1d_selected",
    "pelef_reactive_2d_selected",
    "pelef_reactive_eb_2d_selected",
    "pelef_reactive_eb_amr_2d_selected",
    "pelef_reactive_3d_selected",
    "pelef_reactive_eb_3d_selected",
    "pelef_selected_species_count",
    "pelef_selected_reaction_count",
    "pelef_selected_chemistry_integrator",
    "pelef_load_selected_reactions",
    "pelef_load_selected_thermo",
    "pelef_load_selected_transport",
    "pelef_selected_bundle_sha256",
    "ieee_arithmetic",
    "ieee_is_finite",
    "precision_mod",
    "nasa7_thermo_mod",
    "nasa7_species",
    "elementary_kinetics_mod",
    "elementary_reaction",
    "nspecies",
    "nreactions",
    "species",
    "reactions",
    "transport_names",
    "transport_geometries",
    "transport_values",
    "mass_fractions",
    "production_rates",
    "jacobian",
    "probe_temperature",
    "ok",
    "species_index",
    "names",
    "geometries",
    "values",
    "temperature",
    "density",
    "molar_production_rates",
    "allocated",
    "size",
    "shape",
    "all",
    "maxval",
    "minval",
    "trim",
    "any",
    "abs",
    "real",
    "integer",
    "logical",
    "character",
    "allocatable",
    "intent",
    "in",
    "out",
    "len",
    "intrinsic",
    "write",
}
FORTRAN_KEYWORDS = {
    "allocate", "block", "call", "case", "class", "contains", "cycle",
    "deallocate", "do", "else", "elseif", "end", "enddo", "endif",
    "entry", "error", "exit", "function", "if", "implicit", "import",
    "interface", "module", "none", "only", "parameter", "private",
    "procedure", "program", "public", "result", "return", "select",
    "stop", "subroutine", "then", "type", "use", "where", "while",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def fortran_real(value: float) -> str:
    return f"{value:.12e}_dp"


def finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label}: expected a real number")
    try:
        converted = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: expected a real number") from exc
    if not math.isfinite(converted):
        raise ValueError(f"{label}: nonfinite value")
    return converted


def validate_arrhenius(rate: dict[str, Any], equation: str) -> None:
    if finite_float(rate["A"], equation) < 0.0:
        raise ValueError(f"{equation}: negative Arrhenius A")
    for key in ("A", "b", "Ea"):
        finite_float(rate[key], equation)


def validate_species_data(data: dict[str, Any], species: list[str]) -> None:
    thermo = data.get("thermo")
    transport = data.get("transport")
    if thermo is None and transport is None:
        return
    if not isinstance(thermo, list) or not isinstance(transport, list):
        raise ValueError("thermo and transport records must be provided together")
    if len(thermo) != len(species) or len(transport) != len(species):
        raise ValueError("thermo and transport records must match species")
    for loader_key in ("thermo_loader_name", "transport_loader_name"):
        loader = data.get(loader_key)
        if (
            not isinstance(loader, str)
            or not FORTRAN_IDENTIFIER.fullmatch(loader)
            or len(loader) > MAX_FORTRAN_IDENTIFIER_LENGTH
            or loader.lower() in FORTRAN_KEYWORDS
        ):
            raise ValueError("generated species data requires valid loader names")

    for index, name in enumerate(species):
        thermo_record = thermo[index]
        if thermo_record.get("name") != name:
            raise ValueError(f"{name}: thermo order mismatch")
        composition = thermo_record.get("composition")
        if not isinstance(composition, dict) or not composition:
            raise ValueError(f"{name}: missing element composition")
        if any(
            not isinstance(element, str)
            or not element
            or finite_float(amount, name) <= 0.0
            for element, amount in composition.items()
        ):
            raise ValueError(f"{name}: invalid element composition")
        if len({element.lower() for element in composition}) != len(composition):
            raise ValueError(f"{name}: element names collide case-insensitively")
        if len(name) > 24:
            raise ValueError(f"{name}: species name exceeds Fortran storage")
        temperature_min = finite_float(thermo_record["temperature_min"], name)
        temperature_mid = finite_float(thermo_record["temperature_mid"], name)
        temperature_max = finite_float(thermo_record["temperature_max"], name)
        if not 0.0 < temperature_min < temperature_mid < temperature_max:
            raise ValueError(f"{name}: invalid NASA7 temperature ranges")
        if finite_float(thermo_record["molecular_weight"], name) <= 0.0:
            raise ValueError(f"{name}: invalid molecular weight")
        reference_pressure = finite_float(
            thermo_record.get("reference_pressure", 101325.0), name
        )
        if abs(reference_pressure - 101325.0) > 1.0e-9:
            raise ValueError(f"{name}: unsupported NASA7 reference pressure")
        for key in ("low_coefficients", "high_coefficients"):
            coefficients = thermo_record[key]
            if len(coefficients) != 7:
                raise ValueError(f"{name}: {key} must contain seven values")
            for value in coefficients:
                finite_float(value, name)

        transport_record = transport[index]
        if transport_record.get("name") != name:
            raise ValueError(f"{name}: transport order mismatch")
        if transport_record.get("geometry") not in TRANSPORT_GEOMETRIES:
            raise ValueError(f"{name}: unsupported transport geometry")
        if finite_float(transport_record["well_depth"], name) <= 0.0:
            raise ValueError(f"{name}: invalid transport well depth")
        if finite_float(transport_record["diameter"], name) <= 0.0:
            raise ValueError(f"{name}: invalid transport diameter")
        for key in ("dipole", "polarizability", "rotational_relaxation"):
            if finite_float(transport_record[key], name) < 0.0:
                raise ValueError(f"{name}: invalid transport {key}")


def validate_source(source: Any) -> None:
    if source is None:
        return
    if isinstance(source, str):
        if not source:
            raise ValueError("source provenance is empty")
        return
    if not isinstance(source, dict):
        raise ValueError("source provenance must be an object")
    for key in (
        "format", "file", "sha256", "cantera_version",
        "runtime_cantera_version", "runtime_cantera_git_commit", "phase",
    ):
        if not isinstance(source.get(key), str) or not source[key]:
            raise ValueError(f"source provenance is missing {key}")
    digest = source["sha256"]
    if len(digest) != 64 or any(
        character not in "0123456789abcdef" for character in digest
    ):
        raise ValueError("source provenance has invalid sha256")
    for key in (
        "phase", "cantera_version", "runtime_cantera_version",
        "runtime_cantera_git_commit",
    ):
        if any(character in source[key] for character in ('"', "\n", "\r")):
            raise ValueError(f"source provenance {key} is not Fortran-literal safe")
    for key in ("source_units", "target_units"):
        units = source.get(key)
        if not isinstance(units, dict) or any(
            not isinstance(name, str) or not isinstance(value, str)
            for name, value in units.items()
        ):
            raise ValueError(f"source provenance is missing valid {key}")
    expected_target_units = {
        "length": "m",
        "time": "s",
        "quantity": "kmol",
        "activation-energy": "J/kmol",
    }
    if source["target_units"] != expected_target_units:
        raise ValueError("source provenance has unsupported target units")


def validate(data: dict[str, Any]) -> None:
    if "schema_version" in data and data["schema_version"] != 1:
        raise ValueError("unsupported normalized mechanism schema")
    species = data["species"]
    if not species or len(species) != len(set(species)):
        raise ValueError("species must be non-empty and unique")
    if any(
        not isinstance(name, str)
        or not SPECIES_LABEL.fullmatch(name)
        or len(name) > 24
        for name in species
    ):
        raise ValueError("species names must be Fortran-safe identifiers or supported chemical labels")
    if len({name.lower() for name in species}) != len(species):
        raise ValueError("species names collide in case-insensitive Fortran")
    if len(species) > MAX_SPECIES:
        raise ValueError(f"mechanism exceeds {MAX_SPECIES} species")
    if not data.get("reactions"):
        raise ValueError("mechanism must contain reactions")
    identifiers = ["module_name", "loader_name", "kernel_name"]
    identifiers += ["jacobian_name", "symbol_prefix"]
    for required in identifiers:
        if required in {"jacobian_name", "symbol_prefix"} and required not in data:
            continue
        value = data.get(required)
        if value is None or value == "":
            raise ValueError(f"missing {required}")
        if (
            not isinstance(value, str)
            or not FORTRAN_IDENTIFIER.fullmatch(value)
            or len(value) > MAX_FORTRAN_IDENTIFIER_LENGTH
        ):
            raise ValueError(f"{required} is not a valid Fortran identifier")
        if required != "symbol_prefix" and value.lower() in FORTRAN_KEYWORDS:
            raise ValueError(f"{required} is a reserved Fortran keyword")

    module_name = data["module_name"]
    if module_name.lower() in FORTRAN_DEPENDENCY_MODULES:
        raise ValueError(
            f"module_name collides with generator dependency module {module_name}"
        )

    if "jacobian_name" not in data:
        jacobian_name = f"{data['kernel_name']}_jacobian"
        if (
            not FORTRAN_IDENTIFIER.fullmatch(jacobian_name)
            or len(jacobian_name) > MAX_FORTRAN_IDENTIFIER_LENGTH
            or jacobian_name.lower() in FORTRAN_KEYWORDS
        ):
            raise ValueError(
                "derived jacobian_name is not a valid Fortran identifier"
            )
    else:
        jacobian_name = data["jacobian_name"]

    validate_species_data(data, species)
    validate_source(data.get("source"))
    if "thermo" in data or "transport" in data:
        chemistry_integrator = data.get("chemistry_integrator")
        if chemistry_integrator not in CHEMISTRY_INTEGRATORS:
            raise ValueError(
                "complete mechanism bundle requires chemistry_integrator "
                "'explicit' or 'implicit'"
            )

    generated_identifiers = [
        data["module_name"],
        data["loader_name"],
        data["kernel_name"],
        jacobian_name,
    ]
    if "thermo" in data:
        generated_identifiers.extend((
            data["thermo_loader_name"], data["transport_loader_name"]
        ))
    if len({name.lower() for name in generated_identifiers}) != len(
        generated_identifiers
    ):
        raise ValueError("generated Fortran identifiers must be distinct")
    imported_identifiers = set(FORTRAN_IMPORTED_IDENTIFIERS)
    if "thermo" not in data:
        imported_identifiers.discard("valid_nasa7_species")
    for name in generated_identifiers:
        if name.lower() in imported_identifiers:
            raise ValueError(
                f"generated identifier {name} collides with imported Fortran identifier"
            )
    prefix = data.get("symbol_prefix", "h2o2")
    generated_symbols = [
        f"{prefix}_nspecies",
        f"{prefix}_nreactions",
        f"{prefix}_source_sha256",
        f"{prefix}_source_phase",
        f"{prefix}_source_cantera_version",
        f"{prefix}_runtime_cantera_version",
        f"{prefix}_runtime_cantera_git_commit",
        f"{prefix}_source_indices",
        f"{prefix}_duplicate_reactions",
    ]
    if "thermo" in data or "transport" in data:
        generated_symbols.append(f"{prefix}_chemistry_integrator")
    generated_symbols.extend(
        f"{prefix}_{species_symbol_stem(name)}_index" for name in species
    )
    if any(len(name) > MAX_FORTRAN_IDENTIFIER_LENGTH for name in generated_symbols):
        raise ValueError("generated Fortran symbol exceeds identifier length limit")
    all_generated_names = generated_identifiers + generated_symbols
    if len({name.lower() for name in all_generated_names}) != len(
        all_generated_names
    ):
        raise ValueError("generated Fortran names collide")
    probe_collisions = sorted(
        name
        for name in all_generated_names
        if name.lower() in FORTRAN_PROBE_RESERVED_IDENTIFIERS
    )
    if probe_collisions:
        raise ValueError(
            "generated Fortran names collide with probe namespace: "
            + ", ".join(probe_collisions)
        )

    has_reaction_provenance = any(
        "source_index" in reaction or "duplicate" in reaction
        for reaction in data["reactions"]
    )
    if has_reaction_provenance:
        for index, reaction in enumerate(data["reactions"], start=1):
            if reaction.get("source_index") != index:
                raise ValueError("reaction source indices must be ordered and one-based")
            if not isinstance(reaction.get("duplicate"), bool):
                raise ValueError("reaction duplicate flags must be logical")

    for reaction in data["reactions"]:
        equation = reaction.get("equation")
        if not isinstance(equation, str) or not equation:
            raise ValueError("reaction equation must be a non-empty string")
        if len(equation) > MAX_EQUATION_LENGTH:
            raise ValueError(
                f"{equation}: equation exceeds Fortran storage"
            )
        if any(character in equation for character in ('"', "\\", "\n", "\r")):
            raise ValueError(f"{equation}: equation is not Fortran-literal safe")
        if "reversible" in reaction and not isinstance(
            reaction["reversible"], bool
        ):
            raise ValueError(f"{equation}: reversible flag must be logical")
        reaction_type = reaction.get("type", "elementary")
        if reaction_type not in REACTION_TYPES:
            raise ValueError(f"{equation}: unsupported reaction type")
        for side in ("reactants", "products"):
            if not reaction[side]:
                raise ValueError(f"{equation}: empty {side}")
            for name, coefficient in reaction[side].items():
                if name not in species or finite_float(coefficient, equation) <= 0.0:
                    raise ValueError(f"{equation}: invalid {side}")
        if reaction_type == "falloff":
            validate_arrhenius(reaction["low_rate"], equation)
            validate_arrhenius(reaction["high_rate"], equation)
            troe = reaction.get("troe")
            if troe is not None:
                alpha = finite_float(troe["A"], equation)
                if not 0.0 < alpha < 1.0:
                    raise ValueError(f"{equation}: invalid Troe A")
                if (finite_float(troe["T3"], equation) <= 0.0 or
                        finite_float(troe["T1"], equation) <= 0.0):
                    raise ValueError(f"{equation}: invalid Troe temperatures")
                if finite_float(troe.get("T2", 0.0), equation) < 0.0:
                    raise ValueError(f"{equation}: invalid Troe T2")
        else:
            validate_arrhenius(reaction["arrhenius"], equation)
        if reaction_type in {"three-body", "falloff"}:
            default_efficiency = finite_float(
                reaction.get("default_efficiency", 1.0), equation
            )
            if default_efficiency < 0.0:
                raise ValueError(f"{equation}: negative default efficiency")
            for name, efficiency in reaction.get("efficiencies", {}).items():
                if name not in species or finite_float(efficiency, equation) < 0.0:
                    raise ValueError(f"{equation}: invalid third-body efficiency")

    if "thermo" in data:
        thermo_by_name = {
            record["name"]: record for record in data["thermo"]
        }
        elements = sorted({
            element
            for record in data["thermo"]
            for element in record["composition"]
        })
        for reaction in data["reactions"]:
            equation = reaction["equation"]
            reactant_mass = sum(
                float(coefficient)
                * float(thermo_by_name[name]["molecular_weight"])
                for name, coefficient in reaction["reactants"].items()
            )
            product_mass = sum(
                float(coefficient)
                * float(thermo_by_name[name]["molecular_weight"])
                for name, coefficient in reaction["products"].items()
            )
            mass_scale = max(1.0, abs(reactant_mass), abs(product_mass))
            if abs(product_mass - reactant_mass) > 2.0e-12 * mass_scale:
                raise ValueError(f"{equation}: reaction mass is not balanced")
            for element in elements:
                reactant_atoms = sum(
                    float(coefficient)
                    * float(
                        thermo_by_name[name]["composition"].get(element, 0.0)
                    )
                    for name, coefficient in reaction["reactants"].items()
                )
                product_atoms = sum(
                    float(coefficient)
                    * float(
                        thermo_by_name[name]["composition"].get(element, 0.0)
                    )
                    for name, coefficient in reaction["products"].items()
                )
                atom_scale = max(1.0, abs(reactant_atoms), abs(product_atoms))
                if abs(product_atoms - reactant_atoms) > 2.0e-12 * atom_scale:
                    raise ValueError(
                        f"{equation}: element {element} is not balanced"
                    )


def emit_rate(lines: list[str], target: str, rate: dict[str, Any]) -> None:
    lines += [
        f"    {target}%pre_exponential = {fortran_real(float(rate['A']))}",
        f"    {target}%temperature_exponent = {fortran_real(float(rate['b']))}",
        f"    {target}%activation_energy = {fortran_real(float(rate['Ea']))}",
    ]


def emit_fortran_string_assignment(
    lines: list[str], target: str, value: str
) -> None:
    """Emit a bounded-width Fortran character assignment."""

    assignment = f'    {target} = "{value}"'
    if len(assignment) <= MAX_FORTRAN_LINE_LENGTH:
        lines.append(assignment)
        return
    lines.append(f"    {target} = &")
    chunk_size = MAX_FORTRAN_LINE_LENGTH - len('      "" // &')
    chunks = [
        value[offset:offset + chunk_size]
        for offset in range(0, len(value), chunk_size)
    ]
    for index, chunk in enumerate(chunks):
        suffix = " // &" if index + 1 < len(chunks) else ""
        lines.append(f'      "{chunk}"{suffix}')


def emit_fortran_parameter_string(
    lines: list[str], name: str, value: str
) -> None:
    """Emit a bounded-width public character parameter declaration."""

    declaration = f"  character(len=*), parameter, public :: {name}"
    assignment = f"{declaration} ="
    if len(assignment) + 2 <= MAX_FORTRAN_LINE_LENGTH:
        lines.append(f"{assignment} &")
    else:
        lines.append("  character(len=*), parameter, public :: &")
        lines.append(f"    {name} = &")

    chunk_size = MAX_FORTRAN_LINE_LENGTH - len('    "" // &')
    chunks = [
        value[offset:offset + chunk_size]
        for offset in range(0, len(value), chunk_size)
    ]
    for index, chunk in enumerate(chunks):
        suffix = " // &" if index + 1 < len(chunks) else ""
        lines.append(f'    "{chunk}"{suffix}')


def emit_species_data_loaders(lines: list[str], data: dict[str, Any]) -> None:
    thermo: list[dict[str, Any]] = data["thermo"]
    transport: list[dict[str, Any]] = data["transport"]
    thermo_loader = data["thermo_loader_name"]
    transport_loader = data["transport_loader_name"]
    nspecies = len(thermo)

    lines += [
        "",
        f"  subroutine {thermo_loader}(species, ok)",
        "    type(nasa7_species), allocatable, intent(out) :: species(:)",
        "    logical, intent(out) :: ok",
        "    integer :: species_index",
        "",
        f"    allocate(species({nspecies}))",
    ]
    for index, record in enumerate(thermo, start=1):
        lines += [
            "",
            f"    species({index})%name = \"{record['name']}\"",
            f"    species({index})%molecular_weight = "
            f"{fortran_real(float(record['molecular_weight']))}",
            f"    species({index})%temperature_min = "
            f"{fortran_real(float(record['temperature_min']))}",
            f"    species({index})%temperature_mid = "
            f"{fortran_real(float(record['temperature_mid']))}",
            f"    species({index})%temperature_max = "
            f"{fortran_real(float(record['temperature_max']))}",
            f"    species({index})%low_coefficients = [ &",
        ]
        low = record["low_coefficients"]
        high = record["high_coefficients"]
        lines += [
            "      " + ", ".join(fortran_real(float(value)) for value in low[:3]) + ", &",
            "      " + ", ".join(fortran_real(float(value)) for value in low[3:6]) + ", &",
            f"      {fortran_real(float(low[6]))} ]",
            f"    species({index})%high_coefficients = [ &",
            "      " + ", ".join(fortran_real(float(value)) for value in high[:3]) + ", &",
            "      " + ", ".join(fortran_real(float(value)) for value in high[3:6]) + ", &",
            f"      {fortran_real(float(high[6]))} ]",
        ]
    lines += [
        "",
        "    ok = .true.",
        "    do species_index = 1, size(species)",
        "      if (.not. valid_nasa7_species(species(species_index))) then",
        "        ok = .false.",
        "        return",
        "      end if",
        "    end do",
        f"  end subroutine {thermo_loader}",
        "",
        f"  subroutine {transport_loader}( &",
        "      names, geometries, values, ok)",
        "    character(len=*), intent(out) :: names(:)",
        "    integer, intent(out) :: geometries(:)",
        "    real(dp), intent(out) :: values(:, :)",
        "    logical, intent(out) :: ok",
        "",
        f"    ok = size(names) == {nspecies} .and. &",
        f"      size(geometries) == {nspecies} .and. &",
        f"      all(shape(values) == [5, {nspecies}])",
        "    if (.not. ok) return",
        "    names = \"\"",
        "    geometries = 0",
        "    values = 0.0_dp",
    ]
    for index, record in enumerate(transport, start=1):
        values = [
            record["well_depth"],
            record["diameter"],
            record["dipole"],
            record["polarizability"],
            record["rotational_relaxation"],
        ]
        lines += [
            f"    names({index}) = \"{record['name']}\"",
            f"    geometries({index}) = "
            f"{TRANSPORT_GEOMETRIES[record['geometry']]}",
            f"    values(:, {index}) = [ &",
            "      " + ", ".join(fortran_real(float(value)) for value in values[:3]) + ", &",
            "      " + ", ".join(fortran_real(float(value)) for value in values[3:]) + " ]",
        ]
    lines += [
        f"  end subroutine {transport_loader}",
    ]


def emit_parameter_array(
    lines: list[str], declaration: str, values: list[str], chunk_size: int
) -> None:
    lines.append(f"  {declaration} = [ &")
    for offset in range(0, len(values), chunk_size):
        chunk = values[offset:offset + chunk_size]
        suffix = ", &" if offset + chunk_size < len(values) else " ]"
        lines.append("    " + ", ".join(chunk) + suffix)


def generate(data: dict[str, Any]) -> str:
    species: list[str] = data["species"]
    reactions: list[dict[str, Any]] = data["reactions"]
    module_name = data["module_name"]
    loader_name = data["loader_name"]
    kernel_name = data["kernel_name"]
    jacobian_name = data.get("jacobian_name", f"{kernel_name}_jacobian")
    prefix = data.get("symbol_prefix", "h2o2")
    has_species_data = "thermo" in data
    nasa_import = "  use nasa7_thermo_mod, only: nasa7_species"
    if has_species_data:
        nasa_import += ", valid_nasa7_species"

    lines: list[str] = [
        f"module {module_name}",
        "  use precision_mod, only: dp",
        nasa_import,
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
    if "thermo" in data or "transport" in data:
        emit_fortran_parameter_string(
            lines,
            f"{prefix}_chemistry_integrator",
            data["chemistry_integrator"],
        )
    source = data.get("source")
    if isinstance(source, dict):
        for name, value in (
            (f"{prefix}_source_sha256", source["sha256"]),
            (f"{prefix}_source_phase", source["phase"]),
            (f"{prefix}_source_cantera_version", source["cantera_version"]),
            (
                f"{prefix}_runtime_cantera_version",
                source["runtime_cantera_version"],
            ),
            (
                f"{prefix}_runtime_cantera_git_commit",
                source["runtime_cantera_git_commit"],
            ),
        ):
            emit_fortran_parameter_string(lines, name, value)
    if all("source_index" in reaction for reaction in reactions):
        emit_parameter_array(
            lines,
            f"integer, parameter, public :: {prefix}_source_indices({len(reactions)})",
            [str(reaction["source_index"]) for reaction in reactions],
            10,
        )
        emit_parameter_array(
            lines,
            f"logical, parameter, public :: {prefix}_duplicate_reactions({len(reactions)})",
            [".true." if reaction["duplicate"] else ".false." for reaction in reactions],
            7,
        )
    for index, name in enumerate(species, start=1):
        lines.append(
            f"  integer, parameter, public :: {prefix}_{species_symbol_stem(name)}_index = {index}"
        )
    lines += [
        "",
        f"  public :: {loader_name}",
        f"  public :: {kernel_name}",
        f"  public :: {jacobian_name}",
    ]
    if has_species_data:
        lines += [
            f"  public :: {data['thermo_loader_name']}",
            f"  public :: {data['transport_loader_name']}",
        ]
    lines += [
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
        lines.append("")
        emit_fortran_string_assignment(
            lines, f"reactions({r_index})%equation", reaction["equation"]
        )
        lines += [
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
    ]
    if has_species_data:
        emit_species_data_loaders(lines, data)
    lines += [
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


def write_generated(output_path: Path, output: str) -> None:
    """Atomically publish generated source beside its final destination."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(output)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, output_path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    args = parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    validate(data)
    output = generate(data)
    write_generated(args.output, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
