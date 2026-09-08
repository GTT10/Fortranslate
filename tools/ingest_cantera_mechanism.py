#!/usr/bin/env python3
"""Convert a supported Cantera YAML gas phase to the normalized JSON bundle.

The normalized bundle is intentionally a superset of the input accepted by
``generate_elementary_mechanism.py``.  The existing top-level ``species`` and
``reactions`` fields retain their generator schema, while ``thermo`` and
``transport`` carry the records needed by the runtime loaders.  Cantera is
asked for ``input_data`` instead of reparsing YAML scalars; this gives one
well-defined unit conversion path (the generated kinetics uses J/kmol and
Cantera's normalized Arrhenius values).

The accepted mechanism subset is deliberately small:

* NASA7 species with gas Lennard-Jones transport data;
* elementary Arrhenius reactions;
* three-body Arrhenius reactions; and
* Lindemann falloff and falloff reactions with Troe broadening.

PLOG, Chebyshev, SRI, NASA9, non-gas transport, and missing transport are
rejected before a bundle is written.  Reaction list order is never sorted, so
Cantera duplicate reaction ordering is preserved exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence


MAX_SPECIES = 32
MAX_EQUATION_LENGTH = 128
REFERENCE_PRESSURE = 101325.0
SUPPORTED_CANTERA_VERSION = "3.2.0"
SUPPORTED_CANTERA_GIT_COMMIT = "4a8358e"
_SPECIES_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_SPECIES_LABEL = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_()+.*-]*$")
_TRANSPORT_GEOMETRIES = {"atom", "linear", "nonlinear"}
_TRANSPORT_KEYS = {
    "model",
    "geometry",
    "well-depth",
    "diameter",
    "dipole",
    "polarizability",
    "rotational-relaxation",
    "note",
}
_TARGET_UNITS = {
    "length": "m",
    "time": "s",
    "quantity": "kmol",
    "activation-energy": "J/kmol",
}


def _as_float(value: Any, label: str) -> float:
    """Return a finite real value, rejecting booleans and non-numbers."""

    if isinstance(value, bool):
        raise ValueError(f"{label}: expected a real number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: expected a real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label}: non-finite value")
    return result


def _positive(value: Any, label: str) -> float:
    result = _as_float(value, label)
    if result <= 0.0:
        raise ValueError(f"{label}: expected a positive value")
    return result


def _nonnegative(value: Any, label: str) -> float:
    result = _as_float(value, label)
    if result < 0.0:
        raise ValueError(f"{label}: expected a non-negative value")
    return result


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label}: expected a mapping")
    return value


def _key_case_insensitive(mapping: Mapping[str, Any], key: str) -> Any:
    """Look up a YAML/input_data key without depending on Cantera casing."""

    lowered = key.lower()
    for candidate, value in mapping.items():
        if str(candidate).lower() == lowered:
            return value
    return None


def _has_key_case_insensitive(mapping: Mapping[str, Any], key: str) -> bool:
    lowered = key.lower()
    return any(str(candidate).lower() == lowered for candidate in mapping)


def _load_cantera() -> Any:
    try:
        import cantera as ct
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError(
            "Cantera Python package is required for YAML ingestion"
        ) from exc
    return ct


def _validate_cantera_runtime(cantera: Any) -> tuple[str, str]:
    version = str(getattr(cantera, "__version__", ""))
    if version != SUPPORTED_CANTERA_VERSION:
        raise ValueError(
            f"unsupported Cantera runtime {version or '<missing>'}; "
            f"expected {SUPPORTED_CANTERA_VERSION}"
        )
    commit = str(getattr(cantera, "__git_commit__", ""))
    if commit != SUPPORTED_CANTERA_GIT_COMMIT:
        raise ValueError(
            f"unsupported Cantera git commit {commit or '<missing>'}; "
            f"expected {SUPPORTED_CANTERA_GIT_COMMIT}"
        )
    return version, commit


def _validate_phase_models(cantera_phase: Any) -> None:
    name = str(getattr(cantera_phase, "name", "<unnamed>"))
    phase_thermo = str(getattr(cantera_phase, "thermo_model", "")).lower()
    if phase_thermo != "ideal-gas":
        raise ValueError(
            f"phase {name!r}: unsupported thermo model "
            f"{phase_thermo or '<missing>'}; only ideal-gas is supported"
        )
    phase_kinetics = str(getattr(cantera_phase, "kinetics_model", "")).lower()
    if phase_kinetics != "bulk":
        raise ValueError(
            f"phase {name!r}: unsupported kinetics model "
            f"{phase_kinetics or '<missing>'}; only bulk is supported"
        )
    phase_transport = str(getattr(cantera_phase, "transport_model", "")).lower()
    if phase_transport != "mixture-averaged":
        raise ValueError(
            f"phase {name!r}: unsupported transport model "
            f"{phase_transport or '<missing>'}; only mixture-averaged is supported"
        )


def load_cantera_phase(input_path: Path, phase: str, cantera: Any | None = None) -> Any:
    """Load one explicitly selected Cantera phase."""

    if not phase:
        raise ValueError("phase is required; refusing an ambiguous multi-phase YAML")
    ct = cantera or _load_cantera()
    try:
        try:
            # ``name`` is the supported keyword in Cantera 3.x.  Some older
            # releases exposed the same selector as ``phaseid``.
            return ct.Solution(str(input_path), name=phase)
        except TypeError:
            return ct.Solution(str(input_path), phaseid=phase)
    except Exception as exc:
        selected = f" phase {phase!r}" if phase else ""
        raise ValueError(
            f"cannot load Cantera YAML {input_path}{selected}: {exc}"
        ) from exc


def _default_prefix(input_path: Path) -> str:
    stem = re.sub(r"[^A-Za-z0-9_]", "_", input_path.stem.lower())
    stem = stem.strip("_") or "mechanism"
    if stem[0].isdigit():
        stem = f"mechanism_{stem}"
    return stem


def _identifier(value: str, label: str) -> str:
    if not _SPECIES_IDENTIFIER.fullmatch(value):
        raise ValueError(f"{label} {value!r} is not a valid identifier")
    return value


def _validate_species_names(names: Sequence[str]) -> list[str]:
    if not names:
        raise ValueError("phase contains no species")
    if len(names) > MAX_SPECIES:
        raise ValueError(
            f"phase contains {len(names)} species; maximum supported is "
            f"{MAX_SPECIES}"
        )
    result = [str(name) for name in names]
    if len(result) != len(set(result)):
        raise ValueError("phase species names must be unique")
    if len({name.lower() for name in result}) != len(result):
        raise ValueError("phase species names collide in case-insensitive Fortran")
    for name in result:
        if len(name) > 24:
            raise ValueError(
                f"species {name!r} exceeds the 24-character runtime limit"
            )
        if not _SPECIES_LABEL.fullmatch(name):
            raise ValueError(
                f"species {name!r} is not a supported generator-safe chemical label"
            )
    return result


def _yaml_metadata(source_bytes: bytes) -> tuple[str | None, dict[str, str]]:
    """Read portable top-level provenance scalars without a YAML dependency.

    Cantera exposes normalized phase/species data but not the document header
    (notably ``cantera-version`` and original ``units``).  These fields have a
    deliberately simple top-level YAML representation in mechanism files, so
    extracting only those scalars keeps this tool usable in the small build
    Python environment without adding PyYAML as a dependency.
    """

    try:
        text = source_bytes.decode("utf-8")
    except UnicodeDecodeError:
        return None, {}
    version_match = re.search(
        r"(?m)^cantera-version:\s*([^#\n]+?)\s*$", text
    )
    version = version_match.group(1).strip().strip("'\"") if version_match else None
    units_match = re.search(r"(?m)^units:\s*\{([^}]*)\}", text)
    units: dict[str, str] = {}
    if units_match:
        for item in units_match.group(1).split(","):
            if ":" not in item:
                continue
            key, value = item.split(":", 1)
            units[key.strip().strip("'\"")] = value.strip().strip("'\"")
    else:
        block_match = re.search(
            r"(?ms)^units:\s*(?:#[^\n]*)?\n((?:[ \t]+[^\n]*\n?)*)",
            text,
        )
        if block_match:
            for line in block_match.group(1).splitlines():
                content = line.split("#", 1)[0].strip()
                if ":" not in content:
                    continue
                key, value = content.split(":", 1)
                units[key.strip().strip("'\"")] = value.strip().strip("'\"")
    return version, units


def _thermo_record(species: Any, index: int) -> dict[str, Any]:
    name = str(species.name)
    raw_species = _mapping(getattr(species, "input_data", None), f"species {name}")
    thermo = _mapping(raw_species.get("thermo"), f"species {name} thermo")
    charge = _as_float(
        getattr(species, "charge", raw_species.get("charge", 0.0)),
        f"species {name} charge",
    )
    if charge != 0.0:
        raise ValueError(f"species {name}: charged species are unsupported")
    model = str(thermo.get("model", "")).upper()
    if model != "NASA7":
        raise ValueError(
            f"species {name}: unsupported thermo model {model or '<missing>'}; "
            "only NASA7 is supported"
        )
    reference_pressure = _positive(
        getattr(getattr(species, "thermo", None), "reference_pressure", None),
        f"species {name} NASA7 reference pressure",
    )
    if not math.isclose(
        reference_pressure, REFERENCE_PRESSURE, rel_tol=0.0, abs_tol=1.0e-9
    ):
        raise ValueError(
            f"species {name}: unsupported NASA7 reference pressure "
            f"{reference_pressure:g} Pa"
        )

    ranges = thermo.get("temperature-ranges")
    if not isinstance(ranges, Sequence) or isinstance(ranges, (str, bytes)):
        raise ValueError(f"species {name}: NASA7 temperature-ranges is invalid")
    if len(ranges) not in (2, 3):
        raise ValueError(
            f"species {name}: NASA7 temperature-ranges must contain two or three values"
        )
    temperatures = [
        _positive(value, f"species {name} temperature-ranges[{item}]")
        for item, value in enumerate(ranges)
    ]
    single_interval = len(temperatures) == 2
    if single_interval:
        low, high = temperatures
        temperatures = [low, low + 0.5 * (high - low), high]
    data = thermo.get("data")
    if not isinstance(data, Sequence) or isinstance(data, (str, bytes)):
        raise ValueError(f"species {name}: NASA7 data is invalid")
    if single_interval:
        if len(data) != 1:
            raise ValueError(f"species {name}: single-interval NASA7 requires exactly one row")
        # Duplicate the same polynomial on two subintervals, without changing
        # its function, coefficients, reference pressure or validity interval.
        data = [data[0], data[0]]
    elif len(data) != 2:
        raise ValueError(f"species {name}: NASA7 data must contain low and high rows")
    coefficients: list[list[float]] = []
    for row_index, row in enumerate(data):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            raise ValueError(f"species {name}: NASA7 row {row_index} is invalid")
        if len(row) != 7:
            raise ValueError(
                f"species {name}: NASA7 row {row_index} must contain seven values"
            )
        coefficients.append(
            [
                _as_float(value, f"species {name} NASA7[{row_index}][{coefficient}]")
                for coefficient, value in enumerate(row)
            ]
        )

    if not temperatures[0] < temperatures[1] < temperatures[2]:
        # Cantera 3.2 exports a one-interval NASA7 as [Tmin,Tmax,Tmax]
        # with two identical rows. Preserve the function with an interior split;
        # never "repair" unequal polynomials or a reversed/empty range.
        low, mid, high = temperatures
        if low < high and mid in (low, high) and coefficients[0] == coefficients[1]:
            temperatures[1] = low + 0.5 * (high - low)
        if not temperatures[0] < temperatures[1] < temperatures[2]:
            raise ValueError(f"species {name}: NASA7 temperature ranges are not increasing")

    composition_data = raw_species.get("composition")
    composition_mapping = _mapping(
        composition_data, f"species {name} composition"
    )
    if not composition_mapping:
        raise ValueError(f"species {name}: missing element composition")
    composition = {
        str(element): _positive(
            composition_mapping[element], f"species {name} composition {element}"
        )
        for element in sorted(composition_mapping, key=str)
    }

    return {
        "name": name,
        "composition": composition,
        "molecular_weight": _positive(
            getattr(species, "molecular_weight", None),
            f"species {name} molecular_weight",
        ),
        "reference_pressure": reference_pressure,
        "temperature_min": temperatures[0],
        "temperature_mid": temperatures[1],
        "temperature_max": temperatures[2],
        "low_coefficients": coefficients[0],
        "high_coefficients": coefficients[1],
    }


def _transport_record(species: Any, index: int) -> dict[str, Any]:
    name = str(species.name)
    raw_species = _mapping(getattr(species, "input_data", None), f"species {name}")
    raw_transport = raw_species.get("transport")
    if raw_transport is None:
        raise ValueError(f"species {name}: missing transport data")
    transport = _mapping(raw_transport, f"species {name} transport")
    unknown_keys = [
        str(key) for key in transport if str(key).lower() not in _TRANSPORT_KEYS
    ]
    if unknown_keys:
        raise ValueError(
            f"species {name}: unsupported transport attributes "
            f"{', '.join(sorted(unknown_keys))}"
        )

    model = str(transport.get("model", "")).lower()
    if model != "gas":
        raise ValueError(
            f"species {name}: unsupported transport model {model or '<missing>'}; "
            "only gas Lennard-Jones transport is supported"
        )
    geometry = str(transport.get("geometry", "")).lower()
    if geometry not in _TRANSPORT_GEOMETRIES:
        raise ValueError(
            f"species {name}: unsupported transport geometry "
            f"{geometry or '<missing>'}"
        )
    if "well-depth" not in transport or "diameter" not in transport:
        raise ValueError(
            f"species {name}: transport requires well-depth and diameter"
        )

    return {
        "name": name,
        "geometry": geometry,
        "well_depth": _positive(
            transport["well-depth"], f"species {name} transport well-depth"
        ),
        "diameter": _positive(
            transport["diameter"], f"species {name} transport diameter"
        ),
        "dipole": _nonnegative(
            transport.get("dipole", 0.0), f"species {name} transport dipole"
        ),
        "polarizability": _nonnegative(
            transport.get("polarizability", 0.0),
            f"species {name} transport polarizability",
        ),
        "rotational_relaxation": _nonnegative(
            transport.get("rotational-relaxation", 0.0),
            f"species {name} transport rotational-relaxation",
        ),
    }


def _arrhenius_record(rate: Any, label: str) -> dict[str, float]:
    data = _mapping(rate, label)
    missing = [key for key in ("A", "b", "Ea") if key not in data]
    if missing:
        raise ValueError(f"{label}: missing {', '.join(missing)}")
    return {
        "A": _nonnegative(data["A"], f"{label} A"),
        "b": _as_float(data["b"], f"{label} b"),
        "Ea": _as_float(data["Ea"], f"{label} Ea"),
    }


def _reaction_kind(reaction: Any, raw: Mapping[str, Any], equation: str) -> str:
    reaction_type = str(raw.get("type", "")).lower()
    reaction_class = type(getattr(reaction, "rate", None)).__name__.lower()
    reaction_kind = str(getattr(reaction, "reaction_type", "")).lower()
    combined = " ".join((reaction_type, reaction_class, reaction_kind))
    rate = getattr(reaction, "rate", None)
    rate_subtype = str(getattr(rate, "sub_type", "")).lower()

    api_orders = getattr(reaction, "orders", {})
    if api_orders or bool(getattr(reaction, "allow_negative_orders", False)) or bool(
        getattr(reaction, "allow_nonreactant_orders", False)
    ):
        raise ValueError(f"{equation}: custom reaction orders are unsupported")
    if bool(getattr(rate, "chemically_activated", False)):
        raise ValueError(f"{equation}: chemically-activated rates are unsupported")
    if rate_subtype in {"sri", "tsang"}:
        raise ValueError(f"{equation}: {rate_subtype} falloff is unsupported")

    if "plog" in combined or "pressure-dependent-arrhenius" in combined:
        raise ValueError(f"{equation}: PLOG reactions are unsupported")
    if "chebyshev" in combined:
        raise ValueError(f"{equation}: Chebyshev reactions are unsupported")
    if "sri" in combined or _has_key_case_insensitive(raw, "SRI"):
        raise ValueError(f"{equation}: SRI falloff is unsupported")
    if "tsang" in combined or _has_key_case_insensitive(raw, "Tsang"):
        raise ValueError(f"{equation}: Tsang falloff is unsupported")
    if "chemically-activated" in combined:
        raise ValueError(f"{equation}: chemically-activated rates are unsupported")
    if any(_has_key_case_insensitive(raw, key) for key in (
        "orders", "negative-orders", "nonreactant-orders"
    )):
        raise ValueError(f"{equation}: custom reaction orders are unsupported")
    if any(_has_key_case_insensitive(raw, key) for key in (
        "reverse-rate-constant", "reverse_rate_constant"
    )):
        raise ValueError(f"{equation}: explicit reverse rates are unsupported")

    if reaction_type in {"three-body", "three-body-arrhenius"}:
        return "three-body"
    if "three-body" in combined and "falloff" not in combined:
        return "three-body"
    if reaction_type in {"falloff", "falloff-troe", "falloff-lindemann"}:
        return "falloff"
    if "falloff" in combined:
        return "falloff"
    if reaction_type in {"", "elementary", "arrhenius"}:
        return "elementary"
    if "arrhenius" in combined and getattr(reaction, "third_body", None) is None:
        return "elementary"
    raise ValueError(
        f"{equation}: unsupported reaction/rate type "
        f"{raw.get('type', getattr(reaction, 'reaction_type', '<unknown>'))!r}"
    )


def _ordered_stoich(
    values: Any, species_names: Sequence[str], equation: str, side: str
) -> dict[str, float]:
    mapping = _mapping(values, f"{equation} {side}")
    unknown = [name for name in mapping if name not in species_names]
    if unknown:
        raise ValueError(
            f"{equation}: {side} contains species not present in phase: "
            f"{', '.join(map(str, unknown))}"
        )
    result: dict[str, float] = {}
    for name in species_names:
        if name not in mapping:
            continue
        result[name] = _positive(mapping[name], f"{equation} {side} {name}")
    if not result:
        raise ValueError(f"{equation}: {side} is empty")
    return result


def _third_body_data(
    reaction: Any, raw: Mapping[str, Any], species_names: Sequence[str], equation: str
) -> tuple[float, dict[str, float]]:
    third_body = getattr(reaction, "third_body", None)
    if third_body is not None:
        default = getattr(third_body, "default_efficiency", 1.0)
        efficiencies = getattr(third_body, "efficiencies", {})
    else:
        default = raw.get("default-efficiency", raw.get("default_efficiency", 1.0))
        efficiencies = raw.get("efficiencies", {})
    default_value = _nonnegative(default, f"{equation} default third-body efficiency")
    efficiency_mapping = _mapping(
        efficiencies, f"{equation} third-body efficiencies"
    )
    unknown = [name for name in efficiency_mapping if name not in species_names]
    if unknown:
        raise ValueError(
            f"{equation}: third-body efficiency has unknown species "
            f"{', '.join(map(str, unknown))}"
        )
    ordered = {
        name: _nonnegative(
            efficiency_mapping[name], f"{equation} efficiency {name}"
        )
        for name in species_names
        if name in efficiency_mapping
    }
    return default_value, ordered


def _troe_record(raw: Mapping[str, Any], equation: str) -> dict[str, float] | None:
    troe = _key_case_insensitive(raw, "Troe")
    if troe is None:
        return None
    values = _mapping(troe, f"{equation} Troe")
    required = [key for key in ("A", "T3", "T1") if key not in values]
    if required:
        raise ValueError(f"{equation} Troe: missing {', '.join(required)}")
    result = {
        "A": _as_float(values["A"], f"{equation} Troe A"),
        "T3": _positive(values["T3"], f"{equation} Troe T3"),
        "T1": _positive(values["T1"], f"{equation} Troe T1"),
        "T2": _nonnegative(values.get("T2", 0.0), f"{equation} Troe T2"),
    }
    if not 0.0 < result["A"] < 1.0:
        raise ValueError(f"{equation} Troe A: expected a value between zero and one")
    return result


def _reaction_record(
    reaction: Any,
    species_names: Sequence[str],
    index: int,
    species_compositions: Mapping[str, Mapping[str, float]] | None = None,
) -> dict[str, Any]:
    raw = _mapping(
        getattr(reaction, "input_data", None), f"reaction {index + 1} input_data"
    )
    equation = str(getattr(reaction, "equation", raw.get("equation", "")))
    if not equation:
        raise ValueError(f"reaction {index + 1}: missing equation")
    if len(equation) > MAX_EQUATION_LENGTH:
        raise ValueError(
            f"{equation}: equation exceeds the {MAX_EQUATION_LENGTH}-character "
            "Fortran storage limit"
        )
    # The existing Fortran generator emits this value inside a quoted literal.
    # Rejecting quote characters gives a useful ingestion error instead of a
    # later generated-source syntax error.
    if any(character in equation for character in ('"', "\\", "\n", "\r")):
        raise ValueError(f"{equation}: equation contains an unsupported quote or slash")

    kind = _reaction_kind(reaction, raw, equation)
    record: dict[str, Any] = {
        "source_index": index + 1,
        "equation": equation,
        "type": kind,
        "reactants": _ordered_stoich(
            getattr(reaction, "reactants", None), species_names, equation, "reactants"
        ),
        "products": _ordered_stoich(
            getattr(reaction, "products", None), species_names, equation, "products"
        ),
        "reversible": bool(getattr(reaction, "reversible", True)),
        "duplicate": bool(
            getattr(reaction, "duplicate", raw.get("duplicate", False))
        ),
    }
    if species_compositions is not None:
        _validate_element_balance(record, species_compositions)

    if kind == "falloff":
        record["low_rate"] = _arrhenius_record(
            raw.get("low-P-rate-constant"), f"{equation} low-pressure rate"
        )
        record["high_rate"] = _arrhenius_record(
            raw.get("high-P-rate-constant"), f"{equation} high-pressure rate"
        )
        troe = _troe_record(raw, equation)
        if troe is None:
            # Require Cantera's explicit rate classification. Missing Troe data
            # must not silently downgrade a malformed/unknown falloff model.
            subtype = str(getattr(getattr(reaction, "rate", None), "sub_type", "")).lower()
            rate_kind = str(getattr(reaction, "reaction_type", "")).lower()
            if subtype != "lindemann" or rate_kind != "falloff-lindemann":
                raise ValueError(f"{equation}: missing Troe data for non-Lindemann falloff")
            # An absent troe field selects the existing runtime's F=1 path.
        else:
            record["troe"] = troe
    else:
        record["arrhenius"] = _arrhenius_record(
            raw.get("rate-constant"), f"{equation} Arrhenius rate"
        )

    if kind in {"three-body", "falloff"}:
        default_efficiency, efficiencies = _third_body_data(
            reaction, raw, species_names, equation
        )
        record["default_efficiency"] = default_efficiency
        record["efficiencies"] = efficiencies
    return record


def _validate_element_balance(
    reaction: Mapping[str, Any],
    species_compositions: Mapping[str, Mapping[str, float]],
) -> None:
    equation = str(reaction["equation"])
    elements = set()
    for side in ("reactants", "products"):
        for name in reaction[side]:
            if name not in species_compositions:
                raise ValueError(
                    f"{equation}: missing composition for species {name}"
                )
            elements.update(species_compositions[name])
    for element in sorted(elements):
        reactant_total = sum(
            float(coefficient) * float(species_compositions[name].get(element, 0.0))
            for name, coefficient in reaction["reactants"].items()
        )
        product_total = sum(
            float(coefficient) * float(species_compositions[name].get(element, 0.0))
            for name, coefficient in reaction["products"].items()
        )
        if not math.isclose(
            reactant_total, product_total, rel_tol=1.0e-12, abs_tol=1.0e-12
        ):
            raise ValueError(
                f"{equation}: element balance failed for {element} "
                f"({reactant_total:g} != {product_total:g})"
            )


def build_bundle(
    input_path: str | Path,
    *,
    phase: str,
    module_name: str | None = None,
    loader_name: str | None = None,
    kernel_name: str | None = None,
    jacobian_name: str | None = None,
    thermo_loader_name: str | None = None,
    transport_loader_name: str | None = None,
    symbol_prefix: str | None = None,
    chemistry_integrator: str = "implicit",
    description: str | None = None,
    source_file: str | None = None,
) -> dict[str, Any]:
    """Build a deterministic normalized bundle without writing a file."""

    path = Path(input_path)
    if chemistry_integrator not in {"explicit", "implicit"}:
        raise ValueError(
            "chemistry_integrator must be 'explicit' or 'implicit'"
        )
    try:
        source_bytes = path.read_bytes()
    except OSError as exc:
        raise ValueError(f"cannot read Cantera YAML {path}: {exc}") from exc

    cantera = _load_cantera()
    runtime_cantera_version, runtime_cantera_git_commit = (
        _validate_cantera_runtime(cantera)
    )
    yaml_cantera_version, source_units = _yaml_metadata(source_bytes)
    cantera_phase = load_cantera_phase(path, phase, cantera)
    species_names = _validate_species_names(cantera_phase.species_names)
    if not species_names:
        raise ValueError("selected phase contains no species")

    _validate_phase_models(cantera_phase)
    if int(getattr(cantera_phase, "n_reactions", 0)) <= 0:
        raise ValueError("selected phase contains no reactions")

    species = [cantera_phase.species(index) for index in range(len(species_names))]
    thermo = [_thermo_record(item, index) for index, item in enumerate(species)]
    transport = [
        _transport_record(item, index) for index, item in enumerate(species)
    ]
    species_compositions = {
        record["name"]: record["composition"] for record in thermo
    }
    reactions = [
        _reaction_record(
            cantera_phase.reaction(index),
            species_names,
            index,
            species_compositions,
        )
        for index in range(cantera_phase.n_reactions)
    ]

    prefix = _identifier(symbol_prefix or _default_prefix(path), "symbol_prefix")
    module = _identifier(
        module_name or f"{prefix}_mechanism_mod", "module_name"
    )
    loader = _identifier(
        loader_name or f"load_{prefix}_mechanism", "loader_name"
    )
    kernel = _identifier(
        kernel_name or f"{prefix}_production_rates", "kernel_name"
    )
    jacobian = _identifier(
        jacobian_name or f"{prefix}_mass_fraction_jacobian", "jacobian_name"
    )
    thermo_loader = _identifier(
        thermo_loader_name or f"load_{prefix}_thermo", "thermo_loader_name"
    )
    transport_loader = _identifier(
        transport_loader_name or f"load_{prefix}_transport", "transport_loader_name"
    )
    selected_phase = str(getattr(cantera_phase, "name", phase or "gas"))
    bundle = {
        "schema_version": 1,
        "chemistry_integrator": chemistry_integrator,
        "description": description
        or f"Normalized Cantera YAML mechanism for phase {selected_phase}",
        "module_name": module,
        "loader_name": loader,
        "kernel_name": kernel,
        "jacobian_name": jacobian,
        "thermo_loader_name": thermo_loader,
        "transport_loader_name": transport_loader,
        "symbol_prefix": prefix,
        "species": species_names,
        "thermo": thermo,
        "transport": transport,
        "reactions": reactions,
        "source": {
            "format": "Cantera YAML",
            # The content hash makes the basename unambiguous while keeping
            # generated JSON portable across checkout locations.
            "file": source_file or path.name,
            "sha256": hashlib.sha256(source_bytes).hexdigest(),
            "cantera_version": yaml_cantera_version or str(cantera.__version__),
            "phase": selected_phase,
            "runtime_cantera_version": runtime_cantera_version,
            "runtime_cantera_git_commit": runtime_cantera_git_commit,
            "source_units": source_units,
            "target_units": dict(_TARGET_UNITS),
        },
    }
    # This also catches accidental NaN/Infinity values returned by a future
    # Cantera API before anything reaches disk.
    json.dumps(bundle, allow_nan=False)
    return bundle


def write_bundle(bundle: Mapping[str, Any], output_path: str | Path) -> None:
    """Atomically publish a bundle with stable key ordering and a final newline."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(
        bundle,
        ensure_ascii=False,
        allow_nan=False,
        indent=2,
        sort_keys=True,
    )
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(text + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def ingest_cantera_yaml(
    input_path: str | Path, output_path: str | Path, **kwargs: Any
) -> dict[str, Any]:
    """Build and write a bundle; return the in-memory representation."""

    bundle = build_bundle(input_path, **kwargs)
    write_bundle(bundle, output_path)
    return bundle


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ingest a supported Cantera YAML phase into normalized JSON"
    )
    parser.add_argument(
        "--input",
        "--yaml",
        dest="input_path",
        type=Path,
        required=True,
        help="Cantera YAML file",
    )
    parser.add_argument("--output", type=Path, required=True, help="JSON bundle")
    parser.add_argument(
        "--phase",
        "--phase-id",
        dest="phase",
        required=True,
        help="Cantera phase name (required for multi-phase YAML safety)",
    )
    parser.add_argument("--module-name")
    parser.add_argument("--loader-name")
    parser.add_argument("--kernel-name")
    parser.add_argument("--jacobian-name")
    parser.add_argument("--thermo-loader-name")
    parser.add_argument("--transport-loader-name")
    parser.add_argument("--symbol-prefix")
    parser.add_argument(
        "--chemistry-integrator",
        choices=("explicit", "implicit"),
        default="implicit",
        help="generic reactive-flow chemistry integrator policy",
    )
    parser.add_argument("--description")
    parser.add_argument(
        "--source-file",
        help="portable source label (default: input basename)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        ingest_cantera_yaml(
            args.input_path,
            args.output,
            phase=args.phase,
            module_name=args.module_name,
            loader_name=args.loader_name,
            kernel_name=args.kernel_name,
            jacobian_name=args.jacobian_name,
            thermo_loader_name=args.thermo_loader_name,
            transport_loader_name=args.transport_loader_name,
            symbol_prefix=args.symbol_prefix,
            chemistry_integrator=args.chemistry_integrator,
            description=args.description,
            source_file=args.source_file,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ingest_cantera_mechanism: error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
