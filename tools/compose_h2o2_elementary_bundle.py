#!/usr/bin/env python3
"""Compose the fixed seven-species mechanism as a complete selected bundle."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile
from typing import Any


def load_object(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return value


def records_by_name(
    database: dict[str, Any], key: str
) -> dict[str, dict[str, Any]]:
    records = database.get(key)
    if not isinstance(records, list):
        raise ValueError(f"database is missing the {key} array")
    result: dict[str, dict[str, Any]] = {}
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("name"), str):
            raise ValueError(f"database has an invalid {key} record")
        name = record["name"]
        if name in result:
            raise ValueError(f"database has duplicate {key} record {name}")
        result[name] = record
    return result


def compose(
    kinetics: dict[str, Any], database: dict[str, Any]
) -> dict[str, Any]:
    species = kinetics.get("species")
    reactions = kinetics.get("reactions")
    if not isinstance(species, list) or not all(
        isinstance(name, str) for name in species
    ):
        raise ValueError("kinetics bundle has an invalid species array")
    if not isinstance(reactions, list) or not reactions:
        raise ValueError("kinetics bundle has no reactions")
    thermo = records_by_name(database, "thermo")
    transport = records_by_name(database, "transport")
    missing = [
        name for name in species if name not in thermo or name not in transport
    ]
    if missing:
        raise ValueError(f"database is missing selected species: {missing}")
    return {
        "schema_version": 1,
        "chemistry_integrator": "explicit",
        "description": (
            "Fixed seven-species, four-reaction H2/O2 mechanism with "
            "normalized thermodynamics and transport"
        ),
        "module_name": "h2o2_elementary_selected_mechanism_mod",
        "loader_name": "load_h2o2_elementary_selected_mechanism",
        "kernel_name": "h2o2_elementary_selected_production_rates",
        "jacobian_name": "h2o2_elementary_selected_mass_fraction_jacobian",
        "thermo_loader_name": "load_h2o2_elementary_selected_thermo",
        "transport_loader_name": "load_h2o2_elementary_selected_transport",
        "symbol_prefix": "h2o2_elementary_selected",
        "source": (
            "PeleF fixed elementary kinetics combined with pinned H2/O2 "
            "thermodynamic and transport records"
        ),
        "species": species,
        "thermo": [thermo[name] for name in species],
        "transport": [transport[name] for name in species],
        "reactions": reactions,
    }


def write_if_changed(path: Path, content: str) -> None:
    if path.exists() and path.read_text(encoding="utf-8") == content:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kinetics", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    bundle = compose(
        load_object(arguments.kinetics), load_object(arguments.database)
    )
    content = json.dumps(bundle, indent=2, sort_keys=True) + "\n"
    write_if_changed(arguments.output, content)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
