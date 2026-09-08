#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--field",
        choices=(
            "sha256",
            "integrator",
            "composition",
            "species",
            "baseline",
            "density",
            "state_species",
            "near_floor_species",
            "state_closure",
            "temperature",
            "geometry",
        ),
        required=True,
    )
    args = parser.parse_args()

    lines = args.input.read_text(encoding="utf-8").splitlines()
    header = [int(value) for value in lines[1].split()]
    if len(header) != 4 or header[0] != 2:
        raise AssertionError("input is not a selected sparse MPI AMR checkpoint")
    species_count = header[1]
    marker = 2 + species_count
    if lines[marker].strip() != "SELECTED_CONTEXT":
        raise AssertionError("selected context marker is missing")

    if args.field == "species":
        lines[2], lines[3] = lines[3], lines[2]
    elif args.field == "sha256":
        value = lines[marker + 1].strip()
        lines[marker + 1] = ("0" if value[0] != "0" else "1") + value[1:]
    elif args.field == "integrator":
        lines[marker + 2] = (
            "explicit" if lines[marker + 2].strip() == "implicit" else "implicit"
        )
    elif args.field == "composition":
        values = [float(value) for value in lines[marker + 4].split()]
        values[0] += 1.0e-6
        values[1] -= 1.0e-6
        lines[marker + 4] = " ".join(f"{value:.18e}" for value in values)
    elif args.field == "baseline":
        values = [float(value) for value in lines[marker + 7].split()]
        values[5] += 1.0
        lines[marker + 7] = " ".join(f"{value:.18e}" for value in values)
    else:
        payload = marker + 8
        cursor = payload + 4
        level_count = header[3]
        for _ in range(level_count - 1):
            relation = lines[cursor].split()
            relation_patches = int(relation[1])
            cursor += 1 + relation_patches
        first_state = cursor + 1
        values = lines[first_state].split()
        if len(values) != header[2] + 1:
            raise AssertionError("first checkpoint state record is malformed")
        if args.field == "density":
            values[0] = "-1.000000000000000000e+00"
        elif args.field == "state_species":
            values[5] = "-1.000000000000000000e+00"
        elif args.field == "near_floor_species":
            species_transfer = float(values[5]) + 4.0e-11
            values[5] = "-4.000000000000000000e-11"
            values[6] = f"{float(values[6]) + species_transfer:.18e}"
        elif args.field == "state_closure":
            values[5] = f"{float(values[5]) + 1.0:.18e}"
        elif args.field == "temperature":
            values[-1] = "NaN"
        else:
            geometry = lines[payload].split()
            geometry[1] = "NaN"
            lines[payload] = " ".join(geometry)
        if args.field != "geometry":
            lines[first_state] = " ".join(values)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"mutated_selected_mpi_amr_checkpoint={args.output}")


if __name__ == "__main__":
    main()
