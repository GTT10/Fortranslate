#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


MAGIC = "PELEF_PATCH_TREE_REACTIVE_1D_CHECKPOINT"
SPECIES = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"]
COMPOSITION = [0.29570, 1.0e-5, 1.0e-5, 0.14784, 1.0e-5, 0.0, 0.0, 0.0, 0.0, 0.55643]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_checkpoint(
    path: Path,
    expected_bundle_sha256: str,
    expected_checkpoint_sha256: str,
    final_time: float,
) -> tuple[float, int]:
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected_checkpoint_sha256:
        raise AssertionError("selected sparse MPI AMR checkpoint SHA-256 changed")
    lines = data.decode("utf-8").splitlines()
    if not lines or lines[0].strip() != MAGIC:
        raise AssertionError("invalid selected sparse MPI AMR checkpoint magic")
    header = lines[1].split()
    if len(header) != 4:
        raise AssertionError("invalid selected sparse MPI AMR checkpoint header")
    schema, species_count, variable_count, level_count = map(int, header)
    if (schema, species_count, variable_count) != (2, 10, 15):
        raise AssertionError("selected sparse MPI AMR schema contract changed")
    if level_count < 3:
        raise AssertionError("selected checkpoint did not preserve three levels")
    if [line.strip() for line in lines[2 : 2 + species_count]] != SPECIES:
        raise AssertionError("selected checkpoint species order changed")

    cursor = 2 + species_count
    if lines[cursor].strip() != "SELECTED_CONTEXT":
        raise AssertionError("selected checkpoint context marker is missing")
    if lines[cursor + 1].strip() != expected_bundle_sha256:
        raise AssertionError("selected checkpoint bundle SHA-256 mismatch")
    if lines[cursor + 2].strip() != "implicit":
        raise AssertionError("selected checkpoint integrator mismatch")
    if int(lines[cursor + 3]) != species_count:
        raise AssertionError("selected checkpoint composition size mismatch")
    composition = [float(value) for value in lines[cursor + 4].split()]
    if len(composition) != species_count or any(
        abs(left - right) > 2.0e-18
        for left, right in zip(composition, COMPOSITION)
    ):
        raise AssertionError("selected checkpoint composition changed")
    if lines[cursor + 5].strip() != "COMPOSITE_BASELINE":
        raise AssertionError("selected checkpoint baseline marker is missing")
    if int(lines[cursor + 6]) != variable_count:
        raise AssertionError("selected checkpoint baseline size mismatch")
    baseline = [float(value) for value in lines[cursor + 7].split()]
    if len(baseline) != variable_count or not all(map(math.isfinite, baseline)):
        raise AssertionError("selected checkpoint baseline is invalid")
    tolerance = 5.0e-10 * max(1.0, abs(baseline[0]))
    if baseline[0] <= 0.0 or baseline[4] <= 0.0:
        raise AssertionError("selected checkpoint baseline is nonphysical")
    if min(baseline[5:]) < -tolerance:
        raise AssertionError("selected checkpoint baseline species are negative")
    if abs(sum(baseline[5:]) - baseline[0]) > tolerance:
        raise AssertionError("selected checkpoint baseline violates closure")

    geometry = lines[cursor + 8].split()
    if len(geometry) != 3 or int(geometry[0]) != 24:
        raise AssertionError("selected checkpoint root geometry changed")
    metadata = lines[cursor + 9].split()
    if len(metadata) != 5:
        raise AssertionError("selected checkpoint metadata is invalid")
    checkpoint_time = float(metadata[0])
    checkpoint_steps = int(metadata[1])
    if not 0.0 < checkpoint_time < final_time or checkpoint_steps != 1:
        raise AssertionError("selected checkpoint is not an intermediate state")
    if lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError("selected checkpoint terminator is missing")
    return checkpoint_time, checkpoint_steps


def read_rows(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise AssertionError(f"{path}: missing CSV header")
        names = reader.fieldnames
        rows = [
            {name: float(row[name]) for name in names}
            for row in reader
        ]
    return names, rows


def validate_output(path: Path, expected_time: float) -> None:
    names, rows = read_rows(path)
    required = {
        "level",
        "cell_dx",
        "time",
        "x",
        "rho",
        "pressure",
        "temperature",
        "rhoE",
        "Y_HO2",
        "Y_H2O2",
    }
    if not required.issubset(names):
        raise AssertionError(f"{path}: missing selected sparse AMR columns")
    if len(rows) <= 24 or max(int(row["level"]) for row in rows) < 2:
        raise AssertionError(f"{path}: three-level sparse AMR was not active")
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError(f"{path}: output contains non-finite values")
    tolerance = max(5.0e-18, 5.0e-8 * abs(expected_time))
    if any(abs(row["time"] - expected_time) > tolerance for row in rows):
        raise AssertionError(f"{path}: output time mismatch")
    if min(row["rho"] for row in rows) <= 0.0:
        raise AssertionError(f"{path}: non-positive density")
    if min(row["pressure"] for row in rows) <= 0.0:
        raise AssertionError(f"{path}: non-positive pressure")
    if min(row["temperature"] for row in rows) <= 0.0:
        raise AssertionError(f"{path}: non-positive temperature")
    species_names = [name for name in names if name.startswith("Y_")]
    maximum_closure = max(
        abs(sum(row[name] for name in species_names) - 1.0) for row in rows
    )
    if maximum_closure > 5.0e-10:
        raise AssertionError(f"{path}: species closure failure")
    if min(row[name] for row in rows for name in species_names) < -2.0e-12:
        raise AssertionError(f"{path}: negative species mass fraction")


def require_exact(reference: Path, candidates: list[Path]) -> None:
    expected = reference.read_bytes()
    for candidate in candidates:
        if candidate.read_bytes() != expected:
            raise AssertionError(f"{candidate}: exact sparse MPI AMR parity failed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--fixed-reference", type=Path, required=True)
    parser.add_argument(
        "--selected-reference", type=Path, required=True, action="append"
    )
    parser.add_argument("--restarted", type=Path, required=True, action="append")
    parser.add_argument("--stopped", type=Path, required=True)
    parser.add_argument("--expected-bundle-sha256", required=True)
    parser.add_argument("--expected-checkpoint-sha256", required=True)
    parser.add_argument("--final-time", type=float, required=True)
    args = parser.parse_args()

    checkpoint_time, checkpoint_steps = validate_checkpoint(
        args.checkpoint,
        args.expected_bundle_sha256,
        args.expected_checkpoint_sha256,
        args.final_time,
    )
    final_outputs = [
        args.fixed_reference,
        *args.selected_reference,
        *args.restarted,
    ]
    for path in final_outputs:
        validate_output(path, args.final_time)
    validate_output(args.stopped, checkpoint_time)
    require_exact(args.fixed_reference, [*args.selected_reference, *args.restarted])

    _, rows = read_rows(args.fixed_reference)
    if max(row["Y_HO2"] for row in rows) <= 1.0e-10:
        raise AssertionError("selected sparse MPI AMR HO2 chemistry is inactive")
    if max(row["Y_H2O2"] for row in rows) <= 1.0e-14:
        raise AssertionError("selected sparse MPI AMR H2O2 chemistry is inactive")

    print(f"checkpoint_sha256={digest(args.checkpoint)}")
    print(f"checkpoint_time={checkpoint_time:.16e}")
    print(f"checkpoint_steps={checkpoint_steps}")
    print(f"final_output_sha256={digest(args.fixed_reference)}")
    print("selected_mpi_sparse_amr_restart_1d: PASS")


if __name__ == "__main__":
    main()
