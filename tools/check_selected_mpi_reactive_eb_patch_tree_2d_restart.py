#!/usr/bin/env python3
"""Validate selected schema-9 sparse MPI reactive EB restart parity."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path

from check_selected_mpi_reactive_eb_patch_tree_2d import validate_output


MAGIC = "PELEF_REACTIVE_AMR_EB_PATCH_TREE_2D"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require_digest(label: str, data: bytes, expected: str | None) -> str:
    actual = sha256(data)
    if expected is not None and actual != expected.lower():
        raise AssertionError(f"{label} SHA-256 changed: {actual} != {expected}")
    return actual


def parse_checkpoint(
    path: Path,
    species: list[str],
    bundle_sha256: str,
    integrator: str,
    composition: list[float],
) -> tuple[bytes, float]:
    data = path.read_bytes()
    lines = data.decode("utf-8").splitlines()
    if not lines or lines[0] != MAGIC or lines[-1] != "END_CHECKPOINT":
        raise AssertionError("selected MPI EB checkpoint framing is invalid")
    header = [int(value) for value in lines[1].split()]
    expected_header = [9, len(species), len(species) + 5, 4]
    if header != expected_header:
        raise AssertionError(
            f"selected MPI EB schema-9 header changed: {header} != {expected_header}"
        )
    if lines[2 : 2 + len(species)] != species:
        raise AssertionError("selected MPI EB checkpoint species order changed")

    context = 2 + len(species)
    if lines[context] != "SELECTED_CONTEXT":
        raise AssertionError("selected MPI EB context marker is missing")
    if lines[context + 1] != bundle_sha256:
        raise AssertionError("selected MPI EB bundle SHA-256 changed")
    if lines[context + 2] != integrator:
        raise AssertionError("selected MPI EB integrator changed")
    if int(lines[context + 3]) != len(composition):
        raise AssertionError("selected MPI EB composition size changed")
    stored_composition = [float(value) for value in lines[context + 4].split()]
    if len(stored_composition) != len(composition) or any(
        left != right for left, right in zip(stored_composition, composition)
    ):
        raise AssertionError("selected MPI EB composition changed")

    try:
        baseline = lines.index("COMPOSITE_BASELINE", context + 5)
    except ValueError as error:
        raise AssertionError("selected MPI EB baseline marker is missing") from error
    metadata = lines[baseline - 1].split()
    if len(metadata) != 5:
        raise AssertionError("selected MPI EB checkpoint clock is invalid")
    checkpoint_time = float(metadata[0])
    minimum_dt = float(metadata[1])
    minimum_theta = float(metadata[2])
    steps = int(metadata[3])
    regrids = int(metadata[4])
    if (
        not all(map(math.isfinite, [checkpoint_time, minimum_dt, minimum_theta]))
        or checkpoint_time <= 0.0
        or minimum_dt <= 0.0
        or not 0.0 <= minimum_theta <= 1.0
        or steps != 1
        or regrids < 0
    ):
        raise AssertionError("selected MPI EB checkpoint clock changed")
    nvar = int(lines[baseline + 1])
    baseline_values = [float(value) for value in lines[baseline + 2].split()]
    if nvar != len(species) + 5 or len(baseline_values) != nvar:
        raise AssertionError("selected MPI EB conservation baseline size changed")
    if not all(map(math.isfinite, baseline_values)):
        raise AssertionError("selected MPI EB conservation baseline is nonfinite")
    density = baseline_values[0]
    energy = baseline_values[4]
    species_integrals = baseline_values[5:]
    tolerance = 5.0e-10 * max(1.0, abs(density))
    if (
        density <= 0.0
        or energy <= 0.0
        or min(species_integrals) < -tolerance
        or abs(sum(species_integrals) - density) > tolerance
    ):
        raise AssertionError("selected MPI EB conservation baseline is invalid")
    return data, checkpoint_time


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, action="append", default=[])
    parser.add_argument("--stopped", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--expected-levels", nargs="+", type=int, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--bundle-sha256", required=True)
    parser.add_argument("--integrator", choices=["explicit", "implicit"], required=True)
    parser.add_argument("--composition", nargs="+", type=float, required=True)
    parser.add_argument("--expected-checkpoint-sha256")
    parser.add_argument("--expected-final-sha256")
    parser.add_argument("--expected-stopped-sha256")
    args = parser.parse_args()

    checkpoint, checkpoint_time = parse_checkpoint(
        args.checkpoint,
        args.species,
        args.bundle_sha256,
        args.integrator,
        args.composition,
    )
    checkpoint_digest = require_digest(
        "selected MPI EB checkpoint", checkpoint, args.expected_checkpoint_sha256
    )
    reference = validate_output(
        args.reference, args.species, args.expected_levels, args.final_time
    )
    final_digest = require_digest(
        "selected MPI EB final CSV", reference, args.expected_final_sha256
    )
    for candidate_path in args.candidate:
        candidate = validate_output(
            candidate_path, args.species, args.expected_levels, args.final_time
        )
        if candidate != reference:
            raise AssertionError(
                f"{candidate_path}: selected MPI EB restart parity failed"
            )
    stopped = validate_output(
        args.stopped, args.species, args.expected_levels, checkpoint_time
    )
    stopped_digest = require_digest(
        "selected MPI EB stopped CSV", stopped, args.expected_stopped_sha256
    )
    print(f"checkpoint_sha256={checkpoint_digest}")
    print(f"final_csv_sha256={final_digest}")
    print(f"stopped_csv_sha256={stopped_digest}")
    print("selected_mpi_reactive_eb_patch_tree_2d_restart=PASS")


if __name__ == "__main__":
    main()
