#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys

from compare_mpi_sparse_amr_1d import (
    PARITY_TOLERANCE,
    compare_rows,
    read_csv,
    validate_reference,
)

EXPECTED_FINAL_TIME = 5.0e-7
MAGIC = "PELEF_PATCH_TREE_REACTIVE_1D_CHECKPOINT"


def validate_checkpoint(path: str) -> None:
    lines = Path(path).read_text().splitlines()
    if len(lines) < 8 or lines[0].strip() != MAGIC:
        raise AssertionError("invalid sparse MPI AMR checkpoint header")
    header = lines[1].split()
    if len(header) != 4:
        raise AssertionError("invalid sparse MPI AMR checkpoint schema line")
    schema, species_count, variable_count, level_count = map(int, header)
    if schema != 1 or species_count < 1 or variable_count <= species_count:
        raise AssertionError("invalid sparse MPI AMR checkpoint dimensions")
    if level_count < 3:
        raise AssertionError("checkpoint did not preserve the three-level tree")
    metadata = lines[3 + species_count].split()
    if len(metadata) != 5:
        raise AssertionError("invalid sparse MPI AMR checkpoint metadata")
    checkpoint_time = float(metadata[0])
    checkpoint_steps = int(metadata[1])
    if not 0.0 < checkpoint_time < EXPECTED_FINAL_TIME or checkpoint_steps < 1:
        raise AssertionError("checkpoint was not written at an intermediate step")


def main(argv: list[str]) -> None:
    if len(argv) < 4:
        raise SystemExit(
            "usage: compare_mpi_sparse_amr_restart.py checkpoint.chk "
            "reference.csv restarted.csv [...]"
        )
    validate_checkpoint(argv[1])
    names, reference = read_csv(argv[2])
    validate_reference(names, reference, EXPECTED_FINAL_TIME)
    for filename in argv[3:]:
        candidate_names, candidate = read_csv(filename)
        validate_reference(candidate_names, candidate, EXPECTED_FINAL_TIME)
        maximum = compare_rows(names, reference, filename)
        print(f"{filename}: restart_maximum_relative_difference={maximum:.16e}")
        if maximum > PARITY_TOLERANCE:
            raise AssertionError(
                f"{filename}: sparse MPI AMR restart parity failed"
            )


if __name__ == "__main__":
    main(sys.argv)
