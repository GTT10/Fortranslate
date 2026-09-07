#!/usr/bin/env python3
"""Validate the fixed sparse MPI reactive EB patch-tree checkpoint contract."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


MAGIC = "PELEF_REACTIVE_AMR_EB_PATCH_TREE_2D"
SPECIES = ["H2", "H", "O", "O2", "OH", "H2O", "N2"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--expected-sha256")
    args = parser.parse_args()

    data = args.checkpoint.read_bytes()
    actual_sha256 = hashlib.sha256(data).hexdigest()
    if args.expected_sha256 is not None:
        expected = args.expected_sha256.lower()
        if len(expected) != 64 or any(value not in "0123456789abcdef" for value in expected):
            raise AssertionError("invalid expected checkpoint SHA-256")
        if actual_sha256 != expected:
            raise AssertionError(
                "fixed sparse MPI EB checkpoint SHA-256 changed: "
                f"{actual_sha256} != {expected}"
            )

    lines = data.decode("utf-8").splitlines()
    if not lines or lines[0] != MAGIC or lines[-1] != "END_CHECKPOINT":
        raise AssertionError("fixed sparse MPI EB checkpoint framing changed")
    header = lines[1].split()
    if len(header) != 4:
        raise AssertionError("fixed sparse MPI EB checkpoint header is invalid")
    schema, species_count, variable_count, level_count = map(int, header)
    if (schema, species_count, variable_count, level_count) != (8, 7, 12, 4):
        raise AssertionError("fixed sparse MPI EB schema-8 contract changed")
    if lines[2 : 2 + species_count] != SPECIES:
        raise AssertionError("fixed sparse MPI EB species order changed")

    fingerprint = 2 + species_count
    if lines[fingerprint + 5] != "linear":
        raise AssertionError("fixed sparse MPI EB prolongation fingerprint changed")
    print(f"checkpoint_sha256={actual_sha256}")
    print("fixed_mpi_reactive_eb_patch_tree_2d_checkpoint=PASS")


if __name__ == "__main__":
    main()
