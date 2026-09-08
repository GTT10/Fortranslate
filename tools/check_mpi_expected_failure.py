#!/usr/bin/env python3
"""Require an MPI command to fail with one diagnostic and no output artifact."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--launcher", type=Path, required=True)
    parser.add_argument("--numproc-flag", required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--argument", action="append", default=[])
    parser.add_argument("--expected", required=True)
    parser.add_argument("--preserved-input", type=Path)
    parser.add_argument(
        "--forbidden-output", type=Path, required=True, action="append"
    )
    args = parser.parse_args()

    preserved_contents = None
    if args.preserved_input is not None:
        preserved_contents = args.preserved_input.read_bytes()
    for output in args.forbidden_output:
        output.unlink(missing_ok=True)
    environment = os.environ.copy()
    environment.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    result = subprocess.run(
        [
            str(args.launcher),
            args.numproc_flag,
            str(args.ranks),
            str(args.executable),
            *args.argument,
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=environment,
    )
    sys.stdout.write(result.stdout)
    if result.returncode == 0:
        raise AssertionError("MPI command unexpectedly succeeded")
    if args.expected not in result.stdout:
        raise AssertionError(f"missing expected diagnostic: {args.expected}")
    if args.preserved_input is not None:
        if args.preserved_input.read_bytes() != preserved_contents:
            raise AssertionError("failed MPI command modified its input file")
    residual_outputs = [
        str(output) for output in args.forbidden_output if output.exists()
    ]
    if residual_outputs:
        raise AssertionError(
            "failed MPI command created forbidden output: "
            + ", ".join(residual_outputs)
        )
    print("Expected MPI startup failure: PASS")


if __name__ == "__main__":
    main()
