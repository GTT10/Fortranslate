#!/usr/bin/env python3
"""Require a command to fail with one diagnostic and no output artifact."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--expected", required=True)
    parser.add_argument(
        "--forbidden-output", type=Path, required=True, action="append"
    )
    args = parser.parse_args()

    input_contents = args.input.read_bytes()
    for output in args.forbidden_output:
        output.unlink(missing_ok=True)
    result = subprocess.run(
        [str(args.executable), str(args.input)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    sys.stdout.write(result.stdout)
    if result.returncode == 0:
        raise AssertionError("command unexpectedly succeeded")
    if args.expected not in result.stdout:
        raise AssertionError(f"missing expected diagnostic: {args.expected}")
    if args.input.read_bytes() != input_contents:
        raise AssertionError("failed command modified its input file")
    residual_outputs = [
        str(output) for output in args.forbidden_output if output.exists()
    ]
    if residual_outputs:
        raise AssertionError(
            "failed command created forbidden output: "
            + ", ".join(residual_outputs)
        )
    print("Expected startup failure: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
