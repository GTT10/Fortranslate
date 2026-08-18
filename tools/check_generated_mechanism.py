#!/usr/bin/env python3
"""Verify that the committed generated Fortran mechanism is reproducible."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--committed", type=Path, required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as directory:
        generated = Path(directory) / args.committed.name
        subprocess.run(
            [
                "python3",
                str(args.generator),
                "--input",
                str(args.input),
                "--output",
                str(generated),
            ],
            check=True,
        )
        if generated.read_bytes() != args.committed.read_bytes():
            raise AssertionError("committed mechanism is not generator-clean")
    print("generated mechanism: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
