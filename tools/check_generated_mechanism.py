#!/usr/bin/env python3
"""Verify that the committed generated Fortran mechanism is reproducible."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--committed", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    args = parser.parse_args()

    if args.source is not None:
        bundle = json.loads(args.input.read_text(encoding="utf-8"))
        provenance = bundle.get("source")
        if not isinstance(provenance, dict):
            raise AssertionError("mechanism bundle has no source provenance")
        digest = hashlib.sha256(args.source.read_bytes()).hexdigest()
        if provenance.get("sha256") != digest:
            raise AssertionError("mechanism source hash does not match bundle")
        if provenance.get("file") != args.source.name:
            raise AssertionError("mechanism source filename does not match bundle")

    with tempfile.TemporaryDirectory() as directory:
        generated = Path(directory) / args.committed.name
        subprocess.run(
            [
                sys.executable,
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
