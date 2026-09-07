#!/usr/bin/env python3
"""Run one command, mirror its combined output, and save that output atomically."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "--log", dest="output", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a command is required after --")

    try:
        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        output = result.stdout
    except OSError as error:
        output = f"Could not execute {command[0]}: {error}\n"
        result = subprocess.CompletedProcess(command, 127)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=args.output.parent, prefix=f".{args.output.name}.", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(output)
        os.replace(temporary_name, args.output)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    sys.stdout.write(output)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
