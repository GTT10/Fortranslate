#!/usr/bin/env python3
"""Run a command, mirror its combined output, and save that output to a file."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required after --")

    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    args.log.write_text(completed.stdout, encoding="utf-8")
    sys.stdout.write(completed.stdout)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
