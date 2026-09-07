#!/usr/bin/env python3
"""Check selected sparse MPI EB lexical alias rejection before file I/O."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import re
import subprocess


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_assignment(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=).*$")
    replaced, count = pattern.subn(rf'\1 "{value}"', text, count=1)
    if count != 1:
        raise AssertionError(f"missing namelist assignment: {name}")
    return replaced


def add_restart_alias(text: str, value: str) -> str:
    pattern = re.compile(r'(?m)^(\s*checkpoint_file\s*=.*)$')
    replaced, count = pattern.subn(rf'\1\n  restart_file = "{value}"', text, count=1)
    if count != 1:
        raise AssertionError("missing checkpoint assignment for alias case")
    return replaced


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--launcher", required=True)
    parser.add_argument("--numproc-flag", required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--executable", required=True)
    parser.add_argument("--reference-input", type=Path, required=True)
    parser.add_argument("--checkpoint-input", type=Path, required=True)
    parser.add_argument("--restart-input", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args()

    args.work_directory.mkdir(parents=True, exist_ok=True)
    checkpoint_text = args.checkpoint_input.read_text(encoding="utf-8")
    restart_text = args.restart_input.read_text(encoding="utf-8")
    reference_text = args.reference_input.read_text(encoding="utf-8")
    executable = str(Path(args.executable).resolve())
    environment = os.environ.copy()
    environment["OMPI_MCA_rmaps_base_oversubscribe"] = "1"

    specifications: list[tuple[str, str, str, str]] = []
    for name, source, path_kind, expected in [
        (
            "input_checkpoint",
            checkpoint_text,
            "checkpoint_file",
            "Selected sparse MPI EB checkpoint path aliases input or output",
        ),
        (
            "output_checkpoint",
            checkpoint_text,
            "checkpoint_file",
            "Selected sparse MPI EB checkpoint path aliases input or output",
        ),
        (
            "input_restart",
            restart_text,
            "restart_file",
            "Selected sparse MPI EB restart path aliases input or output",
        ),
        (
            "output_restart",
            restart_text,
            "restart_file",
            "Selected sparse MPI EB restart path aliases input or output",
        ),
    ]:
        specifications.append((name, source, path_kind, expected))

    for name, source, path_kind, expected in specifications:
        case_directory = args.work_directory / name
        case_directory.mkdir(parents=True, exist_ok=True)
        input_path = (case_directory / f"{name}.nml").resolve()
        output_path = (case_directory / f"{name}.csv").resolve()
        alias_path = input_path if name.startswith("input_") else output_path
        text = replace_assignment(source, path_kind, str(alias_path))
        input_path.write_text(text, encoding="utf-8")
        before = digest(input_path)
        command = [
            args.launcher,
            args.numproc_flag,
            str(args.ranks),
            executable,
            str(input_path),
            str(output_path),
        ]
        completed = subprocess.run(
            command,
            cwd=case_directory,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=args.timeout,
            check=False,
        )
        if completed.returncode == 0 or expected not in completed.stdout:
            raise AssertionError(f"{name}: alias gate failed\n{completed.stdout}")
        if digest(input_path) != before or output_path.exists():
            raise AssertionError(f"{name}: alias failure mutated or published a file")
        print(f"{name}=PASS")

    case_directory = args.work_directory / "checkpoint_restart"
    case_directory.mkdir(parents=True, exist_ok=True)
    input_path = (case_directory / "checkpoint_restart.nml").resolve()
    output_path = (case_directory / "checkpoint_restart.csv").resolve()
    shared_path = str((case_directory / "shared.chk").resolve())
    text = replace_assignment(checkpoint_text, "checkpoint_file", shared_path)
    text = add_restart_alias(text, shared_path)
    input_path.write_text(text, encoding="utf-8")
    before = digest(input_path)
    command = [
        args.launcher,
        args.numproc_flag,
        str(args.ranks),
        executable,
        str(input_path),
        str(output_path),
    ]
    completed = subprocess.run(
        command,
        cwd=case_directory,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=args.timeout,
        check=False,
    )
    expected = "Selected sparse MPI EB checkpoint and restart paths alias"
    if completed.returncode == 0 or expected not in completed.stdout:
        raise AssertionError(
            f"checkpoint_restart: alias gate failed\n{completed.stdout}"
        )
    if digest(input_path) != before or output_path.exists():
        raise AssertionError("checkpoint_restart: alias failure published a file")
    print("checkpoint_restart=PASS")

    case_directory = args.work_directory / "input_output"
    case_directory.mkdir(parents=True, exist_ok=True)
    input_path = (case_directory / "input_output.nml").resolve()
    input_path.write_text(reference_text, encoding="utf-8")
    before = digest(input_path)
    command = [
        args.launcher,
        args.numproc_flag,
        str(args.ranks),
        executable,
        str(input_path),
        str(input_path),
    ]
    completed = subprocess.run(
        command,
        cwd=case_directory,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=args.timeout,
        check=False,
    )
    expected = "Selected sparse MPI EB input and output paths alias"
    if completed.returncode == 0 or expected not in completed.stdout:
        raise AssertionError(f"input_output: alias gate failed\n{completed.stdout}")
    if digest(input_path) != before:
        raise AssertionError("input_output: alias failure mutated its input")
    print("input_output=PASS")
    print("selected_mpi_eb_alias_gates=6")


if __name__ == "__main__":
    main()
