#!/usr/bin/env python3
"""Exercise transactional selected sparse MPI EB restart rejection gates."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import subprocess


MAGIC = "PELEF_REACTIVE_AMR_EB_PATCH_TREE_2D"
CHECKPOINT_NAME = "selected_mpi_eb_restart.chk"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replace_tokens(lines: list[str], index: int, tokens: list[str]) -> list[str]:
    candidate = lines.copy()
    candidate[index] = " ".join(tokens)
    return candidate


def build_mutations(data: bytes) -> dict[str, bytes]:
    lines = data.decode("utf-8").splitlines()
    if not lines or lines[0] != MAGIC or lines[-1] != "END_CHECKPOINT":
        raise AssertionError("source selected MPI EB checkpoint is invalid")
    schema, species_count, nvar, _ = map(int, lines[1].split())
    if schema != 9 or species_count < 1 or nvar != species_count + 5:
        raise AssertionError("source selected MPI EB checkpoint is not schema 9")
    context = 2 + species_count
    if lines[context] != "SELECTED_CONTEXT":
        raise AssertionError("source selected MPI EB context is missing")
    fingerprint = context + 5
    baseline = lines.index("COMPOSITE_BASELINE", fingerprint)
    state_lines: list[int] = []
    for index in range(baseline + 9, len(lines) - 1):
        tokens = lines[index].split()
        if len(tokens) == nvar + 1:
            state_lines.append(index)
    if not state_lines:
        raise AssertionError("source selected MPI EB state payload is missing")

    mutations: dict[str, list[str]] = {}

    candidate = lines.copy()
    header = candidate[1].split()
    header[0] = "8"
    candidate[1] = " ".join(header)
    mutations["schema"] = candidate

    candidate = lines.copy()
    candidate[2], candidate[3] = candidate[3], candidate[2]
    mutations["species_order"] = candidate

    candidate = lines.copy()
    bundle = candidate[context + 1]
    candidate[context + 1] = ("0" if bundle[0] != "0" else "1") + bundle[1:]
    mutations["bundle"] = candidate

    candidate = lines.copy()
    candidate[context + 2] = (
        "explicit" if candidate[context + 2] == "implicit" else "implicit"
    )
    mutations["integrator"] = candidate

    composition = lines[context + 4].split()
    first = float(composition[0])
    second = float(composition[1])
    composition[0] = f"{first + 1.0e-6:.18e}"
    composition[1] = f"{second - 1.0e-6:.18e}"
    mutations["composition"] = replace_tokens(lines, context + 4, composition)

    candidate = lines.copy()
    candidate[fingerprint + 1] = "full_h2o2"
    mutations["fingerprint"] = candidate

    root_geometry = lines[fingerprint + 14].split()
    root_geometry[6] = f"{2.0 * float(root_geometry[6]):.18e}"
    mutations["geometry"] = replace_tokens(lines, fingerprint + 14, root_geometry)

    baseline_values = lines[baseline + 2].split()
    baseline_values[5] = f"{float(baseline_values[5]) + 1.0e-3:.18e}"
    mutations["baseline"] = replace_tokens(lines, baseline + 2, baseline_values)

    for name, component, replacement in [
        ("density", 0, lambda _: "-1.0"),
        ("raw_species", 5, lambda _: "-1.0"),
        ("near_floor_species", 5, lambda _: "-4.0e-11"),
        ("species_closure", 5, lambda value: f"{float(value) + 1.0e-2:.18e}"),
        ("temperature", -1, lambda _: "-1.0"),
    ]:
        candidate = lines.copy()
        for state_line in state_lines:
            values = candidate[state_line].split()
            values[component] = replacement(values[component])
            candidate[state_line] = " ".join(values)
        mutations[name] = candidate

    candidate = lines.copy()
    candidate[-1] = "BROKEN_END_CHECKPOINT"
    mutations["end_marker"] = candidate
    mutations["truncated"] = lines[:-1]
    mutations["trailing"] = [*lines, "TRAILING_CONTENT"]

    return {
        name: ("\n".join(candidate_lines) + "\n").encode("utf-8")
        for name, candidate_lines in mutations.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--launcher", required=True)
    parser.add_argument("--numproc-flag", required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--executable", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--expected", default="Sparse patch-tree restart failed")
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args()

    source = args.checkpoint.read_bytes()
    mutations = build_mutations(source)
    args.work_directory.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["OMPI_MCA_rmaps_base_oversubscribe"] = "1"

    for name, mutated in mutations.items():
        case_directory = args.work_directory / name
        case_directory.mkdir(parents=True, exist_ok=True)
        checkpoint = case_directory / CHECKPOINT_NAME
        output = case_directory / "forbidden.csv"
        if checkpoint.exists():
            checkpoint.unlink()
        if output.exists():
            output.unlink()
        checkpoint.write_bytes(mutated)
        before = digest(mutated)
        command = [
            args.launcher,
            args.numproc_flag,
            str(args.ranks),
            args.executable,
            str(args.input.resolve()),
            output.name,
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
        if completed.returncode == 0:
            raise AssertionError(f"{name}: corrupted checkpoint was accepted")
        if args.expected not in completed.stdout:
            raise AssertionError(
                f"{name}: expected failure text is missing\n{completed.stdout}"
            )
        if output.exists():
            raise AssertionError(f"{name}: failed restart published output")
        if digest(checkpoint.read_bytes()) != before:
            raise AssertionError(f"{name}: failed restart mutated its input")
        print(f"{name}=PASS")

    print(f"selected_mpi_eb_restart_failure_gates={len(mutations)}")


if __name__ == "__main__":
    main()
