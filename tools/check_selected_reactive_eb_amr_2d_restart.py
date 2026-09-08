#!/usr/bin/env python3
"""Validate selected static two-level reactive EB AMR 2D restart parity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re


MAGIC = "PELEF_REACTIVE_EB_AMR_2D_CHECKPOINT"
SCHEMA = 4
BUNDLE_SHA256 = (
    "f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3"
)
INTEGRATOR = "implicit"
SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2")
EXPECTED_COMPOSITION = (
    0.29570,
    1.0e-5,
    1.0e-5,
    0.14784,
    1.0e-5,
    0.0,
    0.0,
    0.0,
    0.0,
    0.55643,
)
NVAR = 5 + len(SPECIES)
FINAL_TIME = 1.0e-8
BASE_COLUMNS = (
    "time",
    "x",
    "y",
    "volume_fraction",
    "cell_type",
    "boundary_length",
    "boundary_normal_x",
    "boundary_normal_y",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_csv(path: Path, expected_rows: int) -> tuple[bytes, list[dict[str, float]]]:
    data = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected_header = [*BASE_COLUMNS, *(f"Y_{name}" for name in SPECIES)]
        if reader.fieldnames != expected_header:
            raise AssertionError(f"{path}: CSV schema mismatch")
        rows: list[dict[str, float]] = []
        for source_row in reader:
            row = {name: float(value) for name, value in source_row.items()}
            if not all(math.isfinite(value) for value in row.values()):
                raise AssertionError(f"{path}: nonfinite CSV value")
            rows.append(row)
    if len(rows) != expected_rows:
        raise AssertionError(f"{path}: expected {expected_rows} rows, got {len(rows)}")
    for row in rows:
        if abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0) > 5.0e-12:
            raise AssertionError(f"{path}: species closure mismatch")
    return data, rows


def require_exact(left: bytes, right: bytes, label: str) -> None:
    if left != right:
        raise AssertionError(
            f"{label}: left_sha256={digest(left)}, right_sha256={digest(right)}"
        )


def check_hash(data: bytes, expected: str, label: str) -> None:
    actual = digest(data)
    if actual != expected:
        raise AssertionError(f"{label}: {actual} != {expected}")


def check_selected_log(path: Path, expected_steps: int, stopped: bool) -> float:
    text = path.read_text(encoding="utf-8")
    required = (
        rf"^Bundle SHA-256:\s+{BUNDLE_SHA256}\s*$",
        rf"^Chemistry integrator:\s+{INTEGRATOR}\s*$",
        rf"^Species:\s+{len(SPECIES)}\s*$",
        r"^Reactions:\s+29\s*$",
        rf"^Completed coarse steps:\s+{expected_steps}\s*$",
        rf"^Stopped after checkpoint:\s+{'T' if stopped else 'F'}\s*$",
    )
    for pattern in required:
        if re.search(pattern, text, re.MULTILINE) is None:
            raise AssertionError(f"{path}: missing log identity {pattern}")
    match = re.search(r"^Final time:\s+([0-9.Ee+-]+)\s*$", text, re.MULTILINE)
    if match is None:
        raise AssertionError(f"{path}: missing final time")
    value = float(match.group(1))
    if not math.isfinite(value):
        raise AssertionError(f"{path}: nonfinite final time")
    return value


def check_checkpoint(path: Path, stopped_time: float) -> bytes:
    data = path.read_bytes()
    lines = data.decode("utf-8").splitlines()
    if len(lines) < 40 or lines[0].strip() != MAGIC:
        raise AssertionError("selected EB AMR checkpoint magic is missing")
    header = tuple(int(value) for value in lines[1].split())
    expected_header = (SCHEMA, len(SPECIES), NVAR, 1)
    if header != expected_header:
        raise AssertionError(f"checkpoint header {header} != {expected_header}")
    if lines[2].strip() != "SELECTED_CONTEXT":
        raise AssertionError("selected checkpoint context marker is missing")
    if lines[3].strip() != BUNDLE_SHA256:
        raise AssertionError("selected checkpoint bundle SHA-256 is wrong")
    if lines[4].strip() != INTEGRATOR:
        raise AssertionError("selected checkpoint integrator is wrong")
    if int(lines[5]) != len(SPECIES):
        raise AssertionError("selected checkpoint composition size is wrong")
    composition = tuple(float(value) for value in lines[6].split())
    if len(composition) != len(SPECIES) or any(
        not math.isfinite(value) for value in composition
    ):
        raise AssertionError("selected checkpoint composition is invalid")
    if any(
        abs(actual - expected) > 5.0e-15
        for actual, expected in zip(composition, EXPECTED_COMPOSITION, strict=True)
    ):
        raise AssertionError(f"selected checkpoint composition {composition} is wrong")
    if abs(sum(composition) - 1.0) > 5.0e-15:
        raise AssertionError("selected checkpoint composition is not normalized")
    if tuple(line.strip() for line in lines[7 : 7 + len(SPECIES)]) != SPECIES:
        raise AssertionError("selected checkpoint species order is wrong")

    metadata_index = 7 + len(SPECIES) + 16
    metadata = lines[metadata_index].split()
    if len(metadata) != 5:
        raise AssertionError("selected checkpoint time metadata is malformed")
    checkpoint_time = float(metadata[0])
    minimum_dt = float(metadata[1])
    steps = int(metadata[2])
    regrids = int(metadata[3])
    if (
        checkpoint_time != stopped_time
        or minimum_dt <= 0.0
        or minimum_dt > checkpoint_time
        or steps != 1
        or regrids != 0
    ):
        raise AssertionError(f"selected checkpoint metadata is wrong: {metadata}")
    if lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError("selected checkpoint end marker is missing")
    return data


def main() -> int:
    parser = argparse.ArgumentParser()
    for level in ("coarse", "fine"):
        parser.add_argument(f"--fixed-{level}", required=True, type=Path)
        parser.add_argument(f"--reference-{level}", required=True, type=Path)
        parser.add_argument(f"--stopped-{level}", required=True, type=Path)
        parser.add_argument(f"--restarted-{level}", required=True, type=Path)
        parser.add_argument(f"--expected-reference-{level}-sha256", required=True)
        parser.add_argument(f"--expected-stopped-{level}-sha256", required=True)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--expected-checkpoint-sha256", required=True)
    parser.add_argument("--reference-log", required=True, type=Path)
    parser.add_argument("--stopped-log", required=True, type=Path)
    parser.add_argument("--restarted-log", required=True, type=Path)
    args = parser.parse_args()

    datasets: dict[str, tuple[bytes, list[dict[str, float]]]] = {}
    for level, count in (("coarse", 64), ("fine", 144)):
        for stage in ("fixed", "reference", "stopped", "restarted"):
            datasets[f"{stage}_{level}"] = load_csv(
                getattr(args, f"{stage}_{level}"), count
            )
        fixed = datasets[f"fixed_{level}"][0]
        reference = datasets[f"reference_{level}"][0]
        stopped = datasets[f"stopped_{level}"][0]
        restarted = datasets[f"restarted_{level}"][0]
        require_exact(fixed, reference, f"fixed/selected {level} reference mismatch")
        require_exact(reference, restarted, f"selected {level} restart mismatch")
        if stopped == reference:
            raise AssertionError(f"selected {level} checkpoint-stop output is final")
        check_hash(
            reference,
            getattr(args, f"expected_reference_{level}_sha256"),
            f"reference {level} SHA-256",
        )
        check_hash(
            stopped,
            getattr(args, f"expected_stopped_{level}_sha256"),
            f"stopped {level} SHA-256",
        )

    reference_time = datasets["reference_coarse"][1][0]["time"]
    stopped_time = datasets["stopped_coarse"][1][0]["time"]
    restarted_time = datasets["restarted_coarse"][1][0]["time"]
    fine_times = (
        datasets["reference_fine"][1][0]["time"],
        datasets["stopped_fine"][1][0]["time"],
        datasets["restarted_fine"][1][0]["time"],
    )
    if (
        reference_time != FINAL_TIME
        or restarted_time != FINAL_TIME
        or not 0.0 < stopped_time < FINAL_TIME
        or fine_times != (reference_time, stopped_time, restarted_time)
    ):
        raise AssertionError("selected EB AMR restart times are inconsistent")

    if check_selected_log(args.reference_log, 2, False) != FINAL_TIME:
        raise AssertionError("selected reference log final time mismatch")
    if check_selected_log(args.stopped_log, 1, True) != stopped_time:
        raise AssertionError("selected stopped log time mismatch")
    if check_selected_log(args.restarted_log, 2, False) != FINAL_TIME:
        raise AssertionError("selected restarted log final time mismatch")

    checkpoint_data = check_checkpoint(args.checkpoint, stopped_time)
    check_hash(
        checkpoint_data,
        args.expected_checkpoint_sha256,
        "selected checkpoint SHA-256",
    )
    print("coarse_rows=64")
    print("fine_rows=144")
    print(f"checkpoint_time={stopped_time:.16e}")
    print(f"final_time={FINAL_TIME:.16e}")
    print(f"reference_coarse_sha256={digest(datasets['reference_coarse'][0])}")
    print(f"reference_fine_sha256={digest(datasets['reference_fine'][0])}")
    print(f"checkpoint_sha256={digest(checkpoint_data)}")
    print("selected_reactive_eb_amr_2d_restart_exact_match=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
