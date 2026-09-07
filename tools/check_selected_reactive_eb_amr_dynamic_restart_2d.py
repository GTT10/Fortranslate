#!/usr/bin/env python3
"""Validate selected dynamic two-level EB AMR restart and schema-5 identity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re


MAGIC = "PELEF_REACTIVE_EB_AMR_2D_CHECKPOINT"
SCHEMA = 5
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
CHECKPOINT_TIME = 4.9460878884884128e-9
INITIAL_PATCH = (2, 5, 2, 5, 2)
CHECKPOINT_PATCH = (5, 12, 4, 12, 2)
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


def load_csv(
    path: Path,
    expected_rows: int,
    expected_cell_counts: tuple[int, int, int],
) -> tuple[bytes, list[dict[str, float]]]:
    data = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected_header = [*BASE_COLUMNS, *(f"Y_{name}" for name in SPECIES)]
        if reader.fieldnames != expected_header:
            raise AssertionError(f"{path}: CSV schema mismatch")
        rows = [
            {name: float(value) for name, value in source_row.items()}
            for source_row in reader
        ]
    if len(rows) != expected_rows:
        raise AssertionError(f"{path}: expected {expected_rows} rows, got {len(rows)}")
    if any(not all(math.isfinite(value) for value in row.values()) for row in rows):
        raise AssertionError(f"{path}: nonfinite CSV value")
    for row in rows:
        if row["rho"] <= 0.0 or row["pressure"] <= 0.0:
            raise AssertionError(f"{path}: nonpositive state")
        if not 0.0 <= row["volume_fraction"] <= 1.0:
            raise AssertionError(f"{path}: invalid volume fraction")
        if abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0) > 5.0e-12:
            raise AssertionError(f"{path}: species closure mismatch")
    cell_counts = tuple(
        sum(int(row["cell_type"]) == cell_type for row in rows)
        for cell_type in (2, 1, 0)
    )
    if cell_counts != expected_cell_counts:
        raise AssertionError(f"{path}: cell counts {cell_counts} are wrong")
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


def check_selected_log(
    path: Path, expected_steps: int, stopped: bool
) -> tuple[float, float]:
    text = path.read_text(encoding="utf-8")
    required = (
        rf"^Bundle SHA-256:\s+{BUNDLE_SHA256}\s*$",
        rf"^Chemistry integrator:\s+{INTEGRATOR}\s*$",
        rf"^Species:\s+{len(SPECIES)}\s*$",
        r"^Reactions:\s+29\s*$",
        rf"^Completed coarse steps:\s+{expected_steps}\s*$",
        r"^Completed regrids:\s+1\s*$",
        r"^Coarse patch bounds:\s+5 12 4 12\s*$",
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
    conservation_match = re.search(
        r"^Maximum composite conservation error:\s+([0-9.Ee+-]+)\s*$",
        text,
        re.MULTILINE,
    )
    if conservation_match is None:
        raise AssertionError(f"{path}: missing composite conservation error")
    conservation_error = float(conservation_match.group(1))
    if not math.isfinite(conservation_error) or conservation_error < 0.0:
        raise AssertionError(f"{path}: invalid composite conservation error")
    return value, conservation_error


def check_checkpoint(path: Path) -> bytes:
    data = path.read_bytes()
    lines = data.decode("utf-8").splitlines()
    if len(lines) < 200 or lines[0].strip() != MAGIC:
        raise AssertionError("selected dynamic checkpoint magic is missing")
    header = tuple(int(value) for value in lines[1].split())
    if header != (SCHEMA, len(SPECIES), NVAR, 1):
        raise AssertionError(f"checkpoint header {header} is wrong")
    if lines[2].strip() != "SELECTED_CONTEXT":
        raise AssertionError("selected dynamic context marker is missing")
    if lines[3].strip() != BUNDLE_SHA256 or lines[4].strip() != INTEGRATOR:
        raise AssertionError("selected dynamic context identity is wrong")
    if int(lines[5]) != len(SPECIES):
        raise AssertionError("selected dynamic composition size is wrong")
    composition = tuple(float(value) for value in lines[6].split())
    if len(composition) != len(SPECIES) or any(
        not math.isfinite(value) for value in composition
    ):
        raise AssertionError("selected dynamic composition is invalid")
    if any(
        abs(actual - expected) > 5.0e-15
        for actual, expected in zip(composition, EXPECTED_COMPOSITION, strict=True)
    ):
        raise AssertionError("selected dynamic composition is wrong")
    if lines[7].strip() != "DYNAMIC_BASELINE":
        raise AssertionError("selected dynamic baseline marker is missing")
    if int(lines[8]) != NVAR:
        raise AssertionError("selected dynamic baseline size is wrong")
    baseline = tuple(float(value) for value in lines[9].split())
    if (
        len(baseline) != NVAR
        or any(not math.isfinite(value) for value in baseline)
        or baseline[0] <= 0.0
        or baseline[4] <= 0.0
        or any(value < 0.0 for value in baseline[5:])
    ):
        raise AssertionError("selected dynamic baseline is invalid")
    if tuple(line.strip() for line in lines[10 : 10 + len(SPECIES)]) != SPECIES:
        raise AssertionError("selected dynamic species order is wrong")

    flags = tuple(int(value) for value in lines[33].split())
    expected_flags = (0, 0, 1, 1, 1, 1, 1, 0, 0, 2, 1, 1, 4, 4)
    if flags != expected_flags:
        raise AssertionError(f"selected dynamic policy flags {flags} are wrong")
    patch = tuple(int(value) for value in lines[35].split())
    if patch != CHECKPOINT_PATCH or patch == INITIAL_PATCH:
        raise AssertionError(f"selected dynamic checkpoint patch {patch} is wrong")
    metadata = lines[36].split()
    if len(metadata) != 5:
        raise AssertionError("selected dynamic time metadata is malformed")
    checkpoint_time = float(metadata[0])
    minimum_dt = float(metadata[1])
    steps = int(metadata[2])
    regrids = int(metadata[3])
    if (
        checkpoint_time != CHECKPOINT_TIME
        or minimum_dt != CHECKPOINT_TIME
        or steps != 1
        or regrids != 1
    ):
        raise AssertionError(f"selected dynamic metadata {metadata} is wrong")
    if tuple(int(value) for value in lines[37].split()) != (12, 12):
        raise AssertionError("selected dynamic coarse dimensions are wrong")
    fine_dimensions_index = 38 + 12 * 12
    if tuple(int(value) for value in lines[fine_dimensions_index].split()) != (
        16,
        18,
    ):
        raise AssertionError("selected dynamic fine dimensions are wrong")
    if lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError("selected dynamic checkpoint end marker is missing")
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
    level_contracts = {
        "coarse": (144, (134, 7, 3)),
        "fine": (288, (288, 0, 0)),
    }
    for level, (row_count, cell_counts) in level_contracts.items():
        for stage in ("fixed", "reference", "stopped", "restarted"):
            datasets[f"{stage}_{level}"] = load_csv(
                getattr(args, f"{stage}_{level}"), row_count, cell_counts
            )
        fixed = datasets[f"fixed_{level}"][0]
        reference = datasets[f"reference_{level}"][0]
        stopped = datasets[f"stopped_{level}"][0]
        restarted = datasets[f"restarted_{level}"][0]
        require_exact(fixed, reference, f"fixed/selected dynamic {level} mismatch")
        require_exact(reference, restarted, f"selected dynamic {level} restart mismatch")
        if stopped == reference:
            raise AssertionError(f"selected dynamic {level} stop output is final")
        check_hash(
            reference,
            getattr(args, f"expected_reference_{level}_sha256"),
            f"dynamic reference {level} SHA-256",
        )
        check_hash(
            stopped,
            getattr(args, f"expected_stopped_{level}_sha256"),
            f"dynamic stopped {level} SHA-256",
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
        or stopped_time != CHECKPOINT_TIME
        or restarted_time != FINAL_TIME
        or fine_times != (reference_time, stopped_time, restarted_time)
    ):
        raise AssertionError("selected dynamic restart times are inconsistent")

    reference_log_time, reference_conservation_error = check_selected_log(
        args.reference_log, 3, False
    )
    if reference_log_time != FINAL_TIME:
        raise AssertionError("selected dynamic reference log time mismatch")
    stopped_log_time, _ = check_selected_log(args.stopped_log, 1, True)
    if stopped_log_time != CHECKPOINT_TIME:
        raise AssertionError("selected dynamic stopped log time mismatch")
    restarted_log_time, restarted_conservation_error = check_selected_log(
        args.restarted_log, 3, False
    )
    if restarted_log_time != FINAL_TIME:
        raise AssertionError("selected dynamic restarted log time mismatch")
    if restarted_conservation_error != reference_conservation_error:
        raise AssertionError(
            "selected dynamic cumulative conservation diagnostic mismatch: "
            f"{restarted_conservation_error} != {reference_conservation_error}"
        )

    checkpoint_data = check_checkpoint(args.checkpoint)
    check_hash(
        checkpoint_data,
        args.expected_checkpoint_sha256,
        "selected dynamic checkpoint SHA-256",
    )
    print("coarse_rows=144")
    print("fine_rows=288")
    print(f"checkpoint_time={CHECKPOINT_TIME:.16e}")
    print("checkpoint_patch=5:12,4:12")
    print("completed_regrids=1")
    print(f"maximum_composite_conservation_error={reference_conservation_error:.16e}")
    print(f"reference_coarse_sha256={digest(datasets['reference_coarse'][0])}")
    print(f"reference_fine_sha256={digest(datasets['reference_fine'][0])}")
    print(f"checkpoint_sha256={digest(checkpoint_data)}")
    print("selected_dynamic_reactive_eb_amr_2d_restart_exact_match=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
