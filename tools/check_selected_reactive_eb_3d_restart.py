#!/usr/bin/env python3
"""Validate exact selected-mechanism planar EB 3D restart continuation."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import re

import check_reactive_eb_3d_restart as fixed_restart


BUNDLE_SHA256 = (
    "9b5c133fed08c62c810ccc22ff5834e7a650f4ecab670934e8660bd5f51dd15f"
)
INTEGRATOR = "explicit"
SPECIES_COUNT = 7
REACTION_COUNT = 4
TRANSPORT_COUNT = 7
NVAR = 12
EXPECTED_COMPOSITION = (0.29570, 1.0e-5, 1.0e-5, 0.14784, 1.0e-5, 0.0, 0.55643)


def check_selected_log(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    required = (
        rf"^Bundle SHA-256:\s+{BUNDLE_SHA256}\s*$",
        rf"^Chemistry integrator:\s+{INTEGRATOR}\s*$",
        rf"^Species:\s+{SPECIES_COUNT}\s*$",
        rf"^Reactions:\s+{REACTION_COUNT}\s*$",
    )
    if any(not re.search(pattern, text, re.MULTILINE) for pattern in required):
        raise AssertionError(f"{path}: missing selected-mechanism identity")


def check_selected_checkpoint(path: Path) -> bytes:
    data = path.read_bytes()
    lines = data.decode("utf-8").splitlines()
    if len(lines) < 9 or lines[0].strip() != fixed_restart.MAGIC:
        raise AssertionError("selected checkpoint magic is missing")
    header = tuple(int(value) for value in lines[1].split())
    expected_header = (
        2,
        SPECIES_COUNT,
        REACTION_COUNT,
        TRANSPORT_COUNT,
        NVAR,
    )
    if header != expected_header:
        raise AssertionError(
            f"selected checkpoint header {header} != {expected_header}"
        )
    if lines[2].strip() != "SELECTED_CONTEXT":
        raise AssertionError("selected checkpoint context marker is missing")
    if lines[3].strip() != BUNDLE_SHA256:
        raise AssertionError("selected checkpoint bundle SHA-256 is wrong")
    if lines[4].strip() != INTEGRATOR:
        raise AssertionError("selected checkpoint integrator is wrong")
    if int(lines[5].strip()) != SPECIES_COUNT:
        raise AssertionError("selected checkpoint composition size is wrong")
    composition = tuple(float(value) for value in lines[6].split())
    if len(composition) != SPECIES_COUNT or any(
        not math.isfinite(value) for value in composition
    ):
        raise AssertionError("selected checkpoint composition is invalid")
    if any(
        abs(actual - expected) > 5.0e-15
        for actual, expected in zip(composition, EXPECTED_COMPOSITION)
    ):
        raise AssertionError(
            f"selected checkpoint composition {composition} is wrong"
        )
    if abs(sum(composition) - 1.0) > 5.0e-15:
        raise AssertionError("selected checkpoint composition is not normalized")
    if lines[7].strip() != "SPECIES" or lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError("selected checkpoint body markers are invalid")
    return data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed-reference", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--stopped", type=Path, required=True)
    parser.add_argument("--restarted", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--reference-log", type=Path, required=True)
    parser.add_argument("--stopped-log", type=Path, required=True)
    parser.add_argument("--restarted-log", type=Path, required=True)
    args = parser.parse_args()

    reference_data, reference_rows = fixed_restart.load_csv(args.reference)
    fixed_reference_data = args.fixed_reference.read_bytes()
    if reference_data != fixed_reference_data:
        raise AssertionError(
            "selected/fixed reference CSV mismatch: "
            f"selected_sha256={fixed_restart.digest(reference_data)}, "
            f"fixed_sha256={fixed_restart.digest(fixed_reference_data)}"
        )
    stopped_data, stopped_rows = fixed_restart.load_csv(args.stopped)
    restarted_data, restarted_rows = fixed_restart.load_csv(args.restarted)
    if reference_data != restarted_data:
        raise AssertionError(
            "selected restart CSV mismatch: "
            f"reference_sha256={fixed_restart.digest(reference_data)}, "
            f"restarted_sha256={fixed_restart.digest(restarted_data)}"
        )
    if stopped_data == reference_data:
        raise AssertionError("selected checkpoint-stop output is already final")
    final_time = reference_rows[0]["time"]
    stopped_time = stopped_rows[0]["time"]
    if (
        abs(final_time - fixed_restart.FINAL_TIME) > 5.0e-20
        or not 0.0 < stopped_time < final_time
        or restarted_rows[0]["time"] != final_time
    ):
        raise AssertionError(
            f"invalid selected continuation times {stopped_time}, {final_time}"
        )

    reference_log = fixed_restart.check_log(args.reference_log, "reference")
    stopped_log = fixed_restart.check_log(args.stopped_log, "stop")
    restarted_log = fixed_restart.check_log(args.restarted_log, "restart")
    for path in (args.reference_log, args.stopped_log, args.restarted_log):
        check_selected_log(path)
    if reference_log != restarted_log:
        raise AssertionError(
            f"selected cumulative diagnostics differ: "
            f"{reference_log} != {restarted_log}"
        )
    if reference_log["steps"] < 2.0 or stopped_log["steps"] != 1.0:
        raise AssertionError("selected restart does not exercise continuation")
    if reference_log["time"] != final_time or stopped_log["time"] != stopped_time:
        raise AssertionError("selected log/CSV time mismatch")

    checkpoint_data = check_selected_checkpoint(args.checkpoint)
    print(f"rows={fixed_restart.EXPECTED_ROWS}")
    print(f"checkpoint_time={stopped_time:.16e}")
    print(f"final_time={final_time:.16e}")
    print(f"completed_steps={int(reference_log['steps'])}")
    print(f"csv_sha256={fixed_restart.digest(reference_data)}")
    print("selected_fixed_reference_exact_match=PASS")
    print(f"checkpoint_sha256={fixed_restart.digest(checkpoint_data)}")
    print("selected_reactive_eb_3d_restart_exact_match=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
