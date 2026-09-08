#!/usr/bin/env python3
"""Validate exact planar 3D EB checkpoint/restart continuation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re


NX, NY, NZ = 24, 8, 6
EXPECTED_ROWS = NX * NY * NZ
FINAL_TIME = 2.728e-9
MAGIC = "PELEF_REACTIVE_EB_3D_CHECKPOINT"
X_LENGTH = 2.4e-4
Y_LENGTH = 1.0e-4
Z_LENGTH = 1.0e-4
PLANE_POSITION = 3.95e-5


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(value: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise argparse.ArgumentTypeError(
            "expected a lowercase 64-character SHA-256 digest"
        )
    return value


def summary_value(text: str, label: str) -> float:
    pattern = rf"^{re.escape(label)}\s*([+\-0-9.eEdD]+)\s*$"
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise AssertionError(f"summary has {len(matches)} values for {label}")
    return float(matches[0].replace("D", "E").replace("d", "e"))


def expected_geometry(
    i: int, j: int, k: int
) -> tuple[int, float, tuple[float, float, float], tuple[float, float, float]]:
    dx, dy, dz = X_LENGTH / NX, Y_LENGTH / NY, Z_LENGTH / NZ
    center = ((i - 0.5) * dx, (j - 0.5) * dy, (k - 0.5) * dz)
    if i <= 3:
        return 0, 0.0, center, center
    if i == 4:
        fluid_centroid = (0.5 * (PLANE_POSITION + i * dx), center[1], center[2])
        return 1, (i * dx - PLANE_POSITION) / dx, center, fluid_centroid
    return 2, 1.0, center, center


def load_csv(path: Path) -> tuple[bytes, list[dict[str, float]]]:
    data = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows = [
            {name: float(value.replace("D", "E").replace("d", "e"))
             for name, value in row.items()}
            for row in reader
        ]
    if len(rows) != EXPECTED_ROWS:
        raise AssertionError(f"{path}: expected {EXPECTED_ROWS} rows, got {len(rows)}")
    if any(not all(math.isfinite(value) for value in row.values()) for row in rows):
        raise AssertionError(f"{path}: nonfinite CSV value")
    counts = {0: 0, 1: 0, 2: 0}
    times: set[float] = set()
    for index, row in enumerate(rows):
        expected = (
            index % NX + 1,
            (index // NX) % NY + 1,
            index // (NX * NY) + 1,
        )
        actual = tuple(round(row[name]) for name in ("i", "j", "k"))
        if actual != expected:
            raise AssertionError(f"{path}: row {index} is not x-fastest")
        cell_type, volume_fraction, center, fluid_centroid = expected_geometry(
            *expected
        )
        if row["cell_type"] != float(cell_type):
            raise AssertionError(f"{path}: bad cell type at row {index}")
        if abs(row["volume_fraction"] - volume_fraction) > 3.0e-14:
            raise AssertionError(f"{path}: bad volume fraction at row {index}")
        if max(abs(row[name] - value) for name, value in zip(
            ("x", "y", "z"), center
        )) > 3.0e-18:
            raise AssertionError(f"{path}: bad Cartesian center at row {index}")
        if max(abs(row[name] - value) for name, value in zip(
            ("fluid_centroid_x", "fluid_centroid_y", "fluid_centroid_z"),
            fluid_centroid,
        )) > 3.0e-18:
            raise AssertionError(f"{path}: bad fluid centroid at row {index}")
        counts[cell_type] += 1
        times.add(row["time"])
        if min(row["rho"], row["pressure"], row["temperature"]) <= 0.0:
            raise AssertionError(f"{path}: nonphysical row {index}")
        if abs(sum(row[name] for name in (
            "Y_H2", "Y_H", "Y_O", "Y_O2", "Y_OH", "Y_H2O", "Y_N2"
        )) - 1.0) > 5.0e-12:
            raise AssertionError(f"{path}: species closure failure at row {index}")
    if counts != {0: 144, 1: 48, 2: 960}:
        raise AssertionError(f"{path}: geometry counts {counts}")
    if len(times) != 1:
        raise AssertionError(f"{path}: nonuniform output time")
    return data, rows


def check_log(path: Path, mode: str) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    required = (
        r"^Operator sequence:\s+R\-T\-H\-T\-R\s*$",
        r"^Chemistry:\s+T\s*$",
        r"^Molecular transport:\s+T\s*$",
    )
    if any(not re.search(pattern, text, re.MULTILINE) for pattern in required):
        raise AssertionError(f"{path}: missing coupled-physics summary")
    if mode == "stop":
        if not re.search(r"^Wrote checkpoint:.*step 1, time ", text, re.MULTILINE):
            raise AssertionError(f"{path}: missing first-step checkpoint record")
        if not re.search(r"^Stopped after checkpoint\s*$", text, re.MULTILINE):
            raise AssertionError(f"{path}: missing intentional-stop record")
    elif mode == "restart":
        if not re.search(r"^Restarted from checkpoint:", text, re.MULTILINE):
            raise AssertionError(f"{path}: missing restart record")
        if not re.search(r"^Restored step 1, time ", text, re.MULTILINE):
            raise AssertionError(f"{path}: wrong restored step")
        if not re.search(r"^Restart continuation: complete\s*$", text, re.MULTILINE):
            raise AssertionError(f"{path}: missing continuation completion")
    values = {
        "steps": summary_value(text, "Completed steps:"),
        "time": summary_value(text, "Final time:"),
        "minimum_dt": summary_value(text, "Minimum accepted dt:"),
        "diffusivity": summary_value(text, "Maximum transport diffusivity:"),
        "theta": summary_value(text, "Minimum transport flux theta:"),
        "invariants": summary_value(
            text, "Maximum invariant conservation error:"
        ),
        "elements": summary_value(
            text, "Maximum elemental conservation error:"
        ),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise AssertionError(f"{path}: nonfinite summary")
    if values["invariants"] > 5.0e-11 or values["elements"] > 5.0e-10:
        raise AssertionError(f"{path}: conservation gate exceeded")
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--stopped", type=Path, required=True)
    parser.add_argument("--restarted", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--reference-log", type=Path, required=True)
    parser.add_argument("--stopped-log", type=Path, required=True)
    parser.add_argument("--restarted-log", type=Path, required=True)
    parser.add_argument("--expected-checkpoint-sha256", type=sha256, required=True)
    parser.add_argument("--expected-stopped-sha256", type=sha256, required=True)
    args = parser.parse_args()

    reference_data, reference_rows = load_csv(args.reference)
    stopped_data, stopped_rows = load_csv(args.stopped)
    restarted_data, restarted_rows = load_csv(args.restarted)
    if reference_data != restarted_data:
        raise AssertionError(
            "restart CSV mismatch: "
            f"reference_sha256={digest(reference_data)}, "
            f"restarted_sha256={digest(restarted_data)}"
        )
    if stopped_data == reference_data:
        raise AssertionError("checkpoint-stop output did not precede final output")
    stopped_sha256 = digest(stopped_data)
    if stopped_sha256 != args.expected_stopped_sha256:
        raise AssertionError(
            "checkpoint-stop CSV SHA-256 mismatch: "
            f"{stopped_sha256} != {args.expected_stopped_sha256}"
        )
    final_time = reference_rows[0]["time"]
    stopped_time = stopped_rows[0]["time"]
    if abs(final_time - FINAL_TIME) > 5.0e-20 or not 0.0 < stopped_time < final_time:
        raise AssertionError(
            f"invalid continuation times stopped={stopped_time}, final={final_time}"
        )
    if restarted_rows[0]["time"] != final_time:
        raise AssertionError("restarted CSV time differs from reference")

    reference_log = check_log(args.reference_log, "reference")
    stopped_log = check_log(args.stopped_log, "stop")
    restarted_log = check_log(args.restarted_log, "restart")
    if reference_log != restarted_log:
        raise AssertionError(
            f"cumulative diagnostic mismatch: {reference_log} != {restarted_log}"
        )
    if reference_log["steps"] < 2.0 or stopped_log["steps"] != 1.0:
        raise AssertionError("restart case does not exercise real continuation")
    if reference_log["time"] != final_time or stopped_log["time"] != stopped_time:
        raise AssertionError("log/CSV time mismatch")

    checkpoint_data = args.checkpoint.read_bytes()
    checkpoint_sha256 = digest(checkpoint_data)
    if checkpoint_sha256 != args.expected_checkpoint_sha256:
        raise AssertionError(
            "checkpoint SHA-256 mismatch: "
            f"{checkpoint_sha256} != {args.expected_checkpoint_sha256}"
        )
    lines = checkpoint_data.decode("utf-8").splitlines()
    if len(lines) < 3 or lines[0].strip() != MAGIC:
        raise AssertionError("checkpoint magic is missing")
    header = lines[1].split()
    if not header or int(header[0]) != 1:
        raise AssertionError("checkpoint schema is not version 1")
    if lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError("checkpoint terminal marker is missing")

    print(f"rows={EXPECTED_ROWS}")
    print(f"checkpoint_time={stopped_time:.16e}")
    print(f"final_time={final_time:.16e}")
    print(f"completed_steps={int(reference_log['steps'])}")
    print(f"csv_sha256={digest(reference_data)}")
    print(f"stopped_csv_sha256={stopped_sha256}")
    print(f"checkpoint_sha256={checkpoint_sha256}")
    print("reactive_eb_3d_restart_exact_match=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
