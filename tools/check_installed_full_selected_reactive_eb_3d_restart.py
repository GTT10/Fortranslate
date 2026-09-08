#!/usr/bin/env python3
"""Audit the installed full-bundle selected EB 3D schema-2 restart smoke."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re


MAGIC = "PELEF_REACTIVE_EB_3D_CHECKPOINT"
EXPECTED_HEADER = (2, 10, 29, 10, 15)
EXPECTED_ROWS = 10 * 8 * 6
FINAL_TIME = 1.0e-8
BUNDLE_SHA256 = (
    "f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3"
)
INTEGRATOR = "implicit"
SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2")
COMPOSITION = (
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


def sha256(value: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise argparse.ArgumentTypeError("expected a lowercase SHA-256 digest")
    return value


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def summary_value(text: str, label: str) -> float:
    matches = re.findall(
        rf"^{re.escape(label)}\s*([+\-0-9.eEdD]+)\s*$",
        text,
        flags=re.MULTILINE,
    )
    if len(matches) != 1:
        raise AssertionError(f"summary has {len(matches)} values for {label}")
    return float(matches[0].replace("D", "E").replace("d", "e"))


def load_csv(path: Path) -> tuple[bytes, list[dict[str, float]]]:
    data = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [
            {
                name: float(value.replace("D", "E").replace("d", "e"))
                for name, value in row.items()
            }
            for row in csv.DictReader(stream)
        ]
    if len(rows) != EXPECTED_ROWS:
        raise AssertionError(f"{path}: expected {EXPECTED_ROWS} rows, got {len(rows)}")
    if any(not all(math.isfinite(value) for value in row.values()) for row in rows):
        raise AssertionError(f"{path}: nonfinite CSV value")
    counts = {0: 0, 1: 0, 2: 0}
    for index, row in enumerate(rows):
        expected_indices = (
            index % 10 + 1,
            (index // 10) % 8 + 1,
            index // 80 + 1,
        )
        actual_indices = tuple(round(row[name]) for name in ("i", "j", "k"))
        if actual_indices != expected_indices:
            raise AssertionError(f"{path}: row {index} is not x-fastest")
        cell_type = round(row["cell_type"])
        if cell_type not in counts or row["cell_type"] != float(cell_type):
            raise AssertionError(f"{path}: invalid cell type at row {index}")
        counts[cell_type] += 1
        if min(row["rho"], row["pressure"], row["temperature"]) <= 0.0:
            raise AssertionError(f"{path}: nonphysical row {index}")
        closure = sum(row[f"Y_{name}"] for name in SPECIES)
        if abs(closure - 1.0) > 5.0e-12:
            raise AssertionError(f"{path}: species closure failure at row {index}")
    if counts != {0: 144, 1: 48, 2: 288}:
        raise AssertionError(f"{path}: geometry counts {counts}")
    return data, rows


def check_log(path: Path, mode: str) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    required = (
        rf"^Bundle SHA-256:\s+{BUNDLE_SHA256}\s*$",
        rf"^Chemistry integrator:\s+{INTEGRATOR}\s*$",
        r"^Species:\s+10\s*$",
        r"^Reactions:\s+29\s*$",
        r"^Operator sequence:\s+H\s*$",
        r"^Chemistry:\s+F\s*$",
        r"^Molecular transport:\s+F\s*$",
    )
    if any(not re.search(pattern, text, re.MULTILINE) for pattern in required):
        raise AssertionError(f"{path}: selected full-bundle identity is incomplete")
    if mode == "stop":
        if not re.search(r"^Wrote checkpoint:.*step 1, time ", text, re.MULTILINE):
            raise AssertionError(f"{path}: first-step checkpoint record is missing")
        if not re.search(r"^Stopped after checkpoint\s*$", text, re.MULTILINE):
            raise AssertionError(f"{path}: intentional-stop record is missing")
    elif mode == "restart":
        if not re.search(r"^Restored step 1, time ", text, re.MULTILINE):
            raise AssertionError(f"{path}: restored-step record is missing")
        if not re.search(r"^Restart continuation: complete\s*$", text, re.MULTILINE):
            raise AssertionError(f"{path}: continuation record is missing")
    values = {
        "steps": summary_value(text, "Completed steps:"),
        "time": summary_value(text, "Final time:"),
        "minimum_dt": summary_value(text, "Minimum accepted dt:"),
        "invariants": summary_value(text, "Maximum invariant conservation error:"),
        "elements": summary_value(text, "Maximum elemental conservation error:"),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise AssertionError(f"{path}: nonfinite summary")
    if values["invariants"] > 5.0e-11 or values["elements"] > 5.0e-10:
        raise AssertionError(f"{path}: conservation gate exceeded")
    return values


def check_checkpoint(path: Path, expected_sha256: str) -> None:
    data = path.read_bytes()
    if digest(data) != expected_sha256:
        raise AssertionError(f"{path}: checkpoint SHA-256 mismatch")
    lines = data.decode("utf-8").splitlines()
    if lines[0].strip() != MAGIC:
        raise AssertionError("checkpoint magic mismatch")
    header = tuple(int(value) for value in lines[1].split())
    if header != EXPECTED_HEADER:
        raise AssertionError(f"checkpoint header {header} != {EXPECTED_HEADER}")
    if lines[2].strip() != "SELECTED_CONTEXT":
        raise AssertionError("selected context marker mismatch")
    if lines[3].strip() != BUNDLE_SHA256 or lines[4].strip() != INTEGRATOR:
        raise AssertionError("selected mechanism identity mismatch")
    if int(lines[5]) != len(COMPOSITION):
        raise AssertionError("selected composition size mismatch")
    stored_composition = tuple(float(value) for value in lines[6].split())
    if len(stored_composition) != len(COMPOSITION) or any(
        abs(actual - expected) > 5.0e-15
        for actual, expected in zip(stored_composition, COMPOSITION, strict=True)
    ):
        raise AssertionError("selected checkpoint composition mismatch")
    if lines[7].strip() != "SPECIES" or lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError("selected checkpoint body marker mismatch")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--stopped", required=True, type=Path)
    parser.add_argument("--restarted", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference-log", required=True, type=Path)
    parser.add_argument("--stopped-log", required=True, type=Path)
    parser.add_argument("--restarted-log", required=True, type=Path)
    parser.add_argument("--expected-reference-sha256", required=True, type=sha256)
    parser.add_argument("--expected-stopped-sha256", required=True, type=sha256)
    parser.add_argument("--expected-checkpoint-sha256", required=True, type=sha256)
    args = parser.parse_args()

    reference_data, reference_rows = load_csv(args.reference)
    stopped_data, stopped_rows = load_csv(args.stopped)
    restarted_data, restarted_rows = load_csv(args.restarted)
    if digest(reference_data) != args.expected_reference_sha256:
        raise AssertionError("reference CSV SHA-256 mismatch")
    if digest(stopped_data) != args.expected_stopped_sha256:
        raise AssertionError("stopped CSV SHA-256 mismatch")
    if reference_data != restarted_data:
        raise AssertionError("installed selected restart is not byte-exact")
    final_time = reference_rows[0]["time"]
    stopped_time = stopped_rows[0]["time"]
    if abs(final_time - FINAL_TIME) > 5.0e-20:
        raise AssertionError("reference final time mismatch")
    if not 0.0 < stopped_time < final_time:
        raise AssertionError("checkpoint did not precede final time")
    if restarted_rows[0]["time"] != final_time:
        raise AssertionError("restarted final time mismatch")

    reference_log = check_log(args.reference_log, "reference")
    stopped_log = check_log(args.stopped_log, "stop")
    restarted_log = check_log(args.restarted_log, "restart")
    if reference_log != restarted_log:
        raise AssertionError("reference/restart cumulative diagnostics differ")
    if reference_log["steps"] != 4.0 or stopped_log["steps"] != 1.0:
        raise AssertionError("installed smoke did not exercise 4-step continuation")
    if reference_log["time"] != final_time or stopped_log["time"] != stopped_time:
        raise AssertionError("log/CSV time mismatch")

    check_checkpoint(args.checkpoint, args.expected_checkpoint_sha256)
    print(f"rows={EXPECTED_ROWS}")
    print(f"checkpoint_time={stopped_time:.16e}")
    print(f"final_time={final_time:.16e}")
    print("completed_steps=4")
    print(f"reference_csv_sha256={digest(reference_data)}")
    print(f"stopped_csv_sha256={digest(stopped_data)}")
    print(f"checkpoint_sha256={digest(args.checkpoint.read_bytes())}")
    print("installed_full_selected_reactive_eb_3d_restart=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
