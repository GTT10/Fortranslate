#!/usr/bin/env python3
"""Check public arbitrary-depth EB patch-tree checkpoint/restart parity."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


MAGIC = "PELEF_REACTIVE_AMR_EB_PATCH_TREE_2D"
FINAL_TIME = 3.0e-9
IDENTITY = ("level", "patch", "i", "j")
X_UPPER = 0.012
ROOT_DX = X_UPPER / 12.0
CONSERVATION_PREFIX = "Maximum composite conservation error:"
THETA_PREFIX = "Minimum transport limiter theta:"
COUNTER_PREFIXES = (
    "Chemistry level advances:",
    "Transport level advances:",
    "Hydro level advances:",
)


def load(path: Path) -> dict[tuple[int, int, int, int], dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError(f"{path.name}: empty composite output")
    result: dict[tuple[int, int, int, int], dict[str, str]] = {}
    for row in rows:
        if not all(math.isfinite(float(value)) for value in row.values()):
            raise AssertionError(f"{path.name}: nonfinite value")
        key = tuple(int(row[name]) for name in IDENTITY)
        if key in result:
            raise AssertionError(f"{path.name}: duplicate cell identity {key}")
        result[key] = row
    return result


def compare(
    reference: dict[tuple[int, int, int, int], dict[str, str]],
    restarted: dict[tuple[int, int, int, int], dict[str, str]],
) -> None:
    if reference.keys() != restarted.keys():
        raise AssertionError("reference/restart composite topology mismatch")
    for key in reference:
        expected = reference[key]
        actual = restarted[key]
        if expected.keys() != actual.keys():
            raise AssertionError(f"cell {key}: column mismatch")
        for name in expected:
            lhs = float(actual[name])
            rhs = float(expected[name])
            tolerance = 3.0e-10 * max(1.0, abs(rhs))
            if abs(lhs - rhs) > tolerance:
                raise AssertionError(f"cell {key} {name}: {lhs} != {rhs}")


def check_checkpoint(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != MAGIC:
        raise AssertionError("incorrect patch-tree checkpoint magic")
    if lines[-1] != "END_CHECKPOINT":
        raise AssertionError("incomplete patch-tree checkpoint")
    schema, species, nvar, levels = (int(value) for value in lines[1].split())
    if schema != 7 or species != 7 or nvar <= species or levels != 4:
        raise AssertionError("invalid patch-tree checkpoint header")
    fingerprint = 2 + species
    if lines[fingerprint + 5] != "linear":
        raise AssertionError("checkpoint did not preserve linear prolongation")


def load_diagnostic(path: Path, prefix: str) -> float:
    values = [
        float(line.split(":", maxsplit=1)[1])
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(prefix)
    ]
    if len(values) != 1 or not math.isfinite(values[0]):
        raise AssertionError(f"{path.name}: missing or invalid {prefix}")
    return values[0]


def compare_diagnostic(label: str, expected: float, actual: float) -> None:
    tolerance = 5.0e-12 * max(1.0, abs(expected))
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{label}: {actual} != {expected}")


def load_level_counters(path: Path, prefix: str) -> tuple[int, ...]:
    matches = [
        line.removeprefix(prefix).split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(prefix)
    ]
    if len(matches) != 1 or not matches[0]:
        raise AssertionError(f"{path.name}: missing or invalid {prefix}")
    values = tuple(int(value) for value in matches[0])
    if any(value < 0 for value in values):
        raise AssertionError(f"{path.name}: negative {prefix}")
    return values


def unique_time(rows: dict[tuple[int, int, int, int], dict[str, str]]) -> float:
    times = {float(row["time"]) for row in rows.values()}
    if len(times) != 1:
        raise AssertionError("composite output contains multiple times")
    return next(iter(times))


def check_x_upper_boundary(
    label: str,
    rows: dict[tuple[int, int, int, int], dict[str, str]],
) -> None:
    levels = {key[0] for key in rows}
    boundary_levels = {
        key[0]
        for key, row in rows.items()
        if abs(
            float(row["x"]) + 0.5 * float(row["cell_dx"]) - X_UPPER
        )
        <= 2.0e-14 * ROOT_DX
    }
    if boundary_levels != levels:
        raise AssertionError(
            f"{label}: not every level reaches x-upper {boundary_levels}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--stopped", required=True, type=Path)
    parser.add_argument("--restarted", required=True, nargs="+", type=Path)
    parser.add_argument("--reference-log", type=Path)
    parser.add_argument("--restarted-logs", nargs="+", type=Path)
    args = parser.parse_args()

    if (args.reference_log is None) != (args.restarted_logs is None):
        raise AssertionError("reference and restart logs must be supplied together")
    if (
        args.restarted_logs is not None
        and len(args.restarted) != len(args.restarted_logs)
    ):
        raise AssertionError("restart output/log count mismatch")

    check_checkpoint(args.checkpoint)
    reference = load(args.reference)
    stopped = load(args.stopped)
    restarted_outputs = [load(path) for path in args.restarted]
    labeled_outputs = [("reference", reference), ("stopped", stopped)]
    labeled_outputs.extend(
        (f"restarted[{index}]", rows)
        for index, rows in enumerate(restarted_outputs, start=1)
    )
    for label, rows in labeled_outputs:
        levels = {key[0] for key in rows}
        if levels != {0, 1, 2, 3}:
            raise AssertionError(f"{label}: expected four levels, found {levels}")
        check_x_upper_boundary(label, rows)
    stopped_time = unique_time(stopped)
    if not 0.0 < stopped_time < FINAL_TIME:
        raise AssertionError("checkpoint run did not stop before final time")
    if abs(unique_time(reference) - FINAL_TIME) > 2.0e-20:
        raise AssertionError("reference did not reach final time")
    for index, restarted in enumerate(restarted_outputs, start=1):
        if abs(unique_time(restarted) - FINAL_TIME) > 2.0e-20:
            raise AssertionError(f"restart {index} did not reach final time")
        compare(reference, restarted)
    if args.reference_log is not None:
        reference_conservation = load_diagnostic(
            args.reference_log, CONSERVATION_PREFIX
        )
        reference_theta = load_diagnostic(args.reference_log, THETA_PREFIX)
        reference_counters = {
            prefix: load_level_counters(args.reference_log, prefix)
            for prefix in COUNTER_PREFIXES
        }
        for index, log_path in enumerate(args.restarted_logs, start=1):
            compare_diagnostic(
                f"restart {index} conservation history",
                reference_conservation,
                load_diagnostic(log_path, CONSERVATION_PREFIX),
            )
            compare_diagnostic(
                f"restart {index} transport limiter history",
                reference_theta,
                load_diagnostic(log_path, THETA_PREFIX),
            )
            for prefix in COUNTER_PREFIXES:
                actual = load_level_counters(log_path, prefix)
                if actual != reference_counters[prefix]:
                    raise AssertionError(
                        f"restart {index} {prefix} {actual} != "
                        f"{reference_counters[prefix]}"
                    )
    print("check_reactive_eb_patch_tree_restart_2d: PASS")


if __name__ == "__main__":
    main()
