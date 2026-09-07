#!/usr/bin/env python3
"""Validate selected dynamic-parent three-level EB AMR restart parity."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re

from check_selected_reactive_eb_amr_three_level_restart_2d import (
    BASE_COLUMNS,
    BUNDLE_SHA256,
    EXPECTED_COMPOSITION,
    INTEGRATOR,
    NVAR,
    SPECIES,
    check_hash,
    close,
    digest,
    parse_floats,
    parse_ints,
    require_exact,
)


MAGIC = "PELEF_REACTIVE_EB_AMR_DYNAMIC_THREE_LEVEL_2D_CHECKPOINT"
SCHEMA = 5
FINAL_TIME = 1.0e-8
ROOT_SEED = (2, 11, 2, 11, 2)
FINEST_SEED = (6, 9, 6, 9, 2)
ROOT_PATCH = (5, 10, 3, 10, 2)
FINEST_PATCH = (3, 10, 3, 14, 2)

LEVEL_CONTRACTS = {
    "root": {
        "rows": 144,
        "counts": (84, 12, 48),
        "shape": (12, 12),
        "domain": (0.0, 0.012, 0.0, 0.012),
    },
    "middle": {
        "rows": 192,
        "counts": (176, 16, 0),
        "shape": (12, 16),
        "domain": (0.004, 0.010, 0.002, 0.010),
    },
    "finest": {
        "rows": 384,
        "counts": (384, 0, 0),
        "shape": (16, 24),
        "domain": (0.005, 0.009, 0.003, 0.009),
    },
}


def read_csv(path: Path, level: str) -> tuple[bytes, list[dict[str, str]]]:
    contract = LEVEL_CONTRACTS[level]
    data = path.read_bytes()
    expected_header = [*BASE_COLUMNS, *(f"Y_{name}" for name in SPECIES)]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_header:
            raise AssertionError(f"{path}: dynamic three-level CSV schema mismatch")
        rows = list(reader)
    if len(rows) != contract["rows"]:
        raise AssertionError(
            f"{path}: expected {contract['rows']} rows, got {len(rows)}"
        )
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise AssertionError(f"{path}: incomplete CSV row")
    return data, rows


def check_csv(path: Path, level: str) -> tuple[bytes, list[dict[str, str]], float, float]:
    data, rows = read_csv(path, level)
    contract = LEVEL_CONTRACTS[level]
    nx, ny = contract["shape"]
    x_lower, x_upper, y_lower, y_upper = contract["domain"]
    dx = (x_upper - x_lower) / nx
    dy = (y_upper - y_lower) / ny
    expected_header = [*BASE_COLUMNS, *(f"Y_{name}" for name in SPECIES)]
    cell_counts = {0: 0, 1: 0, 2: 0}
    maximum_closure_error = 0.0

    for index, row in enumerate(rows):
        parse_floats([row[column] for column in expected_header], str(path))
        time = float(row["time"])
        expected_x = x_lower + (index % nx + 0.5) * dx
        expected_y = y_lower + (index // nx + 0.5) * dy
        if not close(float(row["x"]), expected_x, 5.0e-14):
            raise AssertionError(f"{path}: x-coordinate order mismatch at row {index}")
        if not close(float(row["y"]), expected_y, 5.0e-14):
            raise AssertionError(f"{path}: y-coordinate order mismatch at row {index}")

        cell_type_value = float(row["cell_type"])
        cell_type = int(cell_type_value)
        if cell_type_value != cell_type or cell_type not in cell_counts:
            raise AssertionError(f"{path}: invalid cell type at row {index}")
        cell_counts[cell_type] += 1
        volume_fraction = float(row["volume_fraction"])
        if not 0.0 <= volume_fraction <= 1.0:
            raise AssertionError(f"{path}: volume fraction outside [0,1]")
        if cell_type == 0 and volume_fraction != 0.0:
            raise AssertionError(f"{path}: covered cell has nonzero volume")
        if cell_type == 1 and not 0.0 < volume_fraction < 1.0:
            raise AssertionError(f"{path}: cut cell lacks partial volume")
        if cell_type == 2 and volume_fraction != 1.0:
            raise AssertionError(f"{path}: regular cell lacks unit volume")
        if float(row["boundary_length"]) < 0.0:
            raise AssertionError(f"{path}: negative boundary length")
        if (
            float(row["rho"]) <= 0.0
            or float(row["pressure"]) <= 0.0
            or float(row["temperature"]) <= 0.0
            or float(row["rhoE"]) <= 0.0
        ):
            raise AssertionError(f"{path}: nonpositive thermodynamic state")

        fractions = [float(row[f"Y_{name}"]) for name in SPECIES]
        if any(value < -1.0e-12 for value in fractions):
            raise AssertionError(f"{path}: negative species mass fraction")
        closure_error = abs(sum(fractions) - 1.0)
        maximum_closure_error = max(maximum_closure_error, closure_error)
        if closure_error > 5.0e-12:
            raise AssertionError(f"{path}: species closure mismatch")
        if not math.isfinite(time):
            raise AssertionError(f"{path}: nonfinite time")

    expected_counts = {
        2: contract["counts"][0],
        1: contract["counts"][1],
        0: contract["counts"][2],
    }
    if cell_counts != expected_counts:
        raise AssertionError(f"{path}: cell counts {cell_counts} != {expected_counts}")
    times = {float(row["time"]) for row in rows}
    if len(times) != 1:
        raise AssertionError(f"{path}: level contains multiple output times")
    return data, rows, next(iter(times)), maximum_closure_error


def check_log(path: Path, stopped: bool) -> tuple[float, int, int, str]:
    text = path.read_text(encoding="utf-8")
    expected_flag = "T" if stopped else "F"
    required = (
        rf"^Bundle SHA-256:\s+{BUNDLE_SHA256}\s*$",
        rf"^Chemistry integrator:\s+{INTEGRATOR}\s*$",
        rf"^Species:\s+{len(SPECIES)}\s*$",
        r"^Reactions:\s+29\s*$",
        r"^Completed regrids:\s+1\s*$",
        r"^Root patch bounds:\s+5 10 3 10\s*$",
        r"^Middle patch bounds:\s+3 10 3 14\s*$",
        r"^Refinement ratio:\s+2\s*$",
        r"^Chemistry:\s+T\s*$",
        rf"^Stopped after checkpoint:\s+{expected_flag}\s*$",
    )
    for pattern in required:
        if re.search(pattern, text, re.MULTILINE) is None:
            raise AssertionError(f"{path}: missing log identity {pattern}")

    steps_match = re.search(r"^Completed coarse steps:\s+([0-9]+)\s*$", text, re.MULTILINE)
    regrids_match = re.search(r"^Completed regrids:\s+([0-9]+)\s*$", text, re.MULTILINE)
    time_match = re.search(r"^Final time:\s+([0-9.Ee+-]+)\s*$", text, re.MULTILINE)
    conservation_match = re.search(
        r"^Maximum composite conservation error:\s+([0-9.Ee+-]+)\s*$",
        text,
        re.MULTILINE,
    )
    if (
        steps_match is None
        or regrids_match is None
        or time_match is None
        or conservation_match is None
    ):
        raise AssertionError(f"{path}: missing lifecycle diagnostic")
    steps = int(steps_match.group(1))
    regrids = int(regrids_match.group(1))
    final_time = float(time_match.group(1))
    conservation_text = conservation_match.group(1)
    conservation_value = float(conservation_text)
    if steps < 1 or regrids != 1 or not math.isfinite(final_time):
        raise AssertionError(f"{path}: invalid lifecycle diagnostic")
    if not math.isfinite(conservation_value) or conservation_value < 0.0:
        raise AssertionError(f"{path}: invalid conservation diagnostic")
    return final_time, steps, regrids, conservation_text


def checkpoint_field(
    lines: list[str], index: int, nx: int, ny: int, label: str
) -> int:
    dimensions = parse_ints(lines[index].split(), f"{label} dimensions")
    if dimensions != (nx, ny):
        raise AssertionError(f"{label}: dimensions {dimensions} != {(nx, ny)}")
    index += 1
    for row_index in range(nx * ny):
        values = parse_floats(lines[index].split(), f"{label} row {row_index}")
        if len(values) != NVAR + 1:
            raise AssertionError(
                f"{label} row {row_index}: {len(values)} values != {NVAR + 1}"
            )
        if values[0] <= 0.0 or values[4] <= 0.0 or values[-1] <= 0.0:
            raise AssertionError(f"{label} row {row_index}: nonpositive state")
        species_densities = values[5:NVAR]
        negative_tolerance = 5.0e-12 * values[0]
        if any(value < -negative_tolerance for value in species_densities):
            raise AssertionError(f"{label} row {row_index}: negative species density")
        closure_error = abs(sum(species_densities) - values[0]) / values[0]
        if closure_error > 5.0e-12:
            raise AssertionError(
                f"{label} row {row_index}: species closure error {closure_error}"
            )
        index += 1
    return index


def check_checkpoint(path: Path) -> tuple[bytes, float, int, int]:
    data = path.read_bytes()
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise AssertionError(f"{path}: checkpoint is not UTF-8 text") from error
    if len(lines) < 750 or lines[0].strip() != MAGIC:
        raise AssertionError(f"{path}: selected dynamic checkpoint magic is missing")
    if lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError(f"{path}: incomplete selected dynamic checkpoint")

    header = parse_ints(lines[1].split(), "checkpoint header")
    if header != (SCHEMA, len(SPECIES), NVAR):
        raise AssertionError(f"{path}: checkpoint header {header} is wrong")
    if lines[2].strip() != "SELECTED_CONTEXT":
        raise AssertionError(f"{path}: selected context marker is missing")
    if lines[3].strip() != BUNDLE_SHA256:
        raise AssertionError(f"{path}: selected bundle SHA-256 is wrong")
    if lines[4].strip() != INTEGRATOR:
        raise AssertionError(f"{path}: selected chemistry integrator is wrong")
    if parse_ints(lines[5].split(), "checkpoint composition count") != (len(SPECIES),):
        raise AssertionError(f"{path}: selected composition count is wrong")
    composition = parse_floats(lines[6].split(), "checkpoint composition")
    if len(composition) != len(SPECIES) or any(
        abs(actual - expected) > 5.0e-14
        for actual, expected in zip(composition, EXPECTED_COMPOSITION, strict=True)
    ):
        raise AssertionError(f"{path}: selected composition is wrong")
    if abs(sum(composition) - 1.0) > 5.0e-12 or any(value < 0.0 for value in composition):
        raise AssertionError(f"{path}: selected composition is not physical")

    if lines[7].strip() != "COMPOSITE_BASELINE":
        raise AssertionError(f"{path}: composite baseline marker is missing")
    if parse_ints(lines[8].split(), "checkpoint baseline count") != (NVAR,):
        raise AssertionError(f"{path}: composite baseline count is wrong")
    baseline = parse_floats(lines[9].split(), "checkpoint baseline")
    if len(baseline) != NVAR or baseline[0] <= 0.0 or baseline[4] <= 0.0:
        raise AssertionError(f"{path}: composite baseline is invalid")
    if any(value < 0.0 for value in baseline[5:]):
        raise AssertionError(f"{path}: composite baseline species is negative")
    baseline_closure_error = abs(sum(baseline[5:]) - baseline[0]) / baseline[0]
    if baseline_closure_error > 5.0e-12:
        raise AssertionError(
            f"{path}: composite baseline closure error {baseline_closure_error}"
        )

    species_start = 10
    species_end = species_start + len(SPECIES)
    if tuple(line.strip() for line in lines[species_start:species_end]) != SPECIES:
        raise AssertionError(f"{path}: selected checkpoint species order is wrong")

    index = species_end
    if lines[index].strip() != "plane":
        raise AssertionError(f"{path}: checkpoint geometry is not plane")
    index += 1
    if parse_ints(lines[index].split(), "checkpoint root dimensions") != (12, 12, 2):
        raise AssertionError(f"{path}: checkpoint root dimensions are wrong")
    index += 1
    domain = parse_floats(lines[index].split(), "checkpoint domain")
    if len(domain) != 4 or any(
        not close(a, b, 5.0e-14)
        for a, b in zip(domain, (0.0, 0.012, 0.0, 0.012), strict=True)
    ):
        raise AssertionError(f"{path}: checkpoint domain is wrong")
    index += 1
    geometry = parse_floats(lines[index].split(), "checkpoint plane geometry")
    expected_geometry = (1.0, 0.0, 0.00437, 0.0, 0.0, 1.0)
    if len(geometry) != 6 or any(
        not close(a, b, 5.0e-14)
        for a, b in zip(geometry, expected_geometry, strict=True)
    ):
        raise AssertionError(f"{path}: checkpoint plane geometry is wrong")
    index += 1
    if parse_ints(lines[index].split(), "checkpoint circle flag") != (0,):
        raise AssertionError(f"{path}: checkpoint circle flag is wrong")
    index += 1
    if lines[index].strip() != "slip_wall" or lines[index + 1].strip() != "adiabatic":
        raise AssertionError(f"{path}: checkpoint wall identity is wrong")
    wall_values = parse_floats(lines[index + 2].split(), "checkpoint wall values")
    if len(wall_values) != 4 or any(
        not close(a, b, 5.0e-13)
        for a, b in zip(wall_values, (1300.0, 0.0, 0.0, 0.0), strict=True)
    ):
        raise AssertionError(f"{path}: checkpoint wall values are wrong")
    index += 3
    expected_strings = ("selected", "hllc", "pcm", "mc", "linear")
    if tuple(lines[index + offset].strip() for offset in range(5)) != expected_strings:
        raise AssertionError(f"{path}: checkpoint numerical policy is wrong")
    index += 5
    flags = parse_ints(lines[index].split(), "checkpoint policy flags")
    if flags != (1, 0, 1, 1, 1, 1, 2):
        raise AssertionError(f"{path}: checkpoint policy flags {flags} are wrong")
    index += 1
    numerics = parse_floats(lines[index].split(), "checkpoint numerics")
    if len(numerics) != 5 or any(
        not close(actual, expected, 5.0e-13)
        for actual, expected in zip(
            numerics, (0.01, 0.35, 2.0e-7, 1.0e-12, 0.5), strict=True
        )
    ):
        raise AssertionError(f"{path}: checkpoint numerics are wrong")
    index += 1
    controls = parse_ints(lines[index].split(), "checkpoint regrid controls")
    if controls != (1, 1, 1, 0, 4, 4, 0):
        raise AssertionError(f"{path}: checkpoint regrid controls {controls} are wrong")
    index += 1
    regrid_numerics = parse_floats(lines[index].split(), "checkpoint regrid numerics")
    if len(regrid_numerics) != 3 or any(
        not close(actual, expected, 5.0e-13)
        for actual, expected in zip(regrid_numerics, (1.0e-4, 0.1, 1.0), strict=True)
    ):
        raise AssertionError(f"{path}: checkpoint regrid numerics are wrong")
    index += 1
    root_patch = parse_ints(lines[index].split(), "checkpoint root patch")
    if root_patch != ROOT_PATCH or root_patch == ROOT_SEED:
        raise AssertionError(f"{path}: checkpoint root topology {root_patch} is wrong")
    index += 1
    finest_patch = parse_ints(lines[index].split(), "checkpoint finest patch")
    if finest_patch != FINEST_PATCH or finest_patch == FINEST_SEED:
        raise AssertionError(f"{path}: checkpoint finest topology {finest_patch} is wrong")
    index += 1
    metadata = lines[index].split()
    if len(metadata) != 5:
        raise AssertionError(f"{path}: checkpoint lifecycle metadata is malformed")
    checkpoint_time, minimum_dt = parse_floats(
        metadata[:2], "checkpoint lifecycle metadata"
    )
    try:
        steps = int(metadata[2])
        regrids = int(metadata[3])
        base_density = float(metadata[4])
    except ValueError as error:
        raise AssertionError(f"{path}: checkpoint lifecycle metadata is malformed") from error
    if (
        checkpoint_time <= 0.0
        or checkpoint_time >= FINAL_TIME
        or minimum_dt <= 0.0
        or minimum_dt > checkpoint_time
        or steps != 1
        or regrids != 1
        or not math.isfinite(base_density)
        or base_density <= 0.0
    ):
        raise AssertionError(f"{path}: checkpoint lifecycle metadata is invalid")
    index += 1

    index = checkpoint_field(lines, index, 12, 12, "checkpoint root field")
    index = checkpoint_field(lines, index, 12, 16, "checkpoint middle field")
    index = checkpoint_field(lines, index, 16, 24, "checkpoint finest field")
    if index != len(lines) - 1:
        raise AssertionError(f"{path}: unexpected records after checkpoint fields")
    return data, checkpoint_time, steps, regrids


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    for level in LEVEL_CONTRACTS:
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
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    datasets: dict[str, tuple[bytes, list[dict[str, str]], float, float]] = {}
    for level in LEVEL_CONTRACTS:
        for stage in ("fixed", "reference", "stopped", "restarted"):
            datasets[f"{stage}_{level}"] = check_csv(
                getattr(arguments, f"{stage}_{level}"), level
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
            getattr(arguments, f"expected_reference_{level}_sha256"),
            f"reference {level} SHA-256",
        )
        check_hash(
            stopped,
            getattr(arguments, f"expected_stopped_{level}_sha256"),
            f"stopped {level} SHA-256",
        )

    reference_time = datasets["reference_root"][2]
    stopped_time = datasets["stopped_root"][2]
    restarted_time = datasets["restarted_root"][2]
    if not close(reference_time, FINAL_TIME, 5.0e-14):
        raise AssertionError("selected reference output did not reach final time")
    if not close(restarted_time, FINAL_TIME, 5.0e-14):
        raise AssertionError("selected restart output did not reach final time")
    if not close(datasets["fixed_root"][2], FINAL_TIME, 5.0e-14):
        raise AssertionError("fixed reference output did not reach final time")
    if not 0.0 < stopped_time < FINAL_TIME:
        raise AssertionError("selected checkpoint-stop output did not stop early")
    for stage in ("fixed", "reference", "stopped", "restarted"):
        times = [datasets[f"{stage}_{level}"][2] for level in LEVEL_CONTRACTS]
        if any(not close(value, times[0], 5.0e-14) for value in times[1:]):
            raise AssertionError(f"{stage}: three output levels have inconsistent times")

    active_species = max(
        max(float(row["Y_HO2"]), float(row["Y_H2O2"]))
        for level in LEVEL_CONTRACTS
        for row in datasets[f"reference_{level}"][1]
    )
    if active_species <= 1.0e-15:
        raise AssertionError("selected dynamic three-level chemistry is inactive")

    reference_log_time, reference_steps, reference_regrids, reference_conservation = check_log(
        arguments.reference_log, False
    )
    stopped_log_time, stopped_steps, stopped_regrids, _ = check_log(
        arguments.stopped_log, True
    )
    restarted_log_time, restarted_steps, restarted_regrids, restarted_conservation = check_log(
        arguments.restarted_log, False
    )
    if not close(reference_log_time, reference_time, 5.0e-14):
        raise AssertionError("selected reference log/output time mismatch")
    if not close(stopped_log_time, stopped_time, 5.0e-14):
        raise AssertionError("selected checkpoint-stop log/output time mismatch")
    if not close(restarted_log_time, restarted_time, 5.0e-14):
        raise AssertionError("selected restart log/output time mismatch")
    if stopped_steps >= reference_steps or restarted_steps != reference_steps:
        raise AssertionError("selected restart step counts are inconsistent")
    if (reference_regrids, stopped_regrids, restarted_regrids) != (1, 1, 1):
        raise AssertionError("selected restart regrid counts are inconsistent")
    if reference_conservation != restarted_conservation:
        raise AssertionError(
            "selected cumulative conservation diagnostic is not textually identical: "
            f"{reference_conservation} != {restarted_conservation}"
        )

    checkpoint_data, checkpoint_time, checkpoint_steps, checkpoint_regrids = (
        check_checkpoint(arguments.checkpoint)
    )
    if not close(checkpoint_time, stopped_time, 5.0e-14):
        raise AssertionError("checkpoint/output stop time mismatch")
    if checkpoint_steps != stopped_steps or checkpoint_regrids != stopped_regrids:
        raise AssertionError("checkpoint/log lifecycle mismatch")
    check_hash(
        checkpoint_data,
        arguments.expected_checkpoint_sha256,
        "selected dynamic three-level checkpoint SHA-256",
    )

    maximum_closure_error = max(dataset[3] for dataset in datasets.values())
    print("root_rows=144")
    print("middle_rows=192")
    print("finest_rows=384")
    print(f"checkpoint_time={checkpoint_time:.16e}")
    print(f"final_time={FINAL_TIME:.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    print(f"active_ho2_or_h2o2={active_species:.16e}")
    for level in LEVEL_CONTRACTS:
        print(f"reference_{level}_sha256={digest(datasets[f'reference_{level}'][0])}")
    print(f"checkpoint_sha256={digest(checkpoint_data)}")
    print("selected_dynamic_three_level_restart_exact_match=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
