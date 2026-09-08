#!/usr/bin/env python3
"""Validate selected two-level multipatch reactive EB AMR restart parity."""

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


MAGIC = "PELEF_REACTIVE_EB_AMR_PATCH_SET_2D_CHECKPOINT"
SCHEMA = 4
FINAL_TIME = 1.0e-8
PATCHES = ((2, 6, 6, 10, 2), (9, 13, 9, 13, 2))

LEVEL_CONTRACTS = {
    "root": {
        "rows": 196,
        "counts": (130, 21, 45),
        "shape": (14, 14),
        "domain": (0.0, 0.01, 0.0, 0.01),
    },
    "patch1": {
        "rows": 100,
        "counts": (45, 19, 36),
        "shape": (10, 10),
        "domain": (
            0.0007142857142857143,
            0.004285714285714286,
            0.0035714285714285713,
            0.007142857142857143,
        ),
    },
    "patch2": {
        "rows": 100,
        "counts": (100, 0, 0),
        "shape": (10, 10),
        "domain": (
            0.005714285714285714,
            0.009285714285714286,
            0.005714285714285714,
            0.009285714285714286,
        ),
    },
}


def check_csv(
    path: Path, level: str
) -> tuple[bytes, list[dict[str, str]], float, float]:
    contract = LEVEL_CONTRACTS[level]
    data = path.read_bytes()
    expected_header = [*BASE_COLUMNS, *(f"Y_{name}" for name in SPECIES)]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_header:
            raise AssertionError(f"{path}: selected multipatch CSV schema mismatch")
        rows = list(reader)
    if len(rows) != contract["rows"]:
        raise AssertionError(
            f"{path}: expected {contract['rows']} rows, got {len(rows)}"
        )
    if any(
        None in row or any(value is None for value in row.values()) for row in rows
    ):
        raise AssertionError(f"{path}: incomplete CSV row")

    nx, ny = contract["shape"]
    x_lower, x_upper, y_lower, y_upper = contract["domain"]
    dx = (x_upper - x_lower) / nx
    dy = (y_upper - y_lower) / ny
    cell_counts = {0: 0, 1: 0, 2: 0}
    maximum_closure_error = 0.0
    for index, row in enumerate(rows):
        parse_floats([row[column] for column in expected_header], str(path))
        expected_x = x_lower + (index % nx + 0.5) * dx
        expected_y = y_lower + (index // nx + 0.5) * dy
        if not close(float(row["x"]), expected_x, 5.0e-14):
            raise AssertionError(f"{path}: x-coordinate mismatch at row {index}")
        if not close(float(row["y"]), expected_y, 5.0e-14):
            raise AssertionError(f"{path}: y-coordinate mismatch at row {index}")
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

    expected_counts = {
        2: contract["counts"][0],
        1: contract["counts"][1],
        0: contract["counts"][2],
    }
    if cell_counts != expected_counts:
        raise AssertionError(f"{path}: cell counts {cell_counts} != {expected_counts}")
    times = {float(row["time"]) for row in rows}
    if len(times) != 1:
        raise AssertionError(f"{path}: multiple output times")
    return data, rows, next(iter(times)), maximum_closure_error


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
            raise AssertionError(f"{label}: field width is wrong")
        if values[0] <= 0.0 or values[4] <= 0.0 or values[-1] <= 0.0:
            raise AssertionError(f"{label}: nonpositive state at row {row_index}")
        species_densities = values[5:NVAR]
        tolerance = 5.0e-12 * values[0]
        if any(value < -tolerance for value in species_densities):
            raise AssertionError(f"{label}: negative species density")
        if abs(sum(species_densities) - values[0]) > tolerance:
            raise AssertionError(f"{label}: species closure mismatch")
        index += 1
    return index


def check_checkpoint(path: Path) -> tuple[bytes, float, int, int]:
    data = path.read_bytes()
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise AssertionError(f"{path}: checkpoint is not UTF-8 text") from error
    if len(lines) < 400 or lines[0].strip() != MAGIC:
        raise AssertionError(f"{path}: selected multipatch magic is missing")
    if lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError(f"{path}: incomplete selected multipatch checkpoint")
    header = parse_ints(lines[1].split(), "checkpoint header")
    if header != (SCHEMA, len(SPECIES), NVAR, len(PATCHES)):
        raise AssertionError(f"{path}: checkpoint header {header} is wrong")
    if lines[2].strip() != "SELECTED_CONTEXT":
        raise AssertionError(f"{path}: selected context marker is missing")
    if lines[3].strip() != BUNDLE_SHA256 or lines[4].strip() != INTEGRATOR:
        raise AssertionError(f"{path}: selected mechanism identity is wrong")
    if parse_ints(lines[5].split(), "composition count") != (len(SPECIES),):
        raise AssertionError(f"{path}: selected composition count is wrong")
    composition = parse_floats(lines[6].split(), "checkpoint composition")
    if len(composition) != len(SPECIES) or any(
        abs(actual - expected) > 5.0e-14
        for actual, expected in zip(
            composition, EXPECTED_COMPOSITION, strict=True
        )
    ):
        raise AssertionError(f"{path}: selected composition is wrong")
    if any(value < 0.0 for value in composition) or abs(sum(composition) - 1.0) > 5.0e-12:
        raise AssertionError(f"{path}: selected composition is not physical")
    if lines[7].strip() != "COMPOSITE_BASELINE":
        raise AssertionError(f"{path}: composite baseline marker is missing")
    if parse_ints(lines[8].split(), "baseline count") != (NVAR,):
        raise AssertionError(f"{path}: composite baseline count is wrong")
    baseline = parse_floats(lines[9].split(), "checkpoint baseline")
    if len(baseline) != NVAR or baseline[0] <= 0.0 or baseline[4] <= 0.0:
        raise AssertionError(f"{path}: composite baseline is invalid")
    if any(value < 0.0 for value in baseline[5:]):
        raise AssertionError(f"{path}: composite baseline species is negative")
    if abs(sum(baseline[5:]) - baseline[0]) > 5.0e-12 * baseline[0]:
        raise AssertionError(f"{path}: composite baseline closure mismatch")

    species_start = 10
    species_end = species_start + len(SPECIES)
    if tuple(line.strip() for line in lines[species_start:species_end]) != SPECIES:
        raise AssertionError(f"{path}: checkpoint species order is wrong")
    index = species_end
    if lines[index].strip() != "plane":
        raise AssertionError(f"{path}: checkpoint geometry is wrong")
    index += 1
    if parse_ints(lines[index].split(), "root dimensions") != (14, 14, 2):
        raise AssertionError(f"{path}: root dimensions are wrong")
    index += 1
    domain = parse_floats(lines[index].split(), "checkpoint domain")
    if any(
        not close(actual, expected, 5.0e-14)
        for actual, expected in zip(
            domain, (0.0, 0.01, 0.0, 0.01), strict=True
        )
    ):
        raise AssertionError(f"{path}: checkpoint domain is wrong")
    index += 1
    geometry = parse_floats(lines[index].split(), "checkpoint geometry")
    if any(
        not close(actual, expected, 5.0e-14)
        for actual, expected in zip(
            geometry, (1.0, 1.0, 0.0078, 0.0, 0.0, 1.0), strict=True
        )
    ):
        raise AssertionError(f"{path}: checkpoint plane geometry is wrong")
    index += 1
    if parse_ints(lines[index].split(), "circle flag") != (0,):
        raise AssertionError(f"{path}: checkpoint circle flag is wrong")
    index += 1
    if lines[index].strip() != "slip_wall" or lines[index + 1].strip() != "adiabatic":
        raise AssertionError(f"{path}: checkpoint wall identity is wrong")
    wall = parse_floats(lines[index + 2].split(), "checkpoint wall values")
    if wall != (1300.0, 0.0, 0.0, 0.0):
        raise AssertionError(f"{path}: checkpoint wall values are wrong")
    index += 3
    expected_strings = ("selected", "hllc", "characteristic_plm", "mc", "linear")
    if tuple(lines[index + offset].strip() for offset in range(5)) != expected_strings:
        raise AssertionError(f"{path}: checkpoint numerical policy is wrong")
    index += 5
    flags = parse_ints(lines[index].split(), "checkpoint flags")
    if flags != (1, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 0, 5, 5, 0, 1):
        raise AssertionError(f"{path}: checkpoint flags {flags} are wrong")
    index += 1
    numerics = parse_floats(lines[index].split(), "checkpoint numerics")
    expected_numerics = (0.01, 0.35, 2.0e-6, 1.0e-12, 0.5, 1.0e-4, 0.1, 1.0)
    if any(
        not close(actual, expected, 5.0e-13)
        for actual, expected in zip(numerics, expected_numerics, strict=True)
    ):
        raise AssertionError(f"{path}: checkpoint numerics are wrong")
    index += 1
    metadata = lines[index].split()
    if len(metadata) != 5:
        raise AssertionError(f"{path}: lifecycle metadata is malformed")
    checkpoint_time, minimum_dt = parse_floats(metadata[:2], "lifecycle metadata")
    try:
        steps = int(metadata[2])
        regrids = int(metadata[3])
        base_density = float(metadata[4])
    except ValueError as error:
        raise AssertionError(f"{path}: lifecycle metadata is malformed") from error
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
        raise AssertionError(f"{path}: lifecycle metadata is invalid")
    index += 1
    index = checkpoint_field(lines, index, 14, 14, "checkpoint root")
    for child, expected_patch in enumerate(PATCHES, start=1):
        patch = parse_ints(lines[index].split(), f"checkpoint patch {child}")
        if patch != expected_patch:
            raise AssertionError(f"{path}: patch {child} topology is wrong")
        index += 1
        index = checkpoint_field(lines, index, 10, 10, f"checkpoint patch {child}")
    if index != len(lines) - 1:
        raise AssertionError(f"{path}: unexpected record after checkpoint fields")
    return data, checkpoint_time, steps, regrids


def check_log(path: Path, stopped: bool) -> tuple[float, int, int, str]:
    text = path.read_text(encoding="utf-8")
    expected_flag = "T" if stopped else "F"
    required = (
        rf"^Bundle SHA-256:\s+{BUNDLE_SHA256}\s*$",
        rf"^Chemistry integrator:\s+{INTEGRATOR}\s*$",
        rf"^Species:\s+{len(SPECIES)}\s*$",
        r"^Reactions:\s+29\s*$",
        r"^Fine patches:\s+2\s*$",
        r"^Fine patch 1 coarse bounds:\s+2 6 6 10\s*$",
        r"^Fine patch 2 coarse bounds:\s+9 13 9 13\s*$",
        r"^Completed regrids:\s+1\s*$",
        r"^Chemistry:\s+T\s*$",
        rf"^Stopped after checkpoint:\s+{expected_flag}\s*$",
    )
    for pattern in required:
        if re.search(pattern, text, re.MULTILINE) is None:
            raise AssertionError(f"{path}: missing log identity {pattern}")
    steps_match = re.search(r"^Completed coarse steps:\s+([0-9]+)\s*$", text, re.MULTILINE)
    time_match = re.search(r"^Final time:\s+([0-9.Ee+-]+)\s*$", text, re.MULTILINE)
    conservation_match = re.search(
        r"^Maximum composite conservation error:\s+([0-9.Ee+-]+)\s*$",
        text,
        re.MULTILINE,
    )
    if steps_match is None or time_match is None or conservation_match is None:
        raise AssertionError(f"{path}: missing lifecycle diagnostic")
    steps = int(steps_match.group(1))
    final_time = float(time_match.group(1))
    conservation_text = conservation_match.group(1)
    if steps < 1 or not math.isfinite(final_time):
        raise AssertionError(f"{path}: invalid lifecycle diagnostic")
    return final_time, steps, 1, conservation_text


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
        require_exact(fixed, reference, f"fixed/selected {level} mismatch")
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

    for stage in ("fixed", "reference", "restarted"):
        times = [datasets[f"{stage}_{level}"][2] for level in LEVEL_CONTRACTS]
        if any(not close(value, FINAL_TIME, 5.0e-14) for value in times):
            raise AssertionError(f"{stage}: output did not reach final time")
    stopped_times = [datasets[f"stopped_{level}"][2] for level in LEVEL_CONTRACTS]
    if not all(0.0 < value < FINAL_TIME for value in stopped_times):
        raise AssertionError("checkpoint-stop output did not stop early")
    if any(not close(value, stopped_times[0], 5.0e-14) for value in stopped_times[1:]):
        raise AssertionError("checkpoint-stop levels have inconsistent times")

    active_species = max(
        max(float(row["Y_HO2"]), float(row["Y_H2O2"]))
        for level in LEVEL_CONTRACTS
        for row in datasets[f"reference_{level}"][1]
    )
    if active_species <= 1.0e-15:
        raise AssertionError("selected multipatch chemistry is inactive")

    reference_log = check_log(arguments.reference_log, False)
    stopped_log = check_log(arguments.stopped_log, True)
    restarted_log = check_log(arguments.restarted_log, False)
    if not close(reference_log[0], FINAL_TIME, 5.0e-14):
        raise AssertionError("reference log time mismatch")
    if not close(stopped_log[0], stopped_times[0], 5.0e-14):
        raise AssertionError("checkpoint-stop log time mismatch")
    if not close(restarted_log[0], FINAL_TIME, 5.0e-14):
        raise AssertionError("restart log time mismatch")
    if stopped_log[1] >= reference_log[1] or restarted_log[1] != reference_log[1]:
        raise AssertionError("restart step counts are inconsistent")
    if reference_log[3] != restarted_log[3]:
        raise AssertionError("cumulative conservation diagnostic changed on restart")

    checkpoint_data, checkpoint_time, checkpoint_steps, checkpoint_regrids = (
        check_checkpoint(arguments.checkpoint)
    )
    if not close(checkpoint_time, stopped_times[0], 5.0e-14):
        raise AssertionError("checkpoint/output stop time mismatch")
    if checkpoint_steps != stopped_log[1] or checkpoint_regrids != 1:
        raise AssertionError("checkpoint/log lifecycle mismatch")
    check_hash(
        checkpoint_data,
        arguments.expected_checkpoint_sha256,
        "selected multipatch checkpoint SHA-256",
    )

    maximum_closure_error = max(dataset[3] for dataset in datasets.values())
    print("root_rows=196")
    print("patch1_rows=100")
    print("patch2_rows=100")
    print(f"checkpoint_time={checkpoint_time:.16e}")
    print(f"final_time={FINAL_TIME:.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    print(f"active_ho2_or_h2o2={active_species:.16e}")
    for level in LEVEL_CONTRACTS:
        print(f"reference_{level}_sha256={digest(datasets[f'reference_{level}'][0])}")
    print(f"checkpoint_sha256={digest(checkpoint_data)}")
    print("selected_multipatch_restart_exact_match=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
