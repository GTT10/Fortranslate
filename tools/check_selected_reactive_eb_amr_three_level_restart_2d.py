#!/usr/bin/env python3
"""Validate selected static three-level reactive EB AMR restart parity.

The selected three-level envelope deliberately has its own schema (4).  This
checker keeps the fixed schema-3 output as a byte reference, checks the
selected context and composite baseline in the checkpoint, and validates the
three output levels for physicality and species closure.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re


MAGIC = "PELEF_REACTIVE_EB_AMR_THREE_LEVEL_2D_CHECKPOINT"
SCHEMA = 4
BUNDLE_SHA256 = (
    "f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3"
)
INTEGRATOR = "implicit"
SPECIES = (
    "H2",
    "H",
    "O",
    "O2",
    "OH",
    "H2O",
    "HO2",
    "H2O2",
    "AR",
    "N2",
)
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


LEVEL_CONTRACTS = {
    "root": {
        "rows": 64,
        "counts": (32, 8, 24),
        "shape": (8, 8),
        "domain": (0.0, 0.008, 0.0, 0.008),
    },
    "middle": {
        "rows": 144,
        "counts": (84, 12, 48),
        "shape": (12, 12),
        "domain": (0.001, 0.007, 0.001, 0.007),
    },
    "finest": {
        "rows": 256,
        "counts": (160, 16, 80),
        "shape": (16, 16),
        "domain": (0.002, 0.006, 0.002, 0.006),
    },
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def close(actual: float, expected: float, tolerance: float = 2.0e-14) -> bool:
    return abs(actual - expected) <= tolerance


def parse_floats(values: list[str], label: str) -> tuple[float, ...]:
    try:
        result = tuple(float(value) for value in values)
    except ValueError as error:
        raise AssertionError(f"{label}: nonnumeric record") from error
    if not all(math.isfinite(value) for value in result):
        raise AssertionError(f"{label}: nonfinite record")
    return result


def parse_ints(values: list[str], label: str) -> tuple[int, ...]:
    try:
        return tuple(int(value) for value in values)
    except ValueError as error:
        raise AssertionError(f"{label}: noninteger record") from error


def read_csv(path: Path, level: str) -> tuple[bytes, list[dict[str, str]]]:
    contract = LEVEL_CONTRACTS[level]
    data = path.read_bytes()
    expected_header = [*BASE_COLUMNS, *(f"Y_{name}" for name in SPECIES)]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_header:
            raise AssertionError(f"{path}: selected three-level CSV schema mismatch")
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
        values = parse_floats([row[column] for column in expected_header], str(path))
        if not all(math.isfinite(value) for value in values):
            raise AssertionError(f"{path}: nonfinite CSV value")
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

    expected_counts = {2: contract["counts"][0], 1: contract["counts"][1], 0: contract["counts"][2]}
    if cell_counts != expected_counts:
        raise AssertionError(
            f"{path}: cell counts {cell_counts} != {expected_counts}"
        )
    times = {float(row["time"]) for row in rows}
    if len(times) != 1:
        raise AssertionError(f"{path}: levels contain multiple output times")
    return data, rows, next(iter(times)), maximum_closure_error


def require_exact(left: bytes, right: bytes, label: str) -> None:
    if left != right:
        raise AssertionError(
            f"{label}: left_sha256={digest(left)}, right_sha256={digest(right)}"
        )


def check_hash(data: bytes, expected: str, label: str) -> None:
    actual = digest(data)
    if not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise AssertionError(f"{label}: expected hash must be 64 lowercase hex digits")
    if actual != expected:
        raise AssertionError(f"{label}: {actual} != {expected}")


def check_log(path: Path, stopped: bool) -> tuple[float, int, str]:
    text = path.read_text(encoding="utf-8")
    expected_flag = "T" if stopped else "F"
    required = (
        rf"^Bundle SHA-256:\s+{BUNDLE_SHA256}\s*$",
        rf"^Chemistry integrator:\s+{INTEGRATOR}\s*$",
        rf"^Species:\s+{len(SPECIES)}\s*$",
        r"^Reactions:\s+29\s*$",
        r"^Completed regrids:\s+0\s*$",
        r"^Root patch bounds:\s+2 7 2 7\s*$",
        r"^Middle patch bounds:\s+3 10 3 10\s*$",
        r"^Refinement ratio:\s+2\s*$",
        rf"^Stopped after checkpoint:\s+{expected_flag}\s*$",
    )
    for pattern in required:
        if re.search(pattern, text, re.MULTILINE) is None:
            raise AssertionError(f"{path}: missing log identity {pattern}")

    steps_match = re.search(r"^Completed coarse steps:\s+([0-9]+)\s*$", text, re.MULTILINE)
    if steps_match is None or int(steps_match.group(1)) < 1:
        raise AssertionError(f"{path}: invalid completed-step count")
    steps = int(steps_match.group(1))
    time_match = re.search(r"^Final time:\s+([0-9.Ee+-]+)\s*$", text, re.MULTILINE)
    conservation_match = re.search(
        r"^Maximum composite conservation error:\s+([0-9.Ee+-]+)\s*$",
        text,
        re.MULTILINE,
    )
    if time_match is None or conservation_match is None:
        raise AssertionError(f"{path}: missing final-time or conservation diagnostic")
    final_time = float(time_match.group(1))
    conservation_text = conservation_match.group(1)
    conservation_value = float(conservation_text)
    if not math.isfinite(final_time):
        raise AssertionError(f"{path}: nonfinite final time")
    if not math.isfinite(conservation_value) or conservation_value < 0.0:
        raise AssertionError(f"{path}: invalid conservation diagnostic")
    return final_time, steps, conservation_text


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
        if any(not math.isfinite(value) for value in values):
            raise AssertionError(f"{label} row {row_index}: nonfinite state")
        if values[0] <= 0.0 or values[4] <= 0.0 or values[-1] <= 0.0:
            raise AssertionError(f"{label} row {row_index}: nonpositive state")
        species_densities = values[5:NVAR]
        if any(value < 0.0 for value in species_densities):
            raise AssertionError(f"{label} row {row_index}: negative species density")
        closure_error = abs(sum(species_densities) - values[0]) / values[0]
        if closure_error > 5.0e-12:
            raise AssertionError(
                f"{label} row {row_index}: species closure error {closure_error}"
            )
        index += 1
    return index


def check_checkpoint(path: Path) -> tuple[bytes, float, int]:
    data = path.read_bytes()
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise AssertionError(f"{path}: checkpoint is not UTF-8 text") from error
    if len(lines) < 505 or lines[0].strip() != MAGIC:
        raise AssertionError(f"{path}: selected three-level checkpoint magic is missing")
    if lines[-1].strip() != "END_CHECKPOINT":
        raise AssertionError(f"{path}: incomplete selected three-level checkpoint")

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
    if len(composition) != len(SPECIES):
        raise AssertionError(f"{path}: selected composition has wrong size")
    if any(
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
    if len(baseline) != NVAR:
        raise AssertionError(f"{path}: composite baseline has wrong size")
    if any(not math.isfinite(value) for value in baseline):
        raise AssertionError(f"{path}: composite baseline is nonfinite")
    if baseline[0] <= 0.0 or baseline[4] <= 0.0:
        raise AssertionError(f"{path}: composite baseline rho/energy is not positive")
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
    if parse_ints(lines[index].split(), "checkpoint coarse dimensions") != (8, 8, 2):
        raise AssertionError(f"{path}: checkpoint coarse dimensions are wrong")
    index += 1
    domain = parse_floats(lines[index].split(), "checkpoint domain")
    if len(domain) != 4 or any(not close(a, b, 5.0e-14) for a, b in zip(domain, (0.0, 0.008, 0.0, 0.008), strict=True)):
        raise AssertionError(f"{path}: checkpoint domain is wrong")
    index += 1
    geometry = parse_floats(lines[index].split(), "checkpoint plane geometry")
    expected_geometry = (1.0, 0.0, 0.00337, 0.0, 0.0, 1.0)
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
    if len(wall_values) != 4 or any(not close(a, b, 5.0e-13) for a, b in zip(wall_values, (1300.0, 0.0, 0.0, 0.0), strict=True)):
        raise AssertionError(f"{path}: checkpoint wall values are wrong")
    index += 3
    if lines[index].strip() != "selected":
        raise AssertionError(f"{path}: checkpoint chemistry model is wrong")
    if lines[index + 1].strip() != "hllc":
        raise AssertionError(f"{path}: checkpoint Riemann solver is wrong")
    if lines[index + 2].strip() != "pcm" or lines[index + 3].strip() != "mc":
        raise AssertionError(f"{path}: checkpoint reconstruction policy is wrong")
    if lines[index + 4].strip() != "linear":
        raise AssertionError(f"{path}: checkpoint prolongation policy is wrong")
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
    if parse_ints(lines[index].split(), "checkpoint root patch") != (2, 7, 2, 7, 2):
        raise AssertionError(f"{path}: checkpoint root patch is wrong")
    index += 1
    if parse_ints(lines[index].split(), "checkpoint middle patch") != (3, 10, 3, 10, 2):
        raise AssertionError(f"{path}: checkpoint middle patch is wrong")
    index += 1
    metadata = lines[index].split()
    if len(metadata) != 4:
        raise AssertionError(f"{path}: checkpoint time metadata is malformed")
    metadata_values = parse_floats(metadata[:2], "checkpoint time metadata")
    try:
        steps = int(metadata[2])
    except ValueError as error:
        raise AssertionError(f"{path}: checkpoint step metadata is malformed") from error
    base_density = float(metadata[3])
    checkpoint_time, minimum_dt = metadata_values
    if (
        checkpoint_time <= 0.0
        or checkpoint_time >= FINAL_TIME
        or minimum_dt <= 0.0
        or minimum_dt > checkpoint_time
        or steps < 1
        or base_density <= 0.0
    ):
        raise AssertionError(f"{path}: checkpoint time metadata is invalid")
    index += 1

    index = checkpoint_field(lines, index, 8, 8, "checkpoint root field")
    index = checkpoint_field(lines, index, 12, 12, "checkpoint middle field")
    index = checkpoint_field(lines, index, 16, 16, "checkpoint finest field")
    if index != len(lines) - 1:
        raise AssertionError(f"{path}: unexpected records after checkpoint fields")
    return data, checkpoint_time, steps


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
    if not close(reference_time, FINAL_TIME, 5.0e-14) or not close(restarted_time, FINAL_TIME, 5.0e-14):
        raise AssertionError("selected reference/restart output did not reach final time")
    if not 0.0 < stopped_time < FINAL_TIME:
        raise AssertionError("selected checkpoint-stop output did not stop before final time")
    for stage in ("fixed", "reference", "stopped", "restarted"):
        times = [datasets[f"{stage}_{level}"][2] for level in LEVEL_CONTRACTS]
        if any(not close(value, times[0], 5.0e-14) for value in times[1:]):
            raise AssertionError(f"{stage}: three output levels have inconsistent times")
    if not close(datasets["fixed_root"][2], FINAL_TIME, 5.0e-14):
        raise AssertionError("fixed reference output did not reach final time")

    reference_log_time, reference_steps, reference_conservation = check_log(
        arguments.reference_log, False
    )
    stopped_log_time, stopped_steps, _ = check_log(arguments.stopped_log, True)
    restarted_log_time, restarted_steps, restarted_conservation = check_log(
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
    if reference_conservation != restarted_conservation:
        raise AssertionError(
            "selected cumulative conservation diagnostic is not textually identical: "
            f"{reference_conservation} != {restarted_conservation}"
        )

    checkpoint_data, checkpoint_time, checkpoint_steps = check_checkpoint(
        arguments.checkpoint
    )
    if not close(checkpoint_time, stopped_time, 5.0e-14):
        raise AssertionError("checkpoint/output stop time mismatch")
    if checkpoint_steps != stopped_steps:
        raise AssertionError("checkpoint/log stop step mismatch")
    check_hash(
        checkpoint_data,
        arguments.expected_checkpoint_sha256,
        "selected three-level checkpoint SHA-256",
    )

    maximum_closure_error = max(dataset[3] for dataset in datasets.values())
    print("root_rows=64")
    print("middle_rows=144")
    print("finest_rows=256")
    print(f"checkpoint_time={checkpoint_time:.16e}")
    print(f"final_time={FINAL_TIME:.16e}")
    print(f"maximum_closure_error={maximum_closure_error:.16e}")
    for level in LEVEL_CONTRACTS:
        print(f"reference_{level}_sha256={digest(datasets[f'reference_{level}'][0])}")
    print(f"checkpoint_sha256={digest(checkpoint_data)}")
    print("selected_reactive_eb_amr_three_level_restart_exact_match=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
