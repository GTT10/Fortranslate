#!/usr/bin/env python3
"""Require exact serial/MPI and changed-rank restart parity for 3D AMR."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
from typing import Optional


CHECKPOINT_MAGIC = "PELEF_AMR_REACTIVE_3D_CHECKPOINT"
SELECTED_CONTEXT_MARKER = "SELECTED_CONTEXT"
TRANSPORT_CONTEXT_MARKER = "TRANSPORT_CONTEXT"
TRANSPORT_OPERATOR = "STATIC_AMR_3D_R_T_H_T_R_V1"
TRANSPORT_CONVENTION = (
    "EPSILON_OVER_K_K;SIGMA_ANGSTROM;DIPOLE_DEBYE;"
    "POLARIZABILITY_ANGSTROM3;ROT_RELAX_DIMENSIONLESS"
)
TRANSPORT_PHASE = "POST_ACCEPTED_COARSE_STEP"
TRANSPORT_LOG_FIELDS = (
    "Completed coarse steps",
    "Completed hydro fine substeps",
    "Final time",
    "Maximum reflux correction",
    "Completed transport fine substeps",
    "Maximum transport diffusivity",
    "Minimum transport theta",
    "Composite conservation error",
    "Maximum elemental conservation error",
    "Maximum species integral change",
    "Average-down synchronization error",
    "Maximum species closure error",
)
HEX_DIGITS = frozenset("0123456789abcdef")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def level_path(prefix: Path, level: str) -> Path:
    return Path(f"{prefix}_{level}.csv")


def require_equal(reference: Path, candidate: Path, label: str) -> str:
    reference_data = reference.read_bytes()
    candidate_data = candidate.read_bytes()
    reference_digest = digest(reference_data)
    candidate_digest = digest(candidate_data)
    if candidate_data != reference_data:
        raise SystemExit(
            f"{label} mismatch: reference_sha256={reference_digest}, "
            f"candidate_sha256={candidate_digest}"
        )
    return reference_digest


def transport_log_diagnostics(path: Path) -> dict[str, str]:
    diagnostics: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition(":")
        if separator and key in TRANSPORT_LOG_FIELDS:
            if key in diagnostics:
                raise SystemExit(f"{path}: duplicate diagnostic {key!r}")
            diagnostics[key] = value.strip()
    missing = [field for field in TRANSPORT_LOG_FIELDS if field not in diagnostics]
    if missing:
        raise SystemExit(f"{path}: missing transport diagnostics: {missing}")
    return diagnostics


def parse_finite_values(line: str, count: int, label: str) -> list[float]:
    try:
        values = [float(value) for value in line.split()]
    except ValueError as error:
        raise SystemExit(f"{label}: invalid real value") from error
    if len(values) != count or not all(math.isfinite(value) for value in values):
        raise SystemExit(f"{label}: invalid real-vector shape or value")
    return values


def validate_selected_context(
    path: Path,
    lines: list[str],
    start: int,
    species_count: int,
    expected_bundle_sha256: Optional[str],
    expected_integrator: Optional[str],
) -> int:
    cursor = start

    def take(label: str) -> str:
        nonlocal cursor
        if cursor >= len(lines):
            raise SystemExit(f"{path}: truncated selected context at {label}")
        value = lines[cursor].strip()
        cursor += 1
        return value

    if take("marker") != SELECTED_CONTEXT_MARKER:
        raise SystemExit(f"{path}: missing selected checkpoint context")
    bundle_sha256 = take("bundle SHA-256")
    if len(bundle_sha256) != 64 or any(
        character not in HEX_DIGITS for character in bundle_sha256
    ):
        raise SystemExit(f"{path}: invalid selected bundle SHA-256")
    if (
        expected_bundle_sha256 is not None
        and bundle_sha256 != expected_bundle_sha256.lower()
    ):
        raise SystemExit(
            f"{path}: selected bundle SHA-256 mismatch: "
            f"expected={expected_bundle_sha256.lower()}, actual={bundle_sha256}"
        )
    integrator = take("chemistry integrator")
    if integrator not in {"explicit", "implicit"}:
        raise SystemExit(f"{path}: invalid selected chemistry integrator")
    if expected_integrator is not None and integrator != expected_integrator:
        raise SystemExit(
            f"{path}: selected chemistry integrator mismatch: "
            f"expected={expected_integrator}, actual={integrator}"
        )
    try:
        composition_count = int(take("composition size"))
    except ValueError as error:
        raise SystemExit(f"{path}: invalid selected composition size") from error
    if composition_count != species_count:
        raise SystemExit(f"{path}: selected composition size mismatch")
    composition = parse_finite_values(
        take("composition"), composition_count, f"{path}: selected composition"
    )
    if min(composition) < 0.0 or not math.isclose(
        sum(composition), 1.0, rel_tol=5.0e-10, abs_tol=5.0e-10
    ):
        raise SystemExit(f"{path}: invalid selected composition")
    if take("reaction marker") != "REACTIONS":
        raise SystemExit(f"{path}: missing selected reaction marker")
    try:
        reaction_count = int(take("reaction count"))
    except ValueError as error:
        raise SystemExit(f"{path}: invalid selected reaction count") from error
    if reaction_count < 1:
        raise SystemExit(f"{path}: empty selected mechanism")
    for reaction_index in range(reaction_count):
        if not take(f"reaction {reaction_index} equation"):
            raise SystemExit(f"{path}: empty selected reaction equation")
        try:
            header = [
                int(value)
                for value in take(f"reaction {reaction_index} header").split()
            ]
        except ValueError as error:
            raise SystemExit(f"{path}: invalid selected reaction header") from error
        if len(header) != 6:
            raise SystemExit(f"{path}: invalid selected reaction header shape")
        kind, reversible, efficiencies_allocated, efficiency_size, nr, np = header
        if (
            kind not in {1, 2, 3}
            or reversible not in {0, 1}
            or efficiencies_allocated not in {0, 1}
            or nr != species_count
            or np != species_count
            or efficiency_size != efficiencies_allocated * species_count
        ):
            raise SystemExit(f"{path}: invalid selected reaction dimensions")
        reactants = parse_finite_values(
            take(f"reaction {reaction_index} reactants"),
            nr,
            f"{path}: selected reactants",
        )
        products = parse_finite_values(
            take(f"reaction {reaction_index} products"),
            np,
            f"{path}: selected products",
        )
        if min(reactants) < 0.0 or min(products) < 0.0:
            raise SystemExit(f"{path}: negative selected stoichiometry")
        parse_finite_values(
            take(f"reaction {reaction_index} rates"),
            9,
            f"{path}: selected rates",
        )
        try:
            troe_enabled = int(take(f"reaction {reaction_index} Troe flag"))
        except ValueError as error:
            raise SystemExit(f"{path}: invalid selected Troe flag") from error
        if troe_enabled not in {0, 1}:
            raise SystemExit(f"{path}: invalid selected Troe flag")
        parse_finite_values(
            take(f"reaction {reaction_index} Troe data"),
            4,
            f"{path}: selected Troe data",
        )
        if efficiencies_allocated:
            efficiencies = parse_finite_values(
                take(f"reaction {reaction_index} efficiencies"),
                efficiency_size,
                f"{path}: selected efficiencies",
            )
            if min(efficiencies) < 0.0:
                raise SystemExit(f"{path}: negative selected efficiency")
    if take("chemistry controls marker") != "CHEMISTRY_CONTROLS":
        raise SystemExit(f"{path}: missing selected chemistry controls")
    controls = parse_finite_values(
        take("chemistry controls"), 2, f"{path}: selected chemistry controls"
    )
    if min(controls) <= 0.0:
        raise SystemExit(f"{path}: invalid selected chemistry controls")
    if take("selected context end marker") != "END_SELECTED_CONTEXT":
        raise SystemExit(f"{path}: missing selected context end marker")
    return cursor


def validate_transport_context(
    path: Path,
    lines: list[str],
    start: int,
    species_count: int,
) -> tuple[int, tuple[float, float]]:
    cursor = start

    def take(label: str) -> str:
        nonlocal cursor
        if cursor >= len(lines):
            raise SystemExit(f"{path}: truncated transport context at {label}")
        value = lines[cursor].strip()
        cursor += 1
        return value

    if take("marker") != TRANSPORT_CONTEXT_MARKER:
        raise SystemExit(f"{path}: missing transport checkpoint context")
    if take("operator") != TRANSPORT_OPERATOR:
        raise SystemExit(f"{path}: transport operator mismatch")
    if take("parameter convention") != TRANSPORT_CONVENTION:
        raise SystemExit(f"{path}: transport parameter convention mismatch")
    if take("checkpoint phase") != TRANSPORT_PHASE:
        raise SystemExit(f"{path}: transport checkpoint phase mismatch")
    try:
        transport_count = int(take("record count"))
    except ValueError as error:
        raise SystemExit(f"{path}: invalid transport record count") from error
    if transport_count != species_count:
        raise SystemExit(f"{path}: transport record count mismatch")
    names: list[str] = []
    for species_index in range(transport_count):
        name = take(f"record {species_index} name")
        if not name or name in names:
            raise SystemExit(f"{path}: invalid transport species order")
        names.append(name)
        try:
            geometry = int(take(f"record {species_index} geometry"))
        except ValueError as error:
            raise SystemExit(f"{path}: invalid transport geometry") from error
        if geometry not in {0, 1, 2}:
            raise SystemExit(f"{path}: invalid transport geometry")
        values = parse_finite_values(
            take(f"record {species_index} values"),
            5,
            f"{path}: transport record",
        )
        if values[0] <= 0.0 or values[1] <= 0.0 or min(values[2:]) < 0.0:
            raise SystemExit(f"{path}: invalid transport record")
    if take("controls marker") != "TRANSPORT_CONTROLS":
        raise SystemExit(f"{path}: missing transport controls")
    try:
        controls = [int(value) for value in take("controls").split()]
    except ValueError as error:
        raise SystemExit(f"{path}: invalid transport controls") from error
    if len(controls) != 4 or any(value not in {0, 1} for value in controls):
        raise SystemExit(f"{path}: invalid transport controls")
    if controls[3] and not controls[2]:
        raise SystemExit(f"{path}: barodiffusion lacks species diffusion")
    transport_cfl = parse_finite_values(
        take("transport CFL"), 1, f"{path}: transport CFL"
    )[0]
    if transport_cfl <= 0.0:
        raise SystemExit(f"{path}: invalid transport CFL")
    if take("diagnostics marker") != "TRANSPORT_DIAGNOSTICS":
        raise SystemExit(f"{path}: missing transport diagnostics")
    maximum_diffusivity, minimum_theta = parse_finite_values(
        take("diagnostics"), 2, f"{path}: transport diagnostics"
    )
    if maximum_diffusivity < 0.0 or not 0.0 <= minimum_theta <= 1.0:
        raise SystemExit(f"{path}: invalid transport diagnostics")
    if take("context end marker") != "END_TRANSPORT_CONTEXT":
        raise SystemExit(f"{path}: missing transport context end marker")
    return cursor, (maximum_diffusivity, minimum_theta)


def validate_checkpoint(
    path: Path,
    final_time: float,
    reconstruction: str,
    limiter: Optional[str],
    expected_schema: Optional[int],
    expected_bundle_sha256: Optional[str],
    expected_integrator: Optional[str],
) -> None:
    lines = path.read_text().splitlines()
    if len(lines) < 10 or lines[0].strip() != CHECKPOINT_MAGIC:
        raise SystemExit(f"{path}: invalid 3D AMR checkpoint header")
    header = lines[1].split()
    if len(header) != 3:
        raise SystemExit(f"{path}: invalid 3D AMR checkpoint schema line")
    schema, species_count, variable_count = map(int, header)
    if schema not in {2, 3, 4, 5} or not 1 <= species_count <= 32 or (
        variable_count != species_count + 5
    ):
        raise SystemExit(f"{path}: invalid 3D AMR checkpoint dimensions")
    if expected_schema is not None and schema != expected_schema:
        raise SystemExit(
            f"{path}: checkpoint schema mismatch: "
            f"expected={expected_schema}, actual={schema}"
        )
    species_start = 2
    if schema in {3, 5}:
        species_start = validate_selected_context(
            path,
            lines,
            species_start,
            species_count,
            expected_bundle_sha256,
            expected_integrator,
        )
    elif expected_bundle_sha256 is not None or expected_integrator is not None:
        raise SystemExit(f"{path}: selected context expected from fixed checkpoint")
    transport_diagnostics: Optional[tuple[float, float]] = None
    if schema in {4, 5}:
        species_start, transport_diagnostics = validate_transport_context(
            path, lines, species_start, species_count
        )
    thermo_index = species_start + 2 * species_count
    reconstruction_index = thermo_index + 3
    limiter_index = reconstruction_index + 1
    if limiter_index >= len(lines):
        raise SystemExit(f"{path}: truncated 3D AMR checkpoint fingerprint")
    if schema in {3, 5} and lines[thermo_index].strip() != "selected":
        raise SystemExit(f"{path}: invalid selected thermodynamic fingerprint")
    if (
        lines[reconstruction_index].strip() != reconstruction
        or lines[limiter_index].strip() not in {"minmod", "mc"}
        or (limiter is not None and lines[limiter_index].strip() != limiter)
    ):
        raise SystemExit(f"{path}: invalid 3D AMR reconstruction fingerprint")
    flags_index = thermo_index + 7
    if flags_index >= len(lines):
        raise SystemExit(f"{path}: truncated 3D AMR checkpoint flags")
    try:
        flags = [int(value) for value in lines[flags_index].split()]
    except ValueError as error:
        raise SystemExit(f"{path}: invalid 3D AMR checkpoint flags") from error
    expected_transport_flag = 1 if schema in {4, 5} else 0
    if len(flags) != 2 or flags[1] != expected_transport_flag:
        raise SystemExit(f"{path}: invalid 3D AMR transport flag")
    metadata_index = thermo_index + 8
    if metadata_index >= len(lines):
        raise SystemExit(f"{path}: truncated 3D AMR checkpoint")
    metadata = lines[metadata_index].split()
    if len(metadata) != 3:
        raise SystemExit(f"{path}: invalid 3D AMR checkpoint metadata")
    checkpoint_time = float(metadata[0])
    maximum_reflux = float(metadata[1])
    steps = int(metadata[2])
    if not 0.0 < checkpoint_time < final_time or steps < 1:
        raise SystemExit(f"{path}: checkpoint is not at an intermediate step")
    if maximum_reflux < 0.0 or lines[-1].strip() != "END_CHECKPOINT":
        raise SystemExit(f"{path}: invalid 3D AMR checkpoint payload")
    print(
        f"checkpoint_sha256={digest(path.read_bytes())}, "
        f"time={checkpoint_time:.16e}, steps={steps}"
    )
    if transport_diagnostics is not None:
        print(
            "checkpoint_maximum_transport_diffusivity="
            f"{transport_diagnostics[0]:.16e}, "
            f"minimum_transport_theta={transport_diagnostics[1]:.16e}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-prefix", type=Path, required=True)
    parser.add_argument(
        "--candidate-prefix", type=Path, action="append", required=True
    )
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--expected-checkpoint-sha256")
    parser.add_argument("--expected-schema", type=int, choices=(2, 3, 4, 5))
    parser.add_argument("--expected-bundle-sha256")
    parser.add_argument(
        "--expected-integrator", choices=("explicit", "implicit")
    )
    parser.add_argument("--expected-coarse-sha256")
    parser.add_argument("--expected-fine-sha256")
    parser.add_argument("--reference-log", type=Path)
    parser.add_argument("--candidate-log", type=Path, action="append")
    parser.add_argument("--final-time", type=float, default=1.0e-6)
    parser.add_argument(
        "--reconstruction",
        choices=("pcm", "characteristic_plm"),
        default="pcm",
    )
    parser.add_argument("--limiter", choices=("minmod", "mc"))
    args = parser.parse_args()

    if (args.reference_log is None) != (args.candidate_log is None):
        raise SystemExit("--reference-log and --candidate-log must be used together")

    if args.final_time <= 0.0:
        raise SystemExit("--final-time must be positive")
    if args.checkpoint is not None:
        validate_checkpoint(
            args.checkpoint,
            args.final_time,
            args.reconstruction,
            args.limiter,
            args.expected_schema,
            args.expected_bundle_sha256,
            args.expected_integrator,
        )
        checkpoint_digest = digest(args.checkpoint.read_bytes())
        if (
            args.expected_checkpoint_sha256 is not None
            and checkpoint_digest != args.expected_checkpoint_sha256.lower()
        ):
            raise SystemExit(
                "checkpoint SHA-256 mismatch: "
                f"expected={args.expected_checkpoint_sha256.lower()}, "
                f"actual={checkpoint_digest}"
            )
    elif any(
        value is not None
        for value in (
            args.expected_checkpoint_sha256,
            args.expected_schema,
            args.expected_bundle_sha256,
            args.expected_integrator,
        )
    ):
        raise SystemExit("checkpoint expectations require --checkpoint")

    coarse_reference = level_path(args.reference_prefix, "coarse")
    fine_reference = level_path(args.reference_prefix, "fine")
    coarse_digest = digest(coarse_reference.read_bytes())
    fine_digest = digest(fine_reference.read_bytes())
    if (
        args.expected_coarse_sha256 is not None
        and coarse_digest != args.expected_coarse_sha256.lower()
    ):
        raise SystemExit(
            "coarse SHA-256 mismatch: "
            f"expected={args.expected_coarse_sha256.lower()}, "
            f"actual={coarse_digest}"
        )
    if (
        args.expected_fine_sha256 is not None
        and fine_digest != args.expected_fine_sha256.lower()
    ):
        raise SystemExit(
            "fine SHA-256 mismatch: "
            f"expected={args.expected_fine_sha256.lower()}, "
            f"actual={fine_digest}"
        )
    for prefix in args.candidate_prefix:
        require_equal(
            coarse_reference, level_path(prefix, "coarse"), f"{prefix} coarse"
        )
        require_equal(fine_reference, level_path(prefix, "fine"), f"{prefix} fine")
        print(f"{prefix}: exact_coarse_and_fine=PASS")
    if args.reference_log is not None:
        reference_diagnostics = transport_log_diagnostics(args.reference_log)
        for candidate_log in args.candidate_log:
            candidate_diagnostics = transport_log_diagnostics(candidate_log)
            if candidate_diagnostics != reference_diagnostics:
                raise SystemExit(
                    f"{candidate_log}: cumulative transport diagnostics mismatch: "
                    f"expected={reference_diagnostics}, actual={candidate_diagnostics}"
                )
            print(f"{candidate_log}: exact_transport_diagnostics=PASS")
    print(f"coarse_sha256={coarse_digest}")
    print(f"fine_sha256={fine_digest}")
    print("mpi_amr_reactive_3d_exact_match=PASS")


if __name__ == "__main__":
    main()
