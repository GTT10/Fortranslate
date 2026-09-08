#!/usr/bin/env python3
"""Check generic selected-mechanism static two-level 3D AMR output."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


BASE_FIELDS = (
    "time",
    "x",
    "y",
    "z",
    "rho",
    "u",
    "v",
    "w",
    "pressure",
    "temperature",
    "rhoE",
    "rhou",
    "rhov",
    "rhow",
)
EULER_FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoE")


def close(left: float, right: float, tolerance: float = 2.0e-12) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=1.0e-15)


def relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(1.0, abs(left), abs(right))


def expected_fields(species: tuple[str, ...]) -> tuple[str, ...]:
    return (
        *BASE_FIELDS,
        *(f"Y_{name}" for name in species),
        *(f"rhoY_{name}" for name in species),
    )


def read_rows(
    path: Path, species: tuple[str, ...]
) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != expected_fields(species):
            raise AssertionError(f"unexpected selected AMR CSV schema: {path}")
        rows = [
            {name: float(value) for name, value in row.items()}
            for row in reader
        ]
    if not rows or any(
        not math.isfinite(value) for row in rows for value in row.values()
    ):
        raise AssertionError(f"nonfinite or empty selected AMR CSV: {path}")
    return rows


def validate_grid(
    rows: list[dict[str, float]],
    species: tuple[str, ...],
    shape: tuple[int, int, int],
    final_time: float,
    label: str,
) -> tuple[list[float], list[float], list[float], float, float]:
    nx, ny, nz = shape
    if len(rows) != nx * ny * nz:
        raise AssertionError(
            f"{label}: expected {nx * ny * nz} rows, found {len(rows)}"
        )
    axes = [sorted({row[name] for row in rows}) for name in ("x", "y", "z")]
    if tuple(len(axis) for axis in axes) != shape:
        raise AssertionError(f"{label}: unexpected coordinate topology")
    for axis_name, axis in zip(("x", "y", "z"), axes):
        if len(axis) > 1:
            spacing = axis[1] - axis[0]
            if spacing <= 0.0 or any(
                not close(axis[index] - axis[index - 1], spacing)
                for index in range(2, len(axis))
            ):
                raise AssertionError(f"{label}: nonuniform {axis_name} spacing")

    relation_error = 0.0
    closure_error = 0.0
    for index, row in enumerate(rows):
        i = index % nx
        j = (index // nx) % ny
        k = index // (nx * ny)
        if not all(
            close(actual, expected)
            for actual, expected in zip(
                (row["x"], row["y"], row["z"]),
                (axes[0][i], axes[1][j], axes[2][k]),
            )
        ):
            raise AssertionError(f"{label}: rows are not in k-j-i order")
        if not close(row["time"], final_time):
            raise AssertionError(f"{label}: unexpected output time")
        if min(row["rho"], row["pressure"], row["temperature"], row["rhoE"]) <= 0:
            raise AssertionError(f"{label}: nonphysical thermodynamic state")
        relation_error = max(
            relation_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
        )
        mass_fraction_sum = sum(row[f"Y_{name}"] for name in species)
        species_density_sum = sum(row[f"rhoY_{name}"] for name in species)
        if min(row[f"Y_{name}"] for name in species) < 0.0:
            raise AssertionError(f"{label}: negative species mass fraction")
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )
        relation_error = max(
            relation_error,
            *(
                relative_error(
                    row[f"rhoY_{name}"], row["rho"] * row[f"Y_{name}"]
                )
                for name in species
            ),
        )
    if relation_error > 2.0e-11 or closure_error > 2.0e-11:
        raise AssertionError(
            f"{label}: relation={relation_error}, closure={closure_error}"
        )
    return axes[0], axes[1], axes[2], relation_error, closure_error


def validate_fine_alignment(
    coarse_axes: tuple[list[float], list[float], list[float]],
    fine_axes: tuple[list[float], list[float], list[float]],
    lower: tuple[int, int, int],
    upper: tuple[int, int, int],
    ratio: int,
) -> None:
    for direction, (coarse, fine) in enumerate(zip(coarse_axes, fine_axes)):
        if len(coarse) < 2 or len(fine) < 2:
            raise AssertionError("selected AMR checker requires multi-cell axes")
        coarse_spacing = coarse[1] - coarse[0]
        fine_spacing = fine[1] - fine[0]
        if not close(fine_spacing, coarse_spacing / ratio):
            raise AssertionError("fine spacing does not match refinement ratio")
        first = coarse[lower[direction] - 1] - 0.5 * coarse_spacing
        first += 0.5 * fine_spacing
        last = coarse[upper[direction] - 1] + 0.5 * coarse_spacing
        last -= 0.5 * fine_spacing
        if not close(fine[0], first) or not close(fine[-1], last):
            raise AssertionError("fine coordinates do not align with patch faces")


def row_at(
    rows: list[dict[str, float]], nx: int, ny: int, i: int, j: int, k: int
) -> dict[str, float]:
    return rows[(k * ny + j) * nx + i]


def average_down_error(
    coarse: list[dict[str, float]],
    fine: list[dict[str, float]],
    species: tuple[str, ...],
    coarse_shape: tuple[int, int, int],
    lower: tuple[int, int, int],
    upper: tuple[int, int, int],
    ratio: int,
) -> float:
    nx, ny, _ = coarse_shape
    fine_nx = (upper[0] - lower[0] + 1) * ratio
    fine_ny = (upper[1] - lower[1] + 1) * ratio
    fields = (*EULER_FIELDS, *(f"rhoY_{name}" for name in species))
    error = 0.0
    for coarse_k in range(lower[2] - 1, upper[2]):
        local_k = coarse_k - (lower[2] - 1)
        for coarse_j in range(lower[1] - 1, upper[1]):
            local_j = coarse_j - (lower[1] - 1)
            for coarse_i in range(lower[0] - 1, upper[0]):
                local_i = coarse_i - (lower[0] - 1)
                parent = row_at(coarse, nx, ny, coarse_i, coarse_j, coarse_k)
                for field in fields:
                    child_sum = 0.0
                    for child_k in range(local_k * ratio, (local_k + 1) * ratio):
                        for child_j in range(
                            local_j * ratio, (local_j + 1) * ratio
                        ):
                            for child_i in range(
                                local_i * ratio, (local_i + 1) * ratio
                            ):
                                child_sum += row_at(
                                    fine,
                                    fine_nx,
                                    fine_ny,
                                    child_i,
                                    child_j,
                                    child_k,
                                )[field]
                    error = max(
                        error,
                        relative_error(parent[field], child_sum / ratio**3),
                    )
    return error


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--i-lower", type=int, required=True)
    parser.add_argument("--i-upper", type=int, required=True)
    parser.add_argument("--j-lower", type=int, required=True)
    parser.add_argument("--j-upper", type=int, required=True)
    parser.add_argument("--k-lower", type=int, required=True)
    parser.add_argument("--k-upper", type=int, required=True)
    parser.add_argument("--ratio", type=int, required=True)
    parser.add_argument("--final-time", type=float, required=True)
    parser.add_argument("--require-uniform", action="store_true")
    parser.add_argument("--activity-species")
    parser.add_argument("--initial-mass-fraction", type=float)
    parser.add_argument("--minimum-change", type=float, default=0.0)
    args = parser.parse_args()

    species = tuple(args.species)
    if len(species) != len(set(species)):
        raise AssertionError("selected species names must be unique")
    lower = (args.i_lower, args.j_lower, args.k_lower)
    upper = (args.i_upper, args.j_upper, args.k_upper)
    if args.ratio < 2 or any(lo < 2 or hi < lo for lo, hi in zip(lower, upper)):
        raise AssertionError("invalid selected AMR patch contract")
    coarse_shape = (args.nx, args.ny, args.nz)
    fine_shape = tuple((hi - lo + 1) * args.ratio for lo, hi in zip(lower, upper))
    coarse = read_rows(args.coarse, species)
    fine = read_rows(args.fine, species)
    coarse_result = validate_grid(
        coarse, species, coarse_shape, args.final_time, "coarse"
    )
    fine_result = validate_grid(fine, species, fine_shape, args.final_time, "fine")
    validate_fine_alignment(
        (coarse_result[0], coarse_result[1], coarse_result[2]),
        (fine_result[0], fine_result[1], fine_result[2]),
        lower,
        upper,
        args.ratio,
    )
    synchronization_error = average_down_error(
        coarse, fine, species, coarse_shape, lower, upper, args.ratio
    )
    if synchronization_error > 2.0e-11:
        raise AssertionError(f"average-down mismatch {synchronization_error}")

    uniform_error = 0.0
    if args.require_uniform:
        fields = (*EULER_FIELDS, *(f"rhoY_{name}" for name in species))
        for rows in (coarse, fine):
            reference = rows[0]
            uniform_error = max(
                uniform_error,
                *(
                    relative_error(row[field], reference[field])
                    for row in rows
                    for field in fields
                ),
            )
        if uniform_error > 2.0e-11:
            raise AssertionError(f"uniform selected AMR error {uniform_error}")

    activity = 0.0
    if args.activity_species is not None:
        if args.activity_species not in species:
            raise AssertionError("activity species is absent")
        if args.initial_mass_fraction is None:
            raise AssertionError("activity check requires initial mass fraction")
        field = f"Y_{args.activity_species}"
        activity = max(
            abs(row[field] - args.initial_mass_fraction)
            for rows in (coarse, fine)
            for row in rows
        )
        if activity <= args.minimum_change:
            raise AssertionError(f"selected chemistry activity is too small: {activity}")

    print(
        "Selected reactive 3D AMR: PASS "
        f"(coarse={len(coarse)}, fine={len(fine)}, "
        f"sync={synchronization_error:.3e}, "
        f"relation={max(coarse_result[3], fine_result[3]):.3e}, "
        f"closure={max(coarse_result[4], fine_result[4]):.3e}, "
        f"uniform={uniform_error:.3e}, activity={activity:.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
