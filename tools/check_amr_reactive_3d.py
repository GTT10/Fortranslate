#!/usr/bin/env python3
"""Check the public static two-level reactive 3D AMR output contract."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise AssertionError(f"missing CSV header: {path}")
        return [{name: float(value) for name, value in row.items()} for row in reader]


def relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(1.0, abs(left), abs(right))


def validate_level(
    rows: list[dict[str, float]], nx: int, ny: int, nz: int, time: float
) -> tuple[list[str], float, float]:
    if len(rows) != nx * ny * nz:
        raise AssertionError(f"expected {nx * ny * nz} rows, found {len(rows)}")
    species = sorted(name[2:] for name in rows[0] if name.startswith("Y_"))
    if not species:
        raise AssertionError("no species columns")
    relationship_error = 0.0
    closure_error = 0.0
    for row in rows:
        if not all(math.isfinite(value) for value in row.values()):
            raise AssertionError("nonfinite AMR CSV value")
        if relative_error(row["time"], time) > 2.0e-13:
            raise AssertionError("unexpected AMR output time")
        if min(row["rho"], row["pressure"], row["temperature"], row["rhoE"]) <= 0:
            raise AssertionError("nonphysical AMR output state")
        relationship_error = max(
            relationship_error,
            relative_error(row["rhou"], row["rho"] * row["u"]),
            relative_error(row["rhov"], row["rho"] * row["v"]),
            relative_error(row["rhow"], row["rho"] * row["w"]),
        )
        mass_fraction_sum = sum(row[f"Y_{name}"] for name in species)
        species_density_sum = sum(row[f"rhoY_{name}"] for name in species)
        if min(row[f"Y_{name}"] for name in species) < 0.0:
            raise AssertionError("negative AMR species mass fraction")
        closure_error = max(
            closure_error,
            abs(mass_fraction_sum - 1.0),
            relative_error(species_density_sum, row["rho"]),
        )
        for name in species:
            relationship_error = max(
                relationship_error,
                relative_error(row[f"rhoY_{name}"], row["rho"] * row[f"Y_{name}"]),
            )
    if relationship_error > 2.0e-11 or closure_error > 2.0e-11:
        raise AssertionError(
            f"AMR state contract error: relation={relationship_error}, "
            f"closure={closure_error}"
        )
    return species, relationship_error, closure_error


def row_at(rows: list[dict[str, float]], nx: int, ny: int, i: int, j: int, k: int):
    return rows[(k * ny + j) * nx + i]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
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
    parser.add_argument("--time", type=float, required=True)
    args = parser.parse_args()

    coarse = read_rows(args.coarse)
    fine = read_rows(args.fine)
    fine_nx = (args.i_upper - args.i_lower + 1) * args.ratio
    fine_ny = (args.j_upper - args.j_lower + 1) * args.ratio
    fine_nz = (args.k_upper - args.k_lower + 1) * args.ratio
    coarse_species, coarse_relation, coarse_closure = validate_level(
        coarse, args.nx, args.ny, args.nz, args.time
    )
    fine_species, fine_relation, fine_closure = validate_level(
        fine, fine_nx, fine_ny, fine_nz, args.time
    )
    if coarse_species != fine_species:
        raise AssertionError("coarse/fine species layouts differ")

    conserved = [
        "rho",
        "rhou",
        "rhov",
        "rhow",
        "rhoE",
        *(f"rhoY_{name}" for name in coarse_species),
    ]
    synchronization_error = 0.0
    ratio_children = args.ratio**3
    for coarse_k in range(args.k_lower - 1, args.k_upper):
        local_k = coarse_k - (args.k_lower - 1)
        for coarse_j in range(args.j_lower - 1, args.j_upper):
            local_j = coarse_j - (args.j_lower - 1)
            for coarse_i in range(args.i_lower - 1, args.i_upper):
                local_i = coarse_i - (args.i_lower - 1)
                parent = row_at(
                    coarse, args.nx, args.ny, coarse_i, coarse_j, coarse_k
                )
                for name in conserved:
                    child_sum = 0.0
                    for child_k in range(
                        local_k * args.ratio, (local_k + 1) * args.ratio
                    ):
                        for child_j in range(
                            local_j * args.ratio, (local_j + 1) * args.ratio
                        ):
                            for child_i in range(
                                local_i * args.ratio, (local_i + 1) * args.ratio
                            ):
                                child_sum += row_at(
                                    fine,
                                    fine_nx,
                                    fine_ny,
                                    child_i,
                                    child_j,
                                    child_k,
                                )[name]
                    child_average = child_sum / ratio_children
                    synchronization_error = max(
                        synchronization_error,
                        relative_error(parent[name], child_average),
                    )
    if synchronization_error > 2.0e-11:
        raise AssertionError(
            f"coarse/fine average-down mismatch {synchronization_error}"
        )
    density_span = max(row["rho"] for row in fine) - min(
        row["rho"] for row in fine
    )
    if density_span <= 1.0e-3:
        raise AssertionError("public AMR case lost its nonuniform flow signal")
    print(
        "Reactive 3D AMR: PASS "
        f"(coarse={len(coarse)}, fine={len(fine)}, "
        f"sync={synchronization_error:.3e}, "
        f"relation={max(coarse_relation, fine_relation):.3e}, "
        f"closure={max(coarse_closure, fine_closure):.3e}, "
        f"density_span={density_span:.3e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
