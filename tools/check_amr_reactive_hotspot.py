#!/usr/bin/env python3
"""Structural and coverage checks for the 1D AMR reactive hotspot."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--domain-length", type=float, default=0.012)
    parser.add_argument("--refinement-ratio", type=int, default=2)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if len(rows) <= 32:
        raise AssertionError("AMR output does not contain refined coverage")
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise AssertionError("AMR hotspot contains non-finite data")
    levels = {int(row["level"]) for row in rows}
    if levels != {0, 1}:
        raise AssertionError(f"expected coarse and fine rows, found {levels}")
    coarse_dx = max(row["cell_dx"] for row in rows)
    fine_dx = min(row["cell_dx"] for row in rows)
    ratio = coarse_dx / fine_dx
    if abs(ratio - args.refinement_ratio) > 1.0e-12:
        raise AssertionError(f"unexpected refinement ratio {ratio}")
    coverage = math.fsum(row["cell_dx"] for row in rows)
    if abs(coverage - args.domain_length) > 2.0e-13:
        raise AssertionError(f"composite coverage is {coverage}")
    positions = [row["x"] for row in rows]
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise AssertionError("composite output positions are not ordered")

    min_rho = min(row["rho"] for row in rows)
    min_pressure = min(row["pressure"] for row in rows)
    min_temperature = min(row["temperature"] for row in rows)
    temperature_span = max(row["temperature"] for row in rows) - min_temperature
    pressure_span = max(row["pressure"] for row in rows) - min_pressure
    max_velocity = max(abs(row["u"]) for row in rows)
    closure = max(
        abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0)
        for row in rows
    )
    min_fraction = min(row[f"Y_{name}"] for row in rows for name in SPECIES)
    if min_rho <= 0.0 or min_pressure <= 0.0 or min_temperature <= 0.0:
        raise AssertionError("AMR hotspot lost thermodynamic positivity")
    if closure > 2.0e-10 or min_fraction < -2.0e-12:
        raise AssertionError("AMR hotspot composition is invalid")
    if temperature_span < 10.0:
        raise AssertionError("AMR hotspot lost its temperature structure")
    if pressure_span < 1.0 or max_velocity < 1.0e-3:
        raise AssertionError("AMR hotspot did not develop a pressure wave")

    print(f"rows={len(rows)}")
    print(f"coarse_dx={coarse_dx:.16e}")
    print(f"fine_dx={fine_dx:.16e}")
    print(f"coverage={coverage:.16e}")
    print(f"minimum_density={min_rho:.16e}")
    print(f"minimum_pressure={min_pressure:.16e}")
    print(f"minimum_temperature={min_temperature:.16e}")
    print(f"temperature_span={temperature_span:.16e}")
    print(f"pressure_span={pressure_span:.16e}")
    print(f"maximum_velocity={max_velocity:.16e}")
    print(f"maximum_closure_error={closure:.16e}")
    print("AMR reactive hotspot regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
