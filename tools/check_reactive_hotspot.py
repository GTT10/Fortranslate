#!/usr/bin/env python3
"""Structural regression checks for the reactive one-dimensional hotspot."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

SPECIES = ("H2", "H", "O", "O2", "OH", "H2O", "N2")

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    args = parser.parse_args()
    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(handle)]
    if len(rows) < 16:
        raise AssertionError("reactive hotspot output is too short")
    values = [value for row in rows for value in row.values()]
    if not all(math.isfinite(value) for value in values):
        raise AssertionError("reactive hotspot contains non-finite data")
    min_rho = min(row["rho"] for row in rows)
    min_p = min(row["pressure"] for row in rows)
    min_t = min(row["temperature"] for row in rows)
    max_t = max(row["temperature"] for row in rows)
    max_u = max(abs(row["u"]) for row in rows)
    pressure_span = max(row["pressure"] for row in rows) - min_p
    closure = max(abs(sum(row[f"Y_{name}"] for name in SPECIES) - 1.0) for row in rows)
    min_y = min(row[f"Y_{name}"] for row in rows for name in SPECIES)
    print(f"minimum_density={min_rho:.16e}")
    print(f"minimum_pressure={min_p:.16e}")
    print(f"minimum_temperature={min_t:.16e}")
    print(f"maximum_temperature={max_t:.16e}")
    print(f"maximum_velocity={max_u:.16e}")
    print(f"pressure_span={pressure_span:.16e}")
    print(f"maximum_closure_error={closure:.16e}")
    print(f"minimum_mass_fraction={min_y:.16e}")
    if min_rho <= 0.0 or min_p <= 0.0 or min_t <= 0.0:
        raise AssertionError("reactive hotspot lost thermodynamic positivity")
    if min_y < -2.0e-12 or closure > 2.0e-10:
        raise AssertionError("reactive hotspot composition is invalid")
    if max_u < 1.0e-3 or pressure_span < 1.0:
        raise AssertionError("reactive hotspot did not generate a pressure wave")
    if max_t - min_t < 10.0:
        raise AssertionError("reactive hotspot lost its thermal structure")
    print("reactive hotspot regression: PASS")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
