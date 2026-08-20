#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--nx", type=int, default=24)
    parser.add_argument("--ny", type=int, default=24)
    args = parser.parse_args()
    rows = list(csv.DictReader(Path(args.input).open(newline="")))
    if len(rows) != args.nx * args.ny:
        raise SystemExit(f"unexpected row count: {len(rows)}")
    required = ["rho", "pressure", "temperature", "u", "v", "Y_H2", "Y_O2", "Y_OH", "Y_H2O", "Y_N2"]
    for name in required:
        if name not in rows[0]:
            raise SystemExit(f"missing column {name}")
    values = {name: [float(row[name]) for row in rows] for name in required}
    if any(not math.isfinite(value) for column in values.values() for value in column):
        raise SystemExit("non-finite reactive 2D output")
    if min(values["rho"]) <= 0.0 or min(values["pressure"]) <= 0.0 or min(values["temperature"]) <= 0.0:
        raise SystemExit("reactive 2D positivity failure")
    species_names = ["Y_H2", "Y_H", "Y_O", "Y_O2", "Y_OH", "Y_H2O", "Y_N2"]
    maximum_closure = 0.0
    for row in rows:
        total = sum(float(row[name]) for name in species_names)
        maximum_closure = max(maximum_closure, abs(total - 1.0))
    if maximum_closure > 2.0e-11:
        raise SystemExit(f"composition closure failure: {maximum_closure}")
    pressure_span = max(values["pressure"]) - min(values["pressure"])
    max_speed = max(math.hypot(u, v) for u, v in zip(values["u"], values["v"]))
    if pressure_span < 1.0:
        raise SystemExit(f"pressure response too small: {pressure_span}")
    if max_speed < 1.0e-3:
        raise SystemExit(f"velocity response too small: {max_speed}")
    if max(values["Y_H2O"]) < 1.0e-4:
        raise SystemExit("water-production signature missing")
    if max(values["Y_OH"]) < 1.0e-5:
        raise SystemExit("OH-production signature missing")
    print(f"rows={len(rows)} pressure_span={pressure_span:.12e} max_speed={max_speed:.12e} max_closure={maximum_closure:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
