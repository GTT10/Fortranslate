#!/usr/bin/env python3
from pathlib import Path
import csv
import sys

if len(sys.argv) < 3:
    raise SystemExit('usage: compare_mpi_multispecies.py reference.csv candidate.csv [...]')

def read(path):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames
        rows = [[float(row[name]) for name in names] for row in reader]
    return names, rows

names, reference = read(sys.argv[1])
for filename in sys.argv[2:]:
    candidate_names, candidate = read(filename)
    if candidate_names != names or len(candidate) != len(reference):
        raise AssertionError('MPI output shape mismatch')
    maximum = 0.0
    for left, right in zip(reference, candidate):
        for a, b in zip(left, right):
            maximum = max(maximum, abs(a - b) / max(1.0, abs(a), abs(b)))
    print(f'{filename}: maximum_relative_difference={maximum:.16e}')
    if maximum > 5.0e-13:
        raise AssertionError('MPI multispecies rank-count parity failed')
