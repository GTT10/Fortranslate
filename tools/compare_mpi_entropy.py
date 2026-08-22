#!/usr/bin/env python3
from pathlib import Path
import csv
import math
import sys

if len(sys.argv) < 3:
    raise SystemExit('usage: compare_mpi_entropy.py reference.csv candidate.csv [...]')

def read(path):
    with Path(path).open(newline='') as handle:
        return [[float(value) for value in row.values()] for row in csv.DictReader(handle)]

reference = read(sys.argv[1])
for name in sys.argv[2:]:
    candidate = read(name)
    if len(candidate) != len(reference):
        raise AssertionError('row count differs')
    maximum = 0.0
    for left, right in zip(reference, candidate):
        for a, b in zip(left, right):
            scale = max(1.0, abs(a), abs(b))
            maximum = max(maximum, abs(a - b) / scale)
    print(f'{name}: maximum_relative_difference={maximum:.16e}')
    if maximum > 5.0e-13:
        raise AssertionError('MPI rank-count parity failed')
