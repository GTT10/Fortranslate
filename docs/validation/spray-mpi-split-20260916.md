# Collective computational particle splitting — 2026-09-16

Tested parent: `9ad99f3` plus this change. Source/test hashes and build/test/preflight
logs are in `build/spray-mpi-split-checks/`.

The MPI split API constructs local multiplicity-capped candidates, collectively
validates them, reserves globally disjoint child ID ranges and commits all local
particle/component arrays with the shared high-water mark. Per-drop component
columns inherit from parents; the existing splitter partitions extensive KH state.
Any rank's input/capacity rejection leaves every caller state and ID mark unchanged.
Input IDs must be locally increasing, globally unique and bounded by the mark.
This is computational splitting, not aerodynamic breakup or a distributed CFD driver.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; configuration/library boundary is retained
from spray-mpi-migration-20260916.md:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) passed in 0.81 s total, no failures/skips. New checks
cover one-rank capacity rollback, disjoint child IDs, preserved parent IDs,
multiplicity/KH count/KH shed-volume conservation and component inheritance.
Existing migration and ID reservation cases also pass. Provenance, contract and
whitespace checks passed; no Python source changed. No full matrix, integrated MPI
spray application or physical/scaling qualification was run.
