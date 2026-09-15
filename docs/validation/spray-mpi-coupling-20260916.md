# Periodic MPI mixed-particle coupling — 2026-09-16

Tested parent: `3da2f64` plus this change. Source hashes and build/test/preflight
logs: `build/spray-mpi-coupling-checks/`.

The distributed coupler now advances mixed particles through owner sampling,
common rate/CFL limits, frozen-gas drag/heating/evaporation and collective gas
feedback. Failed trials halve the common step. All interval state stays private
until surviving particles migrate successfully. Synchronous feedback differs from
serial particle-by-particle feedback; trajectory parity is not claimed.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; existing compiler/library configuration
is recorded in `spray-mpi-migration-20260916.md`.

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) passed in 0.84 s, no failures/skips. Added cases
use an initially empty particle rank, check drag-induced velocity reduction,
heating/evaporation with gas vapor production, and global gas-plus-particle mass,
momentum, total energy and species inventories. Conservation tolerance is
`2e-12 * max(1e-12, abs(initial inventory))` per component. A one-substep limit
rejects a longer interval and preserves initial gas, temperature, positions,
velocities and component masses. Existing primitive checks pass. Provenance,
project contract and whitespace checks passed; no Python source changed.

This does not independently force the feedback-halving branch, qualify changed-rank
trajectories, or provide a complete MPI reactive spray CLI. Distributed wall and
breakup events, physical/scaling validation and full configuration tests remain open.
