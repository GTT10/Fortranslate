# Distributed periodic spray gas advection — 2026-09-17

Tested parent: `af2e58e` plus this change. Build/test/preflight logs and source
hashes: `build/spray-mpi-hydro-checks/`.

The x-slab gas update exchanges two layers with adjacent ranks, applies existing
characteristic PLM fluxes and SSPRK2, and commits only after all owner temperatures
recover. Optional hydro arguments in the particle/reaction wrapper enable
P/2-C/2-H-C/2-P/2 ordering with complete state rollback. Gas molecular transport,
SGS, distributed wall/breakup events and the MPI CLI remain open.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; unchanged configuration described in
`spray-mpi-migration-20260916.md`.

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) pass in 9.45 s, zero failures/skips. Added uneven
x-block tests compare nonuniform 3D periodic gas fields and recovered temperatures
against the serial characteristic-PLM SSPRK2 solver at absolute tolerance 1e-10,
using HLLC/MC. They verify actual evolution, collective rejection of a negative
density on one rank, and mass/momentum/energy conservation for the combined
particle/reaction/hydro split at relative 2e-12 with 1e-12 inventory scale floor.
Earlier MPI particle/reaction cases also pass. Provenance, contract and whitespace
checks pass; no Python source changed.

This is focused periodic x-slab evidence, not arbitrary Cartesian gas decomposition,
full configuration/solver-choice coverage, a complete MPI spray application, or
physical/scaling qualification. The combined example has cold nitrogen bulk gas;
active 1400 K reaction evidence remains the separate particle/reaction case.
