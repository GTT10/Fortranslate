# Periodic MPI spray molecular transport — 2026-09-17

Parent `0bac6eb` plus this change. Source hashes and build/test/preflight logs:
`build/spray-mpi-transport-checks/`.

The x-slab transport operator exchanges one gas layer and species-outflow limiter
factors, computes common limited interface fluxes, advances SSPRK2 and recovers
temperatures before collective commit. Four flags select viscosity, conduction,
species diffusion and barodiffusion. Optional transport flags in the coupled
wrapper enable P/2-C/2-T/2-H-T/2-C/2-P/2 ordering with complete interval rollback.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; existing MPI configuration retained
from `spray-mpi-migration-20260916.md`.

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) pass in 9.82 s, no failures/skips. Added nonuniform
water-vapor/nitrogen tests enable all four transport processes, compare conserved
fields against serial transport at `2e-12 * max(1, abs(reference))`, temperatures
within 1e-8 K and minimum limiter factor within 1e-14, and require actual evolution.
The combined particle/reaction/transport/hydro split preserves global mass,
momentum and energy at `2e-12 * max(1e-12, abs(initial inventory))`. One rank with
negative density causes collective rejection without changing fields/temperatures.
Earlier MPI cases also pass. Provenance, contract and whitespace checks pass;
no Python source changed.

The focused fixture does not force an active species limiter, qualify independent
flag combinations, prove MPI CVODE/scaling, or complete the MPI spray application.
No full matrix or physical qualification is claimed.
