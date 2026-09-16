# MPI SGS with physical boundaries — 2026-09-17

Source: parent `7079423` plus this change. Logs and source hashes:
`build/spray-mpi-boundary-les-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

The supplied-ghost SGS kernel now accepts boundary labels, wall temperatures,
a boundary-rate output and a face-inclusion mask. The mask excludes internal
MPI x interfaces from the surface budget without changing their fluxes or
cell updates. Existing serial callers retain all six faces.

`advance_mpi_spray_les` accepts optional physical faces, inlet data, wall
velocity/temperature and a boundary-change output. It exchanges two ghost layers
at each SSPRK2 stage, reuses serial isothermal-wall heat flux and viscous ghost
velocities, retains the diffusion-step bound, and commits gas, temperature and
global ledger only after collective success. Ledger values are replicated and
positive into gas. Omitting faces retains the existing periodic slab path.
The caller provides consistent physical model/boundary data across ranks.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests passed 2/2 in 15.88 s with no failures or skips. Provenance,
contract and whitespace checks passed; no Python source changed.

New comparisons use uneven x partitions, nonuniform H2/O2/N2 composition,
periodic, wall/outflow and inlet/wall boundaries, 600 K wall temperature and
tangential wall motion. Conserved fields match serial within 2e-12 scaled
error, temperature within 1e-8 K and boundary ledgers within 1e-10 relative error
with a 1e-12 floor. Independent volume changes match surface budgets within
2e-12 of initial absolute inventories. Invalid wall temperature on one rank rejects
collectively, preserves gas/temperature exactly and returns zero ledger.
The first comparison failed because its serial API call swapped time-step and
SGS-coefficient positions; correcting the test call resolved that discrepancy.
The diagnostic failure is retained in `diagnostic.log`; tolerances were unchanged.

Nonperiodic particle stepping and application integration remain open. This is
numerical implementation evidence, not physical LES validation or a full matrix.
