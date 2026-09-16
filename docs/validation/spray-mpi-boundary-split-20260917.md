# Nonperiodic MPI coupled splitting — 2026-09-17

Source: parent `577523b` plus this change. Source hashes and logs:
`build/spray-mpi-boundary-split-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`advance_mpi_mixture_reaction` accepts physical boundary labels, inlet state and
temperature, wall temperature and velocity, and a boundary-change output.
The same boundary data now reach particle half-steps, hydro, molecular transport
and SGS in the existing P/2-C/2-T/2-L/2-H-L/2-T/2-C/2-P/2 order. The boundary
ledger sums already-global stage increments without reducing them a second time.
It is positive into gas plus particles and replicated across ranks. Chemical
species conversion is not a boundary flux. The default periodic branch keeps its
previous stage paths; physical-boundary mode is agreed collectively.

The outer candidate arrays and ID watermark commit only after every stage succeeds.
Late failure leaves caller gas, temperature, particles, component masses and IDs
unchanged and returns zero substeps and boundary change. Mixed deposition/splash
and film evolution are not yet included; stationary reflection/open-face removal
remain the particle boundary behavior. Application inputs/checkpoints still need
the nonperiodic integration.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

2/2 focused entries passed in 16.98 s with no failures or skips. The new case uses
uneven x slabs, pure N2 gas, open x faces, no-slip/isothermal y walls, active
hydro/transport/SGS and an escaping two-component parcel. Pure N2 isolates the
boundary inventory test while still invoking chemistry. All gas-plus-liquid
inventories follow the summed ledger within 2e-12 of initial absolute inventories.
Gas changes measurably and the escaped parcel is removed. An invalid wall
temperature on the final rank triggers failure after the first particle half-step;
all outer state, including the already-escaped candidate parcel, is restored exactly.
Existing reacting periodic and individual nonperiodic API cases also pass.
Provenance, contract and whitespace checks passed; no Python source changed.

No full configuration matrix or physical LES/spray validation is claimed.
