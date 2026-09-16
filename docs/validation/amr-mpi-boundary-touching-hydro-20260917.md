# Sparse MPI boundary-touching hydro — 2026-09-17

Source: parent `cf27a5d` plus this change. Logs and source hashes:
`build/amr-mpi-boundary-touching-hydro-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

The sparse two-level hydro API now accepts patches touching nonperiodic domain
faces. Fine physical ghost coordinates are reflected or clamped into fine
state, or filled from the prescribed inlet. Remaining coarse/fine ghost
coordinates use the existing time-interpolated PCM/PLM coarse reconstruction.
The two-layer neighboring-rank fine exchange remains local; there is no
full-level state gather in stepping. Face and inlet metadata are collectively
validated through the coarse physical solver before fine communication.

The physical reflux path corrects uncovered interior coarse cells and replaces
covered physical coarse surface fluxes in the boundary ledger. Each surface
piece belongs to one rank, and its correction is summed once globally. Fine
subcycling, average-down, temperature recovery and transactional hierarchy
commit remain connected. The reaction/hydro split inherits this support and
keeps both levels unchanged if the hydro boundary input fails after reaction.
Touching periodic seams and touching patches without physical face data remain
unsupported. Molecular transport and MPI application inputs retain their
interior-patch restriction.

The existing physical hierarchy regression now also checks lower and upper
three-face corners using PCM/PLM hydro and active reaction/hydro splitting,
and a full-domain patch using PCM/PLM hydro. It compares all coarse/fine states,
temperatures and maximum reflux corrections exactly with serial stepping,
compares global boundary ledgers at the existing tolerance, checks composite
Euler/element budgets and invalid-one-rank-inlet rollback. Existing periodic,
interior physical transport/full splitting and distribution tests remain.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_reactive_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed, no failures or skips (four ranks 104.09 s, eight
ranks 91.16 s; total 195.28 s). Provenance, project contract and diff checks
passed. No Python source changed. No external solver or full configuration
matrix was run.

This is focused MPI/serial numerical evidence, not external solver parity,
production scaling or physical validation. General multipatch/multilevel AMR,
boundary-touching MPI transport, application input/restart and regridding
remain open.
