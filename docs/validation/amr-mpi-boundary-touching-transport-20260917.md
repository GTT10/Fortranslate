# Sparse MPI boundary-touching transport — 2026-09-17

Source: parent `8b3819d` plus this change. Logs and source hashes:
`build/amr-mpi-boundary-touching-transport-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

The MPI AMR molecular transport API now accepts patches touching nonperiodic
physical faces. Fine viscous ghosts reflect all momentum components at no-slip
walls, include prescribed tangential wall velocity, and adjust total energy
to preserve internal energy. The fine transport kernel uses the shared
physical wall heat/work treatment at fine spacing. Internal rank faces retain
neighbor state through the existing two-layer exchange; they are not treated
as physical walls. The physical kernel uses these two layers for its local
transport and positivity stencil.

The coarse/fine interface limiter skips physical sides without an exterior
coarse cell. Physical reflux replaces covered coarse domain-surface fluxes
with the subcycled fine fluxes in the global boundary ledger. The two Euler
candidates of SSPRK2 retain their private budgets until commit. The full
reaction/transport/hydro split inherits the new path. Periodic touching remains
unsupported; the previous interior transport path is retained.

The MPI hierarchy regression extends lower/upper corner and full-domain
fixtures to molecular transport and full reacting splitting, with isothermal
walls, nonzero tangential wall speeds and an inlet. It compares both levels,
temperatures and maximum reflux corrections exactly with serial stepping,
checks the existing boundary-ledger tolerance and composite conservation,
and verifies collective rollback for an invalid inlet on one rank. Four/eight
ranks include zero-fine ownership and one coarse column per rank. Existing
periodic and interior physical tests remain included.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_reactive_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed, no failures or skips (four ranks 255.48 s, eight
ranks 225.59 s; total 481.10 s). Provenance, project contract and diff checks
passed. The build succeeded with existing real-equality warnings. No Python
source changed. No external solver or full configuration matrix was run.

Boundary-touching MPI application inputs, fine-wall timestep integration and
regridding remain open. This is focused numerical evidence, not external
solver parity, production scaling or physical validation.
