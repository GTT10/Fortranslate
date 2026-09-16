# Boundary-touching AMR reaction/hydro split — 2026-09-17

Source: parent `e623edf` plus this change. Logs and source hashes:
`build/amr-boundary-reaction-checks/`. Reused GNU Fortran 13.3.0 Debug build
with bounds and floating-point checks.

The serial AMR chemistry source and reaction/hydro/reaction split now accept
valid boundary-touching patches. The source is cell-local and retains its
existing private two-level candidates, average-down and temperature recovery.
The split forwards physical data to the boundary-touching hydro implementation
and commits its boundary ledger only after the second reaction half-step succeeds.
The state-based CFL helper also accepts these patches. Shape, thermodynamic,
reaction and integrator checks remain intact; the hydro phase still rejects
unsupported periodic seams.

The existing physical split fixture now runs on the interior patch, an inflow
face patch and a three-face corner with inflow and reflecting walls. PCM and
characteristic PLM cases verify nonzero chemical conversion, Euler and H/O/N
balance after subtracting boundary transfers, exact equality with explicit
reaction/hydro/reaction phases, and direct hydro equality with chemistry disabled.
An invalid inlet triggers hydro failure after the first reaction half-step;
both levels and all output diagnostics must roll back. The CFL helper must
return a positive bound for every tested layout. Existing source and hydro
fixtures remain in the same executable.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_reactive_chemistry_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_reactive_chemistry_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 1/1 test passed, no skips (56.66 s; CTest total 56.69 s). Provenance,
project contract and diff checks passed. No Python source changed.

Boundary-touching molecular transport, MPI, regridding, checkpoints and
application inputs remain open. The new reacting fixtures use the fixed
elementary mechanism; they do not establish detailed-fuel ignition or physical
validation. No integrator or solver tolerance changed.
