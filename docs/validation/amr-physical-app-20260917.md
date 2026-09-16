# Serial AMR physical-boundary application and restart — 2026-09-17

Source: parent `1c5af16` plus this change. Source hashes, exact artifact directory
and logs: `build/amr-physical-app-checks/`. Reused GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration. Fixed full-H2/O2 serial frontend exercised;
the selected frontend shares implementation but was not separately qualified.

Serial frontends opt into physical AMR configuration. MPI callers retain the
periodic-only default. The application resolves per-face reservoir states,
including optional named composition, and forwards face/wall settings through
the full split. Inlet acoustic/diffusion and isothermal-wall thermal limits
constrain the timestep. Composite conservation subtracts accumulated boundary
transfer; element diagnostics convert that transfer to H/O/N amounts. The
final log reports all conserved boundary components, while existing level CSVs
retain their format.

Physical static schemas 10–13 and dynamic schemas 14–17 append a canonical
boundary-settings signature and cumulative ledger. The signature includes faces,
wall temperatures/velocities, reservoir temperatures/pressures/velocities,
explicit inlet composition and fixed base composition used when inlet composition
is omitted. Existing selected context also checks its base composition.
Readers retain private state and budget until the whole checkpoint is accepted.
Dynamic readers recover geometry and counters before committing either level.
Periodic schemas 2–9 retain their existing layout. Fixed chemistry-enabled
checkpoint restrictions are unchanged.

The new application regression runs uninterrupted, stopped-after-step-one and
resumed versions of three cases: moving cold wall/outflow; the same boundaries
with a temperature-tagged patch enlarged on a 6-cubed coarse grid; and prescribed
reservoir inflow/outflow. Coarse and fine CSVs match byte-for-byte, as does the
printed cumulative boundary ledger. The stopped wall run reaches
7.7853964487436298e-9 s, before final time 1e-8 s, so restart executes further
steps. Changing wall temperature or inlet temperature rejects restart before
writing CSVs. The sample is `cases/amr_reactive_3d/physical_wall.nml`.

Completed commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^amr_physical_boundary_restart_3d$' --output-on-failure
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_reactive_3d_checkpoint --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_reactive_3d_checkpoint$' --output-on-failure
python3 -m py_compile tools/check_amr_physical_boundary_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Application test: 1/1 passed in 47.29 s. Existing checkpoint unit test: 1/1 passed
in 0.11 s. No selected tests skipped. Provenance, Python syntax, project contract
and diff checks passed. An overlong new Fortran call was wrapped during build.
The first dynamic fixture used a disallowed zero tagging threshold; changing it
to 1 K exercises the intended patch enlargement. The failed fixture log is retained
as `ctest-fixture-failure.log`. No solver tolerance was changed.

This is numerical application/restart evidence, not physical validation or a
full configuration matrix. MPI physical AMR boundaries, fine patches touching
the domain, and general multipatch/multilevel AMR remain open.
