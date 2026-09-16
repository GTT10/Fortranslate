# Serial boundary-touching AMR application/restart — 2026-09-17

Source: parent `a08ce5a` plus this change. Logs and source hashes:
`build/amr-boundary-app-checks/`. Reused GNU Fortran 13.3.0 Debug build.

Serial fixed and selected frontends now opt into boundary-touching static
patches. The configuration reader keeps this permission separate from physical
coarse boundaries, so MPI/default readers retain interior patches. Bounds must
remain inside the domain, touched faces must be nonperiodic, and boundary-touching
dynamic regridding is rejected before stepping. The serial application accepts
the validated geometry and includes fine-state wall thermal timestep limits,
scaled back by the ratio-squared transport subcycle count.

Static checkpoint shape/configuration validation now accepts the same physical
geometry. Existing schemas already store all required bounds and physical and
chemistry context, so no serialization layout changed. Periodic seam and
boundary-touching dynamic checkpoint contexts remain rejected. Reads retain
their private candidate state and final commit behavior.

The application regression checks a full-H2/O2 implicit reacting wall patch
with molecular transport and an elementary explicit reacting three-face corner
without transport. Both run uninterrupted and through stop/restart. Coarse and
fine CSVs and boundary-ledger output match exactly. Inactive-chemistry controls
differ from the reacting output. Changed wall data, chemistry tolerances,
integrator, reaction rates and model/source identity are rejected before output;
attempted dynamic regridding is rejected with its static-patch diagnostic.
Configuration tests retain default-reader rejection and exercise explicit opt-in,
periodic seam, invalid bounds and dynamic restrictions. Existing checkpoint
unit tests are also rerun.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_reactive_3d test_amr_reactive_3d_config test_amr_reactive_3d_checkpoint --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(amr_boundary_touching_restart_3d|unit_amr_reactive_3d_(config|checkpoint))$' --output-on-failure
python3 -m py_compile tools/check_amr_physical_boundary_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 3/3 tests passed, no skips (application 44.79 s, config 0.00 s,
checkpoint 0.10 s; total 44.92 s). Provenance, Python syntax, project contract
and diff checks passed. An initial overlong Fortran call was corrected after
the compiler rejected it; its log is `build-line-length-failure.log`.
Application artifacts:
`build/spray-mpi-migration-gnu-debug/tests/amr_boundary_touching_restart_3d/run-7jifii6m/`.
Runnable example: `cases/amr_reactive_3d/physical_touching_wall.nml`.

MPI boundary-touching stepping and dynamic regridding remain open. The selected
frontend shares the new opt-in/application path but was not independently
qualified in this run. This is numerical restart evidence, not physical validation.
