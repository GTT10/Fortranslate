# MPI gas MOL application — 2026-09-27

Implementation `b8b22fba` adds optional `mol_mode=.true.` to `eb_tree` for the
MPI gas tree applications. Use `time_integrator='sdc'`,
`hydro_method='plm_ctu'`, `riemann_solver='pelec'` and explicit boundary faces.
The existing coupled integrator controls now select the MOL spatial/time
update in both common and recursive branches. The hydro-method setting is
retained for shared configuration validation; MOL does not execute CTU tracing.
Default mode remains false. Serial and particle applications reject enabled
MOL until those paths are connected.

Enabled checkpoints prepend `PELEF_EB_TREE_MOL_V1`; changing mode on restart
is rejected in both directions. Existing disabled-mode envelopes are retained.
The existing coupled hydro/transport timestep selection is used unchanged.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug selected heptane Simple
transport build, MPI/spray enabled, tests/SUNDIALS disabled:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected pelef_reactive_eb_tree_3d_selected --parallel 2 > build/mpi-mol-app-build.log 2>&1
python3 build/mpi-mol-app/run_checks_final.py > build/mpi-mol-app/driver-final.log 2>&1
python3 -m py_compile build/mpi-mol-app/run_checks_final.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

These checks passed. The driver records seven MPI launches in
`build/mpi-mol-app-final/result.json`: continuous MOL, first-step checkpoint,
one-to-two-rank restart, MOL-to-CTU rejection, recursive MOL, CTU checkpoint,
and CTU-to-MOL rejection. Its two negative launches must fail with the specific
mode-mismatch message and no final CSV. The 4³ heated-wall gas fixture uses
`cases/spray_eb_wall_heat_3d/tree.nml` and its level set without particles,
two `1e-7 s` steps, conduction and no chemistry. Final gas CSV bytes match
between continuous/restarted and common/recursive runs. This single-level
application fixture complements the separate nested recursive library check.
MPI binary SHA256:
`5a0835d329dfaf9c7d066f41386f4214cd0172271d5b94127eefd2c4f878c579`.

The first driver in `build/mpi-mol-app/run_checks.py` incorrectly assumed the
particle checkpoint suffix `.step1`; the gas application writes `state.chk`.
That run failed opening the restart and is not counted as passed. Its artifacts
remain; the corrected driver uses a separate output directory. X11 authorization
and executable-stack warnings were nonfatal. No full matrix, particle MOL,
serial MOL, full-duration Spray-EB or upstream temporal comparison was run.
