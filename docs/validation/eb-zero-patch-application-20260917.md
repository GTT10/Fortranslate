# Zero-patch EB application and restart — 2026-09-17

Implementation: `962d490cd3c03f285f5b347d5e8c6abc47cfde2a`, based on
`45ed5f6`. No update kernel or mechanism data changed.

## Behavior

Fixed and selected two-level EB applications accept `patch_count=0` and
scheduled `regrid_patch_count(event)=0`. Coarse-only initialization, complete
fine-grid removal, continued reaction/transport/hydro, later refinement,
automatic timestep selection, cumulative inventories and output all handle
empty topology. First-patch pointers are explicitly null when no fine grid
exists; empty topology never dereferences them.

Static coarse-only checkpoints use a separate `PELEF_EB_COARSE_CHECKPOINT_V1`
through V4 envelope for the existing fixed/automatic and common/subcycling
mode combinations. Scheduled checkpoints retain their existing envelope and
policy but omit fine records at zero patch count. Coarse-only CSV output
contains active level-0 cells without a fine-grid row. Existing one/multiple
patch output layouts and static checkpoints are retained.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`. Fixed elementary H2/O2 is explicit; selected
full H2/O2 is implicit. Cases reuse the scheduled planar geometry, reaction,
transport, timestep and composition controls, replacing the first event count
with zero. All relative tolerances are unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_two_level_3d test_selected_full_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(regression_(selected_)?coarse_gap_eb_two_level_3d|regression_reactive_eb_two_level_3d)$' --output-on-failure --output-log build/eb-zero-app-checks/ctest.log
python3 -m py_compile tools/check_scheduled_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, total 33.58 s, zero failures/skips:

- Existing static two-level application: 7.66 s.
- Fixed coarse-gap application: 8.15 s.
- Selected coarse-gap application: 17.70 s.

Both runtime tests cover fixed common stepping and automatic subcycling for
2→0→2 patch transitions. Restart before removal, while coarse-only, and after
refinement reproduces final state CSV and budget CSV exactly. The empty-stage
output has exactly 48 fluid coarse cells; final refined output has 160 cells.
Changed schedules and event indices reject without overwriting prior outputs.

Each runtime also exercises static coarse-only initialization, V4 checkpoint
continuation and scheduled refinement starting from zero fine patches. Maximum
reported composite residual among the four 2→0→2 executions is
6.6909451336392976e-16; static coarse-only fixed/selected residuals are
7.0662624099902823e-16 and 2.1285011971799971e-16 respectively.

Provenance, Python syntax, contract and whitespace checks passed. Logs:
`build/eb-zero-app-checks/`; per-invocation inputs/checkpoints/CSVs:
`build/selected-tree-mpi-debug/tests/coarse_gap_eb_two_level_3d/` and
`selected_coarse_gap_eb_two_level_3d/`.

Temperature-policy dispatch and policy-bound adaptive restart remain open;
this change supplies their empty-topology endpoint. General multilevel/MPI EB,
external PeleC parity, convergence and physical validation are not established
by these focused application regressions. No installed audit or full matrix
was rerun.
