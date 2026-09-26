# Coupled common-step MPI SDC — 2026-09-26

Implementation `26ee7928`, base `2452cad8`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Behavior

`advance_mpi_amr_tree_sdc_3d` now performs common-step SDC over owner-local
patches: initial restriction/recovery; old external/transport source evaluation;
source-aware hydro prediction; candidate-state new source evaluation; forced
reaction; and final hydro/transport flux synchronization and temperature recovery.
Each distributed stage finishes collectively before the next begins.

Optional reaction-source estimates persist between calls. All updates and
inventories remain private until the entire step succeeds. Iteration count,
chemistry/transport flags, optional estimate presence, species/reaction counts,
forcing, time and tolerance controls agree collectively. Boundary controls
retain the individual stages' collective validation.

External-source inventories use leaf ownership and compensated serial cell
order. Only the small accumulator is broadcast node by node; global field
arrays are not gathered. Boundary and external-source inventories remain
separate, and hydro plus trapezoidal transport fluxes feed the existing reflux.

## Completed focused checks

Artifacts: `build/sdc-mpi-loop-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`; configured Open MPI 4.1.6 through absolute
`/usr/bin/mpiexec`. No compiler/backend change.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Builds passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_split_3d
  --parallel 2`, logs `build.log`, `build-test.log`. Compiler warnings remain.
- New CTest `unit_mpi_amr_patch_tree_sdc_3d_4`: **1/1 passed**, 5.32 s
  (5.37 s total), `ctest.log`. Four ranks, nonuniform three-level tree, ownership
  `[0,1,0]` and two idle ranks. Two SDC iterations at 1e-10 s with chemistry,
  all transport flags and forcing `[20000,-30000,40000]`.
- Owned conserved state, temperature, reaction estimates, boundary/source
  inventories, maximum reflux correction and minimum transport theta match
  the serial coupled implementation exactly. Non-owned fields remain absent.
- An invalid hydro coefficient on an idle rank, encountered after initial
  transport evaluation, rejects collectively and preserves all step outputs.
- Command: `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_sdc_3d_4$'
  --output-on-failure`.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and `git diff --check` passed.

No application, checkpoint/restart, recursive MPI SDC, selected build,
chemistry/transport-disabled coupled case, nonperiodic coupled SDC, full matrix,
install audit, pinned full-solver comparison or physical validation was run.
Serial/MPI agreement is scoped to this periodic small-grid case.
