# MPI SDC hydro stage — 2026-09-26

Implementation `a6ff68dd`, base `5bf7895f`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implemented contract

`compute_sdc_mpi_amr_tree_hydro_3d` communicates packed state/temperature and
signed conserved predictor sources through the existing owner-local tree ghost
exchange. PPM/PeleC CTU includes the source grid; the routine returns hydro
rates and dt-integrated face fluxes only on owners. Physical source halos use
outflow extrapolation. Gas halos retain their existing boundary handling.

The implementation reuses the ordinary MPI hydro validator and face loop via
an explicit stage mode. This mode skips initial restriction, temperature
recovery, state advancement and final reflux. A read-only wrapper uses private
local copies; it does not gather the global tree state. Every rank must agree
on the mode and controls. All rate/flux outputs publish after collective
success. Ordinary stepping keeps its existing behavior.

## Completed focused checks

Artifacts: `build/sdc-mpi-hydro-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`, configured Open MPI 4.1.6 through absolute
`/usr/bin/mpiexec`; no compiler/backend change.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d
  --parallel 2` (`build.log`, `build-test.log`). Compiler warnings remain.
- New CTest `unit_mpi_amr_patch_tree_sdc_hydro_3d_4`: **1/1 passed**, 1.66 s
  (1.70 s total), `ctest.log`. Four ranks, ownership `[0,1,0]`, two idle ranks.
  Nonuniform three-level state and signed momentum sources, periodic/outflow
  cases: rates and all three integrated flux arrays match serial exactly.
  Input state/temperature are unchanged; non-owned outputs remain unallocated.
  One owner's invalid energy rejects collectively and preserves prior outputs.
- Ordinary `unit_mpi_amr_patch_tree_ppm_3d_2`: **1/1 passed**, 4.27 s
  (4.31 s total), `ctest-ppm.log`. This checks the changed shared update path.
- Commands use `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug -R '^TEST_NAME$' --output-on-failure`.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and final `git diff --check` passed. A trailing space
  found by the initial diff check was removed after the build; no numerical
  source change followed the successful tests.

No MPI SDC transport/iteration integration, level-time interpolation, SDC
application, selected build, restart/rank change, full matrix, install audit,
full-solver pinned comparison or physical validation was performed. The new
source-aware MPI checks cover periodic and outflow boundaries specifically.
