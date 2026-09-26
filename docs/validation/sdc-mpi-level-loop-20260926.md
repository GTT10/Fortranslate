# Coupled MPI level SDC — 2026-09-26

Implementation `3dd924b9`, base `ea08e6b5`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`advance_mpi_amr_tree_sdc_3d` accepts level selection, parent state/source
endpoints, start/end fractions, integrated-flux output and old/new predictor
outputs. All ranks agree on selection and optional argument groups. Fine levels
require parent endpoints; level mode requires the integrated-flux output.

Hydro samples parent state at interval start and predictor sources at midpoint.
Transport samples state at start/end. Reaction and temperature recovery select
only the active level. Inactive state, temperature and persistent reaction
estimates are retained. Level mode omits initial restriction and final reflux;
an outer recursive caller owns synchronization. Root-level boundary inventory
and active-leaf external inventory are separate. Exported fluxes include hydro
and trapezoidal transport; predictors include external, transport and reaction
estimates. Outputs publish after collective success. Common-step remains default.

## Focused evidence

Artifacts: `build/sdc-mpi-level-loop-checks/` (local). GNU Fortran 13.3.0,
Open MPI 4.1.6, existing Ninja Debug tests-enabled build
`build/spray-amr-reactive-debug`; compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_split_3d
  --parallel 2` (`build.log`). Initial build failed on a newly extended line
  exceeding 132 columns; wrapped and rebuilt. Retained in
  `initial-build-failure.log`. Compiler warnings remain.
- CTest **1/1 passed**, 8.57 s (8.62 s total):
  `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir
  build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_sdc_3d_4$'
  --output-on-failure` (`ctest.log`).
- Nonuniform three-level tree, owners `[0,1,0]`, four ranks including two idle.
  Two SDC iterations, reactions, all transport flags, external forcing,
  dt=1e-10 s. Existing common-step exact serial/MPI parity and collective
  invalid-hydro rollback passed. Added deepest-level interval uses parent
  fractions `[.25,.5]` to `[.5,1]`, changed gas and signed source endpoints.
- Level state, temperature, reaction estimates, all three integrated fluxes,
  old/new predictors, boundary/external inventories and diagnostics match
  serial exactly. Inactive state/temperature remain exact and fluxes are
  allocated only for the selected owner.
- Python syntax for `tools/check_project_contract.py`, project contract and
  `git diff --check` passed.

No application trajectory, selected build, recursive MPI SDC, full matrix,
physical validation or installed audit was run. Recursive child timesteps and
synchronization still need to be connected to this interval API.
