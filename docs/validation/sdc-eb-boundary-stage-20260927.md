# Physical boundaries for the EB SDC hydro stage — 2026-09-27

Implementation `8707dd86`, base `f5a788d0`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`compute_sdc_eb_boundary_hydro_3d` takes physical state, temperature and conserved
predictor rates with geometry extended by six cells per side. It shares the
existing physical boundary donor mapping and flattening, then forwards working
state and source halos to the SDC hydro stage. It returns physical-grid rates,
time-integrated flux densities and boundary/wall inventories without updating
state. Signed conserved rates use periodic wrapping or nearest physical-cell
extrapolation; they do not use state momentum reflection or prescribed inflow
values. Open geometry cannot use a covered source donor. Existing boundary
advance accepts optional conserved rates as well. All public stage outputs
remain unchanged when validation or source conversion fails.

This source extrapolation follows the existing regular SDC approach; full
pinned boundary-source parity has not been established. The coupled EB SDC
loop and redistribution remain to be connected.

## Focused evidence

Artifacts: `build/sdc-eb-boundary-stage-checks/` (local). GNU Fortran 13.3.0,
existing Ninja Debug tests-enabled build `build/spray-amr-reactive-debug`;
compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- CTest **1/1 passed**, 35.26 s (35.31 s total): `env -u LD_LIBRARY_PATH
  -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R
  '^regression_reactive_eb_hydro_3d$' --output-on-failure` (`ctest.log`).
- Added stage checks for a nonuniform regular 4-cubed physical grid with
  extended geometry. Three axes times five boundary kinds: periodic, slip,
  no-slip, isothermal wall and prescribed inflow; other faces outflow. With
  zero source, state plus dt times rate reproduces the existing update exactly;
  all integrated fluxes and boundary/wall inventories also match exactly.
- Nonzero extrapolated momentum source under outflow matches the source-aware
  boundary advance exactly. NaN source in an active cell rejects without
  changing sentinel rates/inventories or previously returned fluxes. Existing
  cut-cell, redistribution and boundary tests also passed.
- Python syntax for `tools/check_project_contract.py`, project contract,
  touched Fortran line lengths and `git diff --check` passed.

No coupled EB SDC, application trajectory, MPI path, selected build, full
matrix, installed audit or physical validation was run. The full port remains
incomplete.
