# EB SDC reaction stage — 2026-09-27

Implementation `3acfa4ea`, base `61fc10d6`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`react_sdc_eb_state_3d` in `reactive_eb_3d_driver_mod` accepts old state,
temperature, EB geometry and nonreacting rates. It calls the existing forced
SDC reaction kernel on each fluid/cut cell, retaining covered old state and
temperature and the caller's covered reaction estimate. Rates and state are
per fluid volume; the local ODE does not multiply rates by volume fraction.
Geometry, dimensions, mechanism presence and tolerances are validated first.
Results remain private until every active cell succeeds. This is a reaction
stage, not a complete EB SDC timestep.

## Focused evidence

Artifacts: `build/sdc-eb-reaction-checks/` (local). GNU Fortran 13.3.0,
existing Ninja Debug tests-enabled build `build/spray-amr-reactive-debug`;
compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_reactive_eb_chemistry_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- CTest **1/1 passed**, 1.52 s (1.57 s total): `env -u LD_LIBRARY_PATH
  -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R
  '^unit_reactive_eb_chemistry_3d$' --output-on-failure` (`ctest.log`).
- Added a four-cell plane-cut geometry containing covered, cut and regular
  cells. Momentum/energy forcing, dt=1e-9 s, elementary H2/O2 reaction kernel.
  Fluid/cut state, temperature and reaction estimates match the regular-grid
  forced reaction exactly. Covered state and rates contain NaN sentinels;
  covered old state/temperature and reaction estimates are retained exactly.
- A negative-density forcing in the last active cell causes failure after
  earlier active cells were processed. All output arrays remain bitwise equal
  to their pre-call values, including covered sentinels. Existing masked
  chemistry and Strang tests also passed.
- Python syntax for `tools/check_project_contract.py`, project contract and
  `git diff --check` passed.

No coupled EB SDC, hydro/transport SDC stage, application run, MPI path,
selected build, full matrix, installed audit or physical validation was run.
The cellwise adapter reuses the reaction kernel but has not been performance
qualified. The full port remains incomplete.
