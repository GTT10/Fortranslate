# EB SDC hydro rate and flux stage — 2026-09-27

Implementation `f668d017`, base `a090e57e`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`compute_sdc_eb_hydro_3d` uses the extended geometry/ghost layout of the existing
CTU advance and accepts conserved predictor-source rates on its working grid.
It shares source-aware CTU, artificial viscosity, species flux normalization,
face-centroid interpolation, physical geometry extraction, conservative flux
divergence and boundary/wall inventory accounting. It returns rates on the
physical grid, time-integrated flux densities and separate domain boundary and
EB wall inventories. No state update, temperature recovery or redistribution
is performed by the source stage. Rates precede redistribution and already
include the EB wall contribution. All outputs publish after successful checks.
The ordinary CTU advance continues through its existing update when no rate
output is requested.

## Focused evidence

Artifacts: `build/sdc-eb-hydro-stage-checks/` (local). GNU Fortran 13.3.0,
existing Ninja Debug tests-enabled build `build/spray-amr-reactive-debug`;
compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- CTest **1/1 passed**, 25.29 s (25.33 s total): `env -u LD_LIBRARY_PATH
  -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R
  '^regression_reactive_eb_hydro_3d$' --output-on-failure` (`ctest.log`).
- Nonuniform regular and x/y/z cut geometries: 8-cubed working grids cropped
  to 2-cubed physical grids, dt=1e-8 s. With zero predictor source, old state
  plus dt times returned rate equals the existing undistributed Euler update
  exactly. All integrated fluxes and boundary/wall inventories match exactly.
- Volume-weighted rates integrated over dt equal boundary plus wall inventory
  within 2e-12 times the existing inventory scale. Negative artificial
  viscosity rejects without changing rates, fluxes or inventories. Existing
  source-aware CTU and redistribution regression checks also passed.
- Python syntax for `tools/check_project_contract.py`, project contract,
  touched Fortran line lengths and `git diff --check` passed.

No coupled EB SDC, application run, MPI path, selected build, full matrix,
installed audit or physical validation was run. Redistribution must still be
connected to the coupled SDC update. The full port remains incomplete.
