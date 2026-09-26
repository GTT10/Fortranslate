# EB SDC transport source stage — 2026-09-27

Implementation `a32e2e81`, base `9477409f`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`compute_sdc_eb_transport_3d` returns conservative transport rates, time-integrated
face flux densities and the species limiter without advancing state. It reuses
the existing EB flux/RHS calculation, including optional ghost donors, physical
faces and wall temperature. Rates include open-area divergence divided by fluid
cell volume. Exported flux densities include dt; callers apply aperture and
face area in inventories/reflux, matching `eb_integrated_flux_3d` conventions.
Covered rates are zero. Outputs publish only after successful finite checks.
Small-cell redistribution is not performed by this source stage and remains
part of future coupled EB SDC work.

## Focused evidence

Artifacts: `build/sdc-eb-transport-checks/` (local). GNU Fortran 13.3.0,
existing Ninja Debug tests-enabled build `build/spray-amr-reactive-debug`;
compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_reactive_eb_transport_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- CTest **1/1 passed**, 0.99 s (1.04 s total): `env -u LD_LIBRARY_PATH
  -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R
  '^unit_reactive_eb_transport_3d$' --output-on-failure` (`ctest.log`).
- Extended nonuniform 4-cubed plane-cut cases for x/y/z with minimum fluid
  volume fraction 0.05, all transport flags enabled. Rates and integrated
  flux densities match the existing EB flux calculation exactly; state and
  temperature remain bitwise unchanged. Closed-domain volume-weighted rate
  inventories are zero within 5e-12 times their absolute integral scale.
- Invalid transport diameter rejects the stage while retaining rates, all
  flux arrays and limiter. Existing regular, small-cell, curved-geometry and
  full-driver checks in the executable also passed.
- Python syntax for `tools/check_project_contract.py`, project contract and
  `git diff --check` passed. After testing only a flux-convention comment and
  a test label changed; no numerical code changed.

No coupled EB SDC, source-aware hydro, application run, MPI path, selected build,
full matrix, installed audit or physical validation was run. The full port
remains incomplete.
