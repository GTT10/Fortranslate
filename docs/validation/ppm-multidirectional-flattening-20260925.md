# Multidirectional PPM flattening — 2026-09-25

Source: base `95f5088`, implementation `5e9c58f`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts:
`build/ppm-multidirectional-flattening-checks/`.

The optional pressure-traced regular-cell CTU path now computes one flattening
coefficient per donor cell as the minimum over all three spatial directions,
and applies it to every normal reconstruction. The normal cell reconstruction
accepts an optional validated coefficient; callers omitting it retain their
existing directional behavior. The supplied four-layer halo union already
contains all necessary pencils, so unused far-corner padding remains ignored.

## Reference and corrected remaining scope

Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/PPM.cpp::trace_ppm`, inspected via `git show` in
`/mnt/sn850/Projects/gtt2/external/PeleC`. SHA256:
`312dd0ae66bd56b62ec0c3de8d7690ce8181ac8098124188d6576e8a79eab703`;
the inspected file is retained with artifacts.

The reference explicitly minimizes flattening over spatial axes. It computes
internal-energy characteristic amplitudes, but the assignments using those
amplitudes are commented out. Active code instead calls `RYP2E` using traced
density, composition and pressure and assigns `rho*e` on both normal faces.
PeleF's primitive-to-conserved normal predictor follows this active policy.
The earlier statement that independent normal internal-energy tracing was an
additional required reference step was incorrect. This source inspection does
not prove full numerical parity: EOS admissibility, species normalization,
floating-point ordering and complete CTU behavior still require bounded
external comparison. Application binding of the new pressure policy remains open.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure`:
  **1/1 passed**, 4.81 s test / 4.85 s total, no failures/skips.
- Added a compressive transverse shock with a smooth normal density gradient:
  directional reconstruction preserves the gradient, while the transverse
  coefficient collapses both normal states exactly to the cell center. Invalid
  coefficient leaves outputs zero. New-policy periodic and interior-subdomain
  exact flux parity now includes a transverse pressure/velocity discontinuity.
  Existing pressure-work, rollback and six-face boundary inventory checks remain.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected runtime, MPI, full matrix, compiled external oracle, convergence and
physical validation were not run. Compiler warnings remain.
