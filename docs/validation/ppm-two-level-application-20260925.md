# Serial two-level AMR PPM application — 2026-09-25

Source: base `91f8779`, implementation `3d63707`. GNU Fortran 13.3.0,
Open MPI 4.1.6, Ninja Debug, Python 3.12.3. Reused fixed tests-enabled
`build/spray-amr-reactive-debug` and selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug` separately.
Artifacts: `build/ppm-two-level-app-checks/` (exact commands in `commands.txt`); fixed fixture:
`build/spray-amr-reactive-debug/tests/amr_two_level_ppm_restart_3d/run-9ky77zb3`.

Both serial front ends accept `characteristic_ppm`. The transport split now
forwards this selector to hydro. Regridding uses the existing conservative
limited-linear prolongation for PPM, and the checkpoint configuration validator
accepts it while retaining exact stored-selector matching. No schema change or
numerical tolerance change was needed. Two-level MPI input remains disabled
for PPM; it needs a separate distributed hydro implementation.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_amr_reactive_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_amr_reactive_3d_selected --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^amr_two_level_ppm_restart_3d$' --output-on-failure`:
  **1/1 passed**, 48.48 s test / 48.53 s total; no failures or skips.
- The 4x4x4 full H2/O2 fixture uses chemistry, viscosity, heat/species transport,
  barodiffusion, a moving isothermal wall and prescribed-pressure outlet.
  Step-zero regridding expands the coarse patch to `1:4,2:3,2:3`. It completes
  two coarse steps and 32 transport fine substeps to `1e-10` s. Coarse/fine CSV
  and the boundary ledger match uninterrupted execution exactly after restart.
  Changed reconstruction, outlet pressure, wall temperature, chemistry tolerance,
  integrator, rates and mechanism identity reject restart. The chemistry-disabled
  control differs from the reacting result.
- The selected input in `selected-final/reference.nml` preserves the fixed
  fixture and replaces the mechanism selector and composition syntax. From that
  directory, with the library-path variables unset:
  `../../spray-amr-reactive-selected-debug/pelef_amr_reactive_3d_selected reference.nml`:
  passed. Comparing its two CSVs with the fixed reference using
  `tools/compare_amr_reactive_3d_restart.py` passed byte-exactly:
  coarse SHA256 `838fac6aa9bb912d7ec8607373d6c4bc85473bebc7475e8c4dc73479a9d59a79`,
  fine SHA256 `4be2828ed5b407c2d04bcf8823c77bd5fc608ecd31bc53928739d24dbc28a4f4`.
  The initial selected fixture incorrectly removed `x_lower`/`x_upper` with the
  fixed composition entries, causing a comparison mismatch. The input conversion
  was narrowed to species names; the initial artifacts remain under `selected/`.
- `python3 -m py_compile tools/check_amr_physical_boundary_restart.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected restart, MPI two-level PPM, full matrix, convergence, external solver
parity and physical validation were not run. Compiler warnings remain. Reference
primitive-pressure tracing fidelity remains open.
