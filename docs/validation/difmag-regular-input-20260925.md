# Regular gas artificial-viscosity input — 2026-09-25

Source: base `22f50b3`, implementation `2bf0a91`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug` and selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug`. Artifacts:
`build/difmag-regular-input-checks/`.

Regular fixed/selected gas input accepts `difmag`: -1 means automatic (0.1 for
PPM/pelec, zero otherwise), zero disables viscosity and positive values override
PPM strength. Configuration rejects nonfinite values, other negative values
and explicit coefficients with PCM/PLM. A reader capability flag rejects explicit
coefficients in application families whose runtime/restart path is not wired.
The regular application prints the resolved coefficient and forwards it through
full splitting, boundary hydro or periodic Euler stepping to CTU. New optional
library arguments preserve existing call sites; full/boundary stepping rejects
nonzero viscosity on unsupported reconstructions.

These regular gas executables have no checkpoint/restart interface. No restart
format was changed. AMR/spray input and restart binding remain unfinished;
the capability flag prevents silently accepting an ignored coefficient there.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d test_reactive_ppm_step_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^(reactive_3d_difmag_input|unit_reactive_ppm_step_3d)$' --output-on-failure`:
  **2/2 passed**, 3.05 s step / 12.88 s application, 15.97 s total; no failures/skips.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset:
  `python3 tools/check_reactive_difmag_input.py --selected --executable build/spray-amr-reactive-selected-debug/pelef_reactive_3d_selected --work-dir build/difmag-regular-input-checks/selected`: passed.
- Both application variants use the 4-cubed Mach-inlet hotspot/transport fixture.
  Automatic, explicit -1 and explicit 0.1 yield identical CSVs within each
  variant; zero produces different CSVs. Negative -2 and explicit coefficient
  with PCM reject before CSV publication. The printed coefficient is checked.
- Fixed fixture:
  `build/spray-amr-reactive-debug/tests/reactive_3d_difmag_input/run-pzgfdcns/`;
  selected fixture: `build/difmag-regular-input-checks/selected/run-y3wi401v/`.
- `python3 -m py_compile tools/check_reactive_difmag_input.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No cross-variant bitwise comparison, AMR/spray application execution, MPI,
full matrix, external full-CTU parity, convergence or physical validation was
run. Compiler warnings remain.
