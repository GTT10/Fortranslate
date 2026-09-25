# Compiled pinned PPM scalar kernels — 2026-09-25

Source: base `ae5fe9f`, comparison tool `a4e2bc4`. GNU Fortran/C++ 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-pinned-kernel-checks/`.

`tools/check_pinned_ppm_kernels.py` reads the manifest-pinned PeleC
`Source/PPM.H` with `git show`, extracts the unchanged `ppm_reconstruct` and
`ppm_int_profile` function bodies, and compiles a scalar C++ driver. Its shim
only supplies `amrex::Real=double`, standard min/max and stencil index constants.
It does not translate the reference arithmetic into Python or Fortran.
The existing Fortran unit executable exposes an input/output stream mode that
calls the production reconstruction and profile-integration routines directly.
Normal execution still runs its original unit assertions.

Reference commit `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PPM.H SHA256 `294cb223a95bf569b57700cbd9083ae77e562b7124f9cdda6aa338eada890e93`.
The header, compiled-source driver, both output tables, common input and JSON
report are retained under `oracle/` in the artifact directory. The report also
records the compared Fortran executable hash and exact C++ compile command.

## Completed checks

- Provenance preflight:
  `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm --parallel 2`: passed.
- With those library variables unset:
  `python3 tools/check_pinned_ppm_kernels.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm --work-dir build/ppm-pinned-kernel-checks/oracle`: passed.
  **675 cases / 5,400 values**, maximum stencil-scaled absolute error
  **1.2908639827084056e-16**, fixed acceptance limit `3e-14`.
  Inputs cover constant/linear/decreasing/step/extremum/random stencils at
  scales `1e-8`, `1`, `1e6`, flattening `0`, `.37`, `1`, velocity `-2` through
  `2`, sound speed `1` and `dt/dx=.2`. Random seed is `20260925`.
  The eight observables are two reconstructed edges and six wave integrals.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm$' --output-on-failure`:
  **1/1 passed**, 0.00 s test / 0.05 s total; no failures/skips.
- `python3 -m py_compile tools/check_pinned_ppm_kernels.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

This is compiled external scalar-kernel evidence, not a full PeleC executable
comparison. It excludes characteristic EOS projection, Riemann solution,
transverse CTU, boundary treatment and coupled evolution. No numerical changes,
full matrix, convergence or physical validation were needed or performed.
