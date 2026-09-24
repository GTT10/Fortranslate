# Regular MPI spray PPM application — 2026-09-25

Source: base `04cd602`, implementation `203f6c6`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3. Reused fixed
MPI tests-enabled `build/spray-amr-reactive-debug` and selected H2/O2
MPI tests-disabled `build/spray-amr-reactive-selected-debug` separately.
Artifacts: `build/ppm-regular-mpi-app-checks/`; fixture:
`build/spray-amr-reactive-debug/tests/spray_mpi_ppm_restart/run-vftwwjpl`.

`advance_mpi_mixture_reaction` now collectively validates and forwards an
optional hydro reconstruction selector before committing particle/gas work.
The fixed/selected regular MPI application accepts
`hydro_reconstruction='characteristic_ppm'` in `&mpi_spray`; default PLM is
unchanged. PPM adds `mpi-hydro-characteristic-ppm-v1` to the checkpoint context.
The default PLM context does not receive the new marker, preserving its prior
identity. Both ordinary and impact/film coupling branches forward the selector.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2`: passed after repair.
  Initial builds rejected an overlong Fortran controls-array line. It was
  split without semantic change; initial and final build logs are retained.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_ppm_restart$' --output-on-failure`:
  **1/1 passed**, 66.45 s test / 66.49 s total.
- The existing MPI application checker runs with PPM on the 8x4x4 mixed-liquid
  fixture, implicit chemistry and molecular transport. Same-rank continuation
  agrees exactly in parsed gas/particle records; two-to-four-rank continuation
  passes its existing tolerances. Changing PPM to PLM rejects restart. Existing
  output protection, property/context rejection, scheduled injection, LES,
  breakup and physical-boundary/film variants also ran through this checker.
  No tolerances changed.
- `python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No runtime test failures occurred. Compiler warnings remain. Selected runtime,
full matrix, historical PLM checkpoint replay, convergence and external/physical
validation were not run. Reference primitive-pressure fidelity remains open;
these numerical continuation checks do not establish nozzle/diesel validation.
