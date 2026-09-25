# Regular MPI spray hydro inflow masks — 2026-09-25

Source: base `a59c1ab`, implementation `2f9eab6`. GNU Fortran 13.3.0,
OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts:
`build/ppm-regular-mpi-inflow-checks/`.

The regular MPI spray hydro PPM branch passes physical inflow masks to final
PLM flux replacement. X masks are active only on rank zero / the last rank,
which own the global domain ends in this slab decomposition. Y/Z physical
faces span all ranks. Periodic/inter-rank faces retain CTU fluxes. The shared
CTU kernel applies replacement only with the `pelec` pressure policy; `hllc`
retains its existing behavior. Particle stepping and feedback are unchanged.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_ppm --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_ppm_4$' --output-on-failure`:
  **1/1 passed**, 5.15 s test / 5.19 s total; no failures/skips.
- The fixture now covers both `hllc` and `pelec`, periodic and physical-wall
  cases, and each of six Mach-inflow orientations. Four ranks own two X cells
  each, exercising four-layer exchange through short slabs. Every case matches
  serial conserved state and temperature exactly; boundary inventories meet the
  existing `2e-13` relative/absolute-scaled check. Excessive-CFL and mismatched
  reconstruction cases preserve local state and reject collectively.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

This checks the gas-hydro path used by the spray runtime, not particle coupling
or nozzle physics. Selected build/execution, application restart, full matrix,
compiled full-CTU comparison, convergence and physical validation were not run.
Compiler warnings remain. The overall port is still incomplete.
