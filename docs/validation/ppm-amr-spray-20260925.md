# MPI AMR spray PPM selection — 2026-09-25

Source: base `f1a8d11`, implementation `4c0fa2f`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3. Reused fixed
MPI tests-enabled `build/spray-amr-reactive-debug` and selected H2/O2
MPI tests-disabled `build/spray-amr-reactive-selected-debug` separately.
Artifacts: `build/ppm-amr-spray-checks/`; generated fixture:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_ppm_restart/restart_nsho5od2`.

Both MPI AMR spray front ends now accept `characteristic_ppm`.
The existing particle/gas split forwards reconstruction to the newly connected
MPI AMR gas split. Mach inlet timestep sampling now selects four-layer limited
coarse donors for PPM. Existing checkpoint context already includes the gas
configuration; a new checker variant verifies changed reconstruction rejection.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_ppm_restart$' --output-on-failure`:
  **1/1 passed**, 6.04 s test / 6.08 s total.
- Fixed periodic mixed-liquid fixture runs recursive PPM with scheduled
  injection and a topology change. Two-rank uninterrupted output and
  two-to-four-rank continuation agree in gas/particle values within the
  existing 2e-12 relative / 1e-25 absolute tolerances. Three injection events
  produce 12 final parcels. The checker also covers changed PPM-to-PLM
  selection, changed subcycling and injection settings, and its existing
  checkpoint corruption/overwrite rejection cases. No tolerances changed.
- `python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No failed tests occurred; compiler warnings remain. Selected runtime,
chemistry/molecular transport/LES in this spray fixture, physical/Mach boundaries,
common-step spray execution, full matrix, convergence and physical validation
were not run. The regular-grid spray and two-level AMR front ends remain
separate PPM integration work. This is not nozzle or diesel-spray validation.
