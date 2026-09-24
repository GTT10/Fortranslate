# Sparse two-level MPI PPM application — 2026-09-25

Source: hierarchy `172f29a`, application implementation `8c43459`.
GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug,
Python 3.12.3. Reused fixed tests-enabled `build/spray-amr-reactive-debug`
and selected H2/O2 tests-disabled `build/spray-amr-reactive-selected-debug`
separately. Artifacts: `build/ppm-mpi-two-level-app-checks/`; fixture:
`build/spray-amr-reactive-debug/tests/mpi_amr_two_level_ppm_restart_3d/run-2lmcz_57`.

Fixed/selected MPI two-level front ends now accept `characteristic_ppm`, including
transport configuration. Sparse regridding uses the existing conservative
limited-linear prolongation, matching serial PPM remapping. Existing checkpoint
reconstruction matching and rank-neutral state redistribution are reused;
no schema change or weakened tolerance was needed.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_reactive_3d_selected --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^mpi_amr_two_level_ppm_restart_3d$' --output-on-failure`:
  **1/1 passed**, 44.22 s test / 44.26 s total; no failures or skips.
- The 4x4x4 fixed full-H2/O2 fixture uses chemistry, viscosity, heat/species
  transport, barodiffusion, a moving isothermal wall and prescribed-pressure
  outlet. Step-zero regridding expands the patch to `1:4,2:3,2:3`. A four-rank
  run checkpoints after the first coarse step; two ranks finish the second step
  at `1e-10` s. Both final CSVs equal the uninterrupted serial reference byte for
  byte; the boundary ledger passes its existing reduction tolerance.
  Coarse SHA256: `838fac6aa9bb912d7ec8607373d6c4bc85473bebc7475e8c4dc73479a9d59a79`.
  Fine SHA256: `4be2828ed5b407c2d04bcf8823c77bd5fc608ecd31bc53928739d24dbc28a4f4`.
- PPM-to-PLM, outlet pressure, wall temperature, chemistry tolerance, integrator,
  reaction-rate and mechanism-identity changes reject restart. The chemistry-off
  control differs from the reacting result.
- `python3 -m py_compile tools/check_amr_physical_boundary_restart.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected MPI runtime/restart, full matrix, other refinement ratios, convergence,
external solver parity and physical validation were not run. Compiler warnings
remain. Reference primitive-pressure tracing fidelity remains open.
