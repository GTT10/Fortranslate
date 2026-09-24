# MPI common-step AMR PPM — 2026-09-25

Source: base `0973dfd`, implementation `8814d4e`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3. Reused fixed
MPI tests-enabled `build/spray-amr-reactive-debug`. Artifacts:
`build/ppm-mpi-amr-common-checks/`.

`advance_mpi_amr_tree_hydro_3d` accepts collectively agreed
`characteristic_ppm`. It requests four-layer owner-local ghosts with the
existing limited coarse interpolation and physical/Mach/pressure boundaries,
then computes each owned patch's CTU fluxes. Existing collective reflux,
restriction and EOS recovery complete one CTU step before committing state
and diagnostics. PCM/PLM keep their two-stage update. Empty ranks remain
participants in collective validation, ghost exchange and synchronization.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_ppm_3d_2$' --output-on-failure`:
  **1/1 passed**, 6.54 s test / 6.59 s total, two MPI ranks.
- Six cases on root/child/grandchild with refinement ratios 2 and 3 cover
  periodic, outflow/slip wall, reservoir inlet, pressure outlets and Mach
  inlets with extrapolated or prescribed pressure. Ownership includes a rank
  with no patches. Owned U/T, boundary budget and maximum reflux correction
  agree exactly with the serial PPM implementation; nonowned fields remain
  unallocated. Bad Riemann solver, mismatched timestep/pressure/Mach and a
  rank-local invalid wall pressure reject collectively and preserve fields
  and diagnostics. Existing parity and rollback checks were retained.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No failed tests occurred; compiler warnings remain. Three-rank execution,
PCM/PLM variants, selected builds, recursive MPI PPM, application input,
full matrix, convergence and physical/external validation were not run here.
Reference primitive-pressure parity remains open.
