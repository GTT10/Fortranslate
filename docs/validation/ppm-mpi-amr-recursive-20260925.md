# MPI recursive AMR PPM — 2026-09-25

Source: base `769223c`, implementation `9b19d8f`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3. Reused fixed
MPI tests-enabled `build/spray-amr-reactive-debug`. Artifacts:
`build/ppm-mpi-amr-recursive-checks/`.

`advance_mpi_amr_tree_hydro_level_3d` accepts collectively agreed
`characteristic_ppm`. It exchanges four-layer owner-local ghosts, retaining
limited coarse interpolation and starting-time coarse fractions, then
computes one CTU step and dt-weighted fluxes. PCM/PLM retain their two stages
and half-dt weights. Fluxes undergo the collective finite check before
publication. The existing recursive driver accumulates these integrated
fluxes and performs reflux/restriction at synchronization; no recursive
schedule change was needed.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_subcycle_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_ppm_subcycle_3d_2$' --output-on-failure`:
  **1/1 passed**, 19.69 s test / 19.74 s total, two ranks.
- Six cases use root/child/grandchild with refinement ratios 2 and 3,
  including an empty rank. Periodic, outflow/slip wall, reservoir inlet,
  pressure outlets and Mach inlets with extrapolated/prescribed pressure
  are covered. Owned U/T, boundary budget, maximum reflux correction and
  1/2/6 level-step counts agree exactly with serial recursive PPM. Nonowned
  fields remain unallocated. Bad solver, mismatched timestep/pressure/Mach
  and rank-local invalid wall pressure reject collectively with rollback.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No failed tests occurred; existing compiler warnings remain. Three-rank,
selected builds, PCM/PLM test variants, reaction/transport composition,
AMR application input/restart, full matrix and external/physical/convergence
validation were not run in this increment. The conservative transverse EOS
correction still differs from the reference primitive-pressure equation.
