# MPI owner-local SDC reaction stage — 2026-09-26

Implementation `d205c3e6`, base `0e9cd273`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Behavior

`react_sdc_mpi_amr_tree_3d` integrates the regular forced SDC reactor over owned
patches, including covered parents. It gathers no full state. Non-owned fields
remain unallocated. Topology/ownership/layout and time/tolerance/species-count/
reaction-count/level controls are checked collectively. This follows the
existing MPI mechanism-consistency contract; full mechanism contents are not
hashed by this stage.

Optional level selection preserves inactive output state, temperatures and
reaction estimates, requiring valid initial output trees. Whole-tree mode
allows initially unallocated outputs. All numerical outputs publish only after
all ranks succeed. Restriction and temperature synchronization across levels
remain the responsibility of the eventual coupled SDC driver.

## Completed focused checks

Artifacts: `build/sdc-mpi-reaction-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`. Configured `/usr/bin/mpifort` and
`/usr/bin/mpiexec` report Open MPI 4.1.6. The shell's unqualified `mpiexec`
resolves to a different Hydra launcher; the test uses CMake's absolute Open MPI
path. No compiler or MPI configuration was changed.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_chemistry_3d
  --parallel 2` (`build.log`). Warnings remain.
- New focused CTest `unit_mpi_amr_patch_tree_sdc_reaction_3d_4`: **1/1 passed**,
  57.91 s (57.96 s total), `ctest.log`. It invokes `--sdc-only` and does not
  rerun the existing Strang/pressure test branches.
- Four ranks, three-level tree, ownership `[0,1,0]`; ranks 2 and 3 own no cells.
  Nonuniform state with external nonreacting rates yields exact serial/MPI
  conserved state, temperature and reaction estimates. Non-owned arrays remain
  absent. Level-selective integration retains deliberately changed inactive
  fields. A late owner-local invalid density rejects collectively and preserves
  all outputs. An idle-rank time mismatch is also rejected.
- Command: `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug
  -R '^unit_mpi_amr_patch_tree_sdc_reaction_3d_4$' --output-on-failure`.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and `git diff --check` passed.

This is a distributed reaction stage, not the complete MPI SDC integrator.
No MPI SDC hydro/transport coupling, application run, selected build, restart,
rank-changing continuation, full configuration matrix, install audit, pinned
full-solver comparison or physical validation was performed.
