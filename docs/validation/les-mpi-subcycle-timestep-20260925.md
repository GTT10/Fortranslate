# MPI LES subcycle timestep selection — 2026-09-25

Implementation `a881d0d`, based on `ce6ec99`. `compute_mpi_amr_tree_les_timestep_3d` accepts optional
`subcycle=.true.`. Each owner evaluates its existing model/ghost diffusion
bound, divided by the cumulative squared refinement ratio for that depth.
The global minimum supplies the coarse interval. Default/false retains common
time stepping; no full-state gathering is added.

Effective subcycling joins the collective control agreement before halo work.
Rank disagreement rejects without changing dt or diffusivity. Count allocation
and int64 overflow validation are collective, including ranks with no owned
nodes. The existing maximum diffusivity reduction remains unchanged.

## Scope and checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled runtime/test
and selected H2O2 tests-disabled runtime built. Selected execution was not run.
CTest `spray_les_mpi_amr_tree_2`: **1/1 passed**, 118.64 seconds, total 118.68.
No failures, skips or tolerance changes occurred. Provenance preflight, Python
syntax of the contract checker, project contract and whitespace passed.
The existing two-rank AMR LES regression adds periodic root/nested subcycle
selector comparisons for all four models, against the serial selector at
relative tolerance `2e-12`. Both concentrated ownership (empty rank) and
distributed ownership are exercised. A true/false rank disagreement checks
unchanged output sentinels.

Existing common-selector, state, boundary-budget and rollback cases remain in
the same test, including inflow, pressure outflow and moving/isothermal walls.
Subcycle-mode physical-boundary comparisons were not added in this check.
This enables selection only; MPI level stepping, recursive scheduling and
coupled application dispatch remain separate implementation work. No external
solver parity, physical validation or full configuration matrix is claimed.

## Commands and artifacts

Logs: `build/les-mpi-subcycle-timestep-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_amr_tree_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
