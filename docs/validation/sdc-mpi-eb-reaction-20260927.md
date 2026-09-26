# MPI EB tree SDC reaction stage — 2026-09-27

## Source and behavior

Implementation `4b3e1168`, based on `b352779a` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`react_sdc_mpi_eb_tree_3d` evaluates forced EB chemistry only on owned patches.
It collectively agrees timestep/tolerances, species/reaction counts and level
selection, validates owner-local field layouts, and stages state, temperature
and reaction estimates privately. The existing EB reaction kernel handles
fluid/solid masking and refined parents. A collective success decision precedes
publication on every rank. Optional `target_depth` retains inactive levels.
No gather, average-down or reflux is performed by this reaction stage.

## Completed focused checks

GNU Fortran 13.3.0, configured `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_mpi_eb_tree_full_3d` built successfully.
- CTest `unit_mpi_eb_tree_full_3d_2`: 1/1 passed, 51.08 s test / 51.13 s total.
- New checks use three levels, five patches and two ranks, with round-robin
  owners and all patches on rank 0 (rank 1 idle). Forced SDC state, temperature
  and reaction estimates match the serial EB tree reaction exactly. Non-owned
  state remains unallocated.
- Depth-2 updates match serial state/reaction outputs and retain inactive
  state, temperature and reaction sentinels.
- A negative density on the last patch's owner rejects globally. Separately,
  a negative tolerance on rank 1, including when idle, rejects globally.
  State, temperature and reaction outputs remain unchanged on both ranks.
- Existing MPI full split-step tests remain in the same executable.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_full_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-mpi-eb-reaction-checks/`.
No coupled MPI EB SDC hydro/transport, application restart, selected build,
full matrix, installed-binary audit, scaling or physical validation was run.
