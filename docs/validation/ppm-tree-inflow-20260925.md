# Patch-tree PPM inflow masks — 2026-09-25

Source: base `eed5009`, implementation `b75137e`. GNU Fortran 13.3.0,
OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-tree-inflow-checks/`.

The shared tree ghost module now derives each node's physical inflow mask
from the integer layout: zero origin for lower faces, origin plus node shape
equal to the level-domain extent for upper faces. A halo crossing the physical
boundary does not mark an interior node face. Serial and MPI common-step and
single-level CTU callers pass that mask to final PLM inflow flux replacement;
recursive stepping uses the same single-level path. MPI failure still flows
through collective agreement. No particle or EB policy was changed.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d test_mpi_amr_patch_tree_subcycle_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_ppm(_subcycle)?_3d_2$' --output-on-failure`:
  **2/2 passed**, common-step 4.24 s, recursive 13.11 s, total 17.39 s;
  no failures/skips.
- PPM variants now use the `pelec` solver so pressure tracing and inflow
  replacement are exercised. Other variants retain `hllc`. Added direct geometry
  checks: root touches all six faces, child only the three lower faces, interior
  grandchild none. Existing periodic/physical, prescribed inflow, Mach inlet and
  pressure-outlet cases compare exact owner-local state/temperature, boundary
  inventory and reflux diagnostics with serial; recursive counts are exactly
  1/2/6. Failure/collective-context checks remain.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No selected build/execution, application restart, regular MPI spray, full
matrix, external full-CTU oracle, convergence or physical validation was run.
Compiler warnings remain. Regular MPI spray mask propagation remains open.
