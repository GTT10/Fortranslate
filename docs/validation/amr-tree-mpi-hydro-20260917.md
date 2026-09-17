# Synchronous distributed tree hydro — 2026-09-17

## Source and scope

Base: `16f8244cfc50738f6cd5a8539e055b9008d854ec`, local `main`.
Connects owner-only state/temperature storage to common-step SSPRK2 hydro through
distributed halos, peer/parent flux correction, state restriction and EOS recovery.
The serial arithmetic and first-stage final temperature guess are retained.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_hydro_3d_mod.F90`:
  `c028ca5ef429db859fbd884080b16d9226702a3605a7b2c5cd6db546e5ac673e`
- `tests/unit/test_mpi_amr_patch_tree_hydro_3d.F90`:
  `03655209d134c3195eb710bcb7817b90ff327659306f08be5bc851ee4a568fab`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, incremental build in the existing clean
system-MPI configuration; compiler/backend options unchanged.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_hydro_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_hydro_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final build succeeded. Focused CTest **2/2 passed**, 6.85 s total (2 ranks 3.48 s,
3 ranks 3.34 s), no failures/skips. Provenance, contract and whitespace checks passed.
No repository Python changed; syntax checks were not applicable. An initial build
failed on an overlong generated Fortran source line; the line was wrapped before
successful build/testing. Exact-control comparisons emit real-comparison warnings.

Fixtures use elementary H2/O2, spatially varying conserved states, three levels
with ratios 2/3, distributed/all-on-one ownership and empty ranks. Four modes
exercise PCM/characteristic PLM (MC), periodic and physical outflow/slip walls,
and prescribed inflow. After one common step, every owned conserved/temperature
array and the global boundary budget/maximum reflux correction match the serial
API exactly. No nonowner state or temperature arrays are allocated. Invalid flux
solver and mismatched timestep cases preserve prior state, temperature and diagnostics.

Logs: `build/amr-tree-mpi-hydro-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_hydro_3d`.

## Remaining scope

This is hydro operator parity, not a complete MPI reacting/transport application.
Chemistry/transport composition, per-level time endpoints/subcycling, rebalancing,
selected-MPI tree runtime and application integration remain. The caller supplies
replicated thermodynamic records and an appropriate dt. There is no internal retry
or production scaling qualification. No external solver parity or physical
validation is claimed.
