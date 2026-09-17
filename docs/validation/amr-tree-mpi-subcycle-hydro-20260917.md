# MPI recursive hydro subcycling — 2026-09-17

## Source and scope

Base: `b094a9f531b667664d3b0caeca23c696ff4936d7`, local `main`.
Adds owner-local recursive hydro interval execution, coarse-time endpoints,
substep face accumulation, depth-limited distributed reflux and coarse CFL
selection. Fixed/selected runtime libraries include the API. Serial algorithms,
mechanism sources, tolerances and golden hashes remain unchanged.

Verified SHA-256:

- `src/parallel/mpi_amr_patch_tree_subcycle_3d_mod.F90`:
  `420b4c6ce00ff642c209a3bf38bbeb599c20e9c009fa9e05284b5068b28ef0ed`
- `tests/unit/test_mpi_amr_patch_tree_subcycle_3d.F90`:
  `1601c7e844bed21f67125ad00223dfa3fdb2be8ca18d1402c3f9f1f96879a05d`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target test_mpi_amr_patch_tree_subcycle_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_subcycle_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed/selected builds passed. CTest **2/2 passed**, 20.27 s total: 2 ranks
10.07 s, 3 ranks 10.16 s; no failures/skips. Provenance, project contract and
whitespace checks passed. No Python sources changed. Existing compiler warnings
remain.

The nonuniform H2/O2 three-level fixture has refinement ratios 2 then 3. Four
modes exercise PCM/characteristic PLM, periodic boundaries, outflow/slip walls,
reservoir inflow, split owners and all payload on one rank with empty owners.
The MPI scaled coarse CFL equals serial exactly. One coarse interval at 0.1 of
that bound executes exactly **1/2/6** steps by depth on every rank.

Owned conserved states and temperatures, physical boundary-change vectors,
maximum reflux corrections and executed counters match serial exactly. Invalid
solver and rank-dependent dt are collectively rejected while preserving prior
state, temperatures, diagnostics and counters. This focused test does not inject
a failure specifically during a later child interval.

Logs: `build/amr-tree-mpi-subcycle-checks/{runtime-build,build,ctest}.log`.
Binary: `build/selected-tree-mpi-debug/tests/test_mpi_amr_patch_tree_subcycle_3d`.

## Remaining scope

MPI recursive transport, reaction/transport/hydro recursive composition and
application/checkpoint mode selection are not connected yet. Advance does not
compute/enforce CFL or retry dt. Coefficient identity is the caller's contract.
No global field gather is used, but no production scaling study, full matrix,
installed-prefix audit, external solver parity or physical validation was run.
