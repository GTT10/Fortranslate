# Distributed tree temperature tagging — 2026-09-17

## Source and scope

Base: `2f1379aef5060f85d59d5577c72b0b5da09135d2`, local `main`.
Adds owner-local temperature-tag planning, hot-cell box/temperature exchange and
collective candidate publication. Extracts the serial deterministic tile builder
for reuse; numerical threshold/intersection/order semantics are unchanged.
Fixed and isolated selected MPI runtime libraries include the new API.

Verified SHA-256:

- `src/amr/amr_patch_tree_tagging_3d_mod.F90`:
  `d842bb893559af2049d3ecf3666f3c0274beb6b0e1330ddf915ca63a323e4c3c`
- `src/parallel/mpi_amr_patch_tree_tagging_3d_mod.F90`:
  `d46e37dcdc54670f70a293fe72fa588f2820fe574894626b2088181b93b58bd0`
- `tests/unit/test_mpi_amr_patch_tree_tagging_3d.F90`:
  `a71816769afa0d435cad1da2fef03933532464208eea3fc5cb8d1af8e0f15055`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target test_mpi_amr_patch_tree_tagging_3d test_amr_patch_tree_tagging_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_(mpi_amr_patch_tree_tagging_3d_[23]|amr_patch_tree_tagging_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Build passed. CTest **3/3 passed**, 0.94 s total: MPI 2 ranks 0.45 s, MPI 3 ranks
0.45 s, serial tagging 0.01 s; no failures/skips. Provenance, project contract
and whitespace checks passed. No Python sources changed. Existing real-comparison
and integer-division warnings remain. No installed-prefix or broad matrix audit.

Each MPI fixture runs cyclic owners and all payload on the last rank. It compares
complete candidate metadata with serial: node count/order, parents, depths,
shapes, all patch bounds and ratios. Cases cover separated hot islands, nested
thresholds with ratios 2/3, physical versus periodic buffered corners, fine-only
hot cells across ratio changes, ignored covered hot parents, complete coarsening,
and zero requested levels. Existing serial remap/inventory checks also run.

Node budget, decreasing/nonfinite thresholds, lattice overflow, differing thresholds
on one rank and a NaN on one owner are rejected collectively. Failed calls preserve
destination allocations and complete metadata. No golden hashes or tolerances changed.

Logs: `build/amr-tree-mpi-tagging-checks/{build,ctest}.log`.
Binaries: `build/selected-tree-mpi-debug/tests/test_mpi_amr_patch_tree_tagging_3d`
and `test_amr_patch_tree_tagging_3d` in the same directory.

## Remaining scope

The application still rejects automatic adaptation controls; this API is the
geometry-generation layer for that integration. Hot boxes/temperatures are replicated
and tile intersection scans remain quadratic, so no production scaling is claimed.
MPI automatic-adaptation checkpoint policy and subcycling remain open. These are
implementation checks, not external solver parity or physical validation.
