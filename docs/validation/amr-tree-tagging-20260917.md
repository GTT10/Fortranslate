# Temperature-driven tree topology API — 2026-09-17

## Source and scope

Parent commit: `5468589` on local `main`. This change adds a transactional
multi-level temperature tagging API to the fixed and selected tree runtimes.
It selects disjoint parent-cell tiles from exposed old-cell temperatures, with
per-depth thresholds/ratios, cross-parent/periodic buffers and a node budget.
It produces geometry for the existing conservative remapper. Application cadence
and automatic-geometry checkpoint continuation remain separate integration work.

Verified source SHA-256:

- `src/amr/amr_patch_tree_tagging_3d_mod.F90`:
  `8c905508b15f2568a265fc2b07c4b9aec1a18b8016f96c9e4e83524418f048ff`
- `tests/unit/test_amr_patch_tree_tagging_3d.F90`:
  `26968d13c1d4d8889512f67dc87b460b690b64c523dc7153ae1e9e2a72029f1a`

## Completed verification

GNU Fortran 13.3.0 Debug, bounds checks/backtraces/FPE traps, MPI disabled;
existing `build/selected-amr-tree-debug` build directory.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target test_amr_patch_tree_tagging_3d pelef_selected_amr_tree_reactive_3d_runtime -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^unit_amr_patch_tree_tagging_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Both runtime builds succeeded; focused CTest **1/1 passed**, 0.01 s, no skips or
failed test invocations. Provenance, project contract and whitespace checks passed.
No Python source changed, so Python syntax verification was not applicable.

The unit checks two disconnected root regions and a hotter nested region at
ratio 3; exact node ordering/bounds/shapes; conservative remapping onto generated
geometry; physical versus periodic corner buffering; fine-only hot-cell retention
when changing ratio 2 to 3; exclusion of covered hot parents; full coarsening;
root-only maximum depth; node-budget, decreasing-threshold, nonfinite-threshold
and integer-overflow rejection with destination preservation.

Artifacts: `build/amr-tree-tagging-checks/{build,ctest}.log` and
`build/selected-amr-tree-debug/tests/test_amr_patch_tree_tagging_3d`.

This is a focused serial geometry/API check, not automatic application evolution,
MPI, full configuration qualification, external solver parity, physical validation
or a scaling result. Fine-cell/tile intersection uses exact integer coordinates,
but hot-cell scans are not spatially indexed. Rectangular tiling does not claim
Berger-Rigoutsos clustering or proper nesting.
