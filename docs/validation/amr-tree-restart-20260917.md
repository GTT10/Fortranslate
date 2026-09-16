# Serial tree application checkpoint/restart — 2026-09-17

## Source and scope

Base: `b7d407e514367981e9f9e041841d4d07ca9f2025`, local `main`.
The implementation adds optional checkpoint controls to the fixed/selected tree
driver, canonical initial geometry/full schedule identity, continuation validation
and exact restoration of state, temperature, counters and budgets. The shared
lexical-path helper's mechanism-independent configuration module is now compiled
into the fixed core as well as the selected runtime.

Source SHA-256 at verification:

- `src/driver/amr_patch_tree_application_3d_mod.F90`:
  `d5ad8aab5ed9f8024f29c9180d02328534529e64eb1acf4967e34164cc672904`
- `tools/check_amr_tree_restart_3d.py`:
  `8654bbb941a037230c4a24a142a4c149d54cc99cfff790faf3e5492575cabef4`

## Completed checks

GNU Fortran 13.3.0, Debug, bounds checks/backtraces/FPE traps; MPI disabled.
Existing clean-selected-runtime build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^regression_amr_tree_restart_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

- Build succeeded for both frontends.
- Focused CTest: **1/1 passed**, 37.75 s; no skips.
- Provenance, project contract, whitespace and Python `ast.parse` syntax checks passed.
- Each frontend stops/restarts independently at steps 1, 2 and 3: moved ratio-3
  tree, root-only tree and final-time restored tree. Final fields and budgets are
  byte-identical to the uninterrupted baseline; final counters remain 3/3.
- Changed physics, changed event schedule, missing footer, invalid consumed-event
  count, inconsistent mass ledger, lexical output alias and unterminated namelist
  are rejected. Existing output bytes remain unchanged.

Final field SHA-256:
`6f24802d2e6effac50d84cbbe3e37707b1d7a47aed07d2c4c724317d77f20529`.
Budget SHA-256:
`6bc1abefa52fadc15a328bd6d155f4c847ab5f2bd6266fd40cad29cf4f82e1c2`.

Logs: `build/amr-tree-restart-checks/{build,ctest}.log`.
Per-case logs, checkpoints and CSVs:
`build/selected-amr-tree-debug/tests/amr_tree_restart_3d/{fixed,selected}/`.

## Development failures and limits

An initial build failed because the shared lexical-path helper was absent from
the fixed core source list; its module registration corrected the dependency.
The first restart test rejected a valid early snapshot: transverse momentum was
only `2.24e-18` in absolute inventory, with a `4.18e-27` ledger residual. The new
restart check now scales momentum components by the full momentum-vector
magnitude; its `2e-10` coefficient remains unchanged. The final test includes
intentional ledger corruption rejection.

This is focused serial implementation evidence, not an MPI, compiler/configuration
matrix, installed-binary audit, external PeleC parity or physical validation.
Path checks are lexical and do not resolve filesystem aliases. Physical write
failure can leave a partial checkpoint; cross-compiler context compatibility is
not established. Automatic tree tagging, subcycling/PLM donors and tree MPI remain
open.
