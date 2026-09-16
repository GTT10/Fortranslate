# Selected-mechanism regular 3D tree frontend — 2026-09-17

Source parent: `56c86012a4931e2dd93449af5e883da72ae2b45f` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Pinned mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Generated elementary fixture bundle SHA256:
`9b5c133fed08c62c810ccc22ff5834e7a650f4ecab670934e8660bd5f51dd15f`.
Template SHA256: `022b92865ff0b9448557aae69840ba5aec6bc84eb5753c4f9578d48b3a4ae9b5`.
Driver SHA256: `c35961e7f2833efe472e9d8c3bb5bd94971759338699aa934f8d0066d2560d60`.
Regression SHA256: `1670906d7cbce17cfd37627687c6c3dc40b0a1e7431ccf1610426e57badba540`.

## Scope

The bundle factory accepts `AMR_TREE_REACTIVE_3D_NAME` and configures a frontend
from the checked-in template. A dedicated selected tree runtime compiles the
existing topology/ghost/reflux/hydro/source/transport/remap/driver sources against
selected lower-level libraries. It does not link the fixed core. The frontend
loads generated data, checks composition and requested temperature range,
prints bundle identity and forwards the generated reaction integrator through
the shared driver. Fixed callers retain their existing default integrator.

Configure-time production bundles register `pelef_amr_tree_reactive_3d_selected`
and include it in their install target list. This turn verifies the elementary
fixture target; it does not qualify every bundle, install artifacts, MPI tree
execution or checkpoint/restart.

## Completed focused checks

New clean `build/selected-amr-tree-debug`: GNU Fortran 13.3 Debug with bounds and
invalid/zero/overflow checks; MPI off, default spray off, tests enabled.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake -S . -B build/selected-amr-tree-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DPELEF_ENABLE_MPI=OFF
cmake --build build/selected-amr-tree-debug \
  --target test_selected_elementary_amr_tree_3d pelef_amr_tree_reactive_3d -j 4
ctest --test-dir build/selected-amr-tree-debug -R '^regression_selected_amr_tree_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Build and **1/1 focused regression passed, 0 failures/skips, 9.29 s**.
The test compares fixed and generated elementary mechanisms through reacting,
transporting, moving-wall evolution with three scheduled topology changes.
Both exposed-cell CSV and budget CSV match byte-for-byte:

- Field SHA256: `6f24802d2e6effac50d84cbbe3e37707b1d7a47aed07d2c4c724317d77f20529`.
- Budget SHA256: `6bc1abefa52fadc15a328bd6d155f4c847ab5f2bd6266fd40cad29cf4f82e1c2`.

Unknown species and out-of-mechanism-range initial temperature are rejected
without changing existing field output. The selected log identifies the bundle
and `explicit` integrator. `nm` finds no fixed H2/O2 database references in the
tree runtime archive, and the executable link command contains no `pelef_core`.
(The generated elementary fixture itself intentionally uses its generated
mechanism module name.) No test failure or tolerance/golden change occurred.

Provenance, project contract, whitespace and Python AST parsing of the two
checkers plus new regression passed (3/3). Evidence logs:
`build/selected-amr-tree-debug/evidence/`; case artifacts:
`build/selected-amr-tree-debug/tests/selected_amr_tree_3d/`.
