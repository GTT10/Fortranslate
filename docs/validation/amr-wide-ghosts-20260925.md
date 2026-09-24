# Configurable AMR donor halo width — 2026-09-25

Implementation `3e86576ae4ecc89f750d5dfbe6e05805e3e56813`, based on `c70743d`.

Serial tree ghost filling and its shared cell-reader sampler now accept optional
`ghost_width`, default two, range 1--128. Shapes and global coordinates use that
width. Periodic modulo mapping already supports multiple crossings; widths
larger than the level domain are rejected on reflecting-wall axes to preserve
the single-reflection contract. Invalid widths preserve the caller's output.
This supplies the wider donor stencil needed by Dynamic LES; no LES time step,
MPI wrapper extension or variable-width pressure correction is claimed here.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, existing fixed tests-enabled and selected H2O2
tests-disabled configurations. Fixed `test_amr_patch_tree_ghost_3d` and selected
`pelef_selected_amr_tree_reactive_3d_runtime` built successfully. Selected
execution was not run.

One direct `test_amr_patch_tree_ghost_3d` invocation passed. New seven-layer
periodic halos on the existing multilevel/multipatch fixture agree exactly with
an independently prolonged dense grid over every interior and halo cell. This
includes cross-parent peer donors, ancestor fallback and multiple periodic wraps.
Zero width rejects without changing output. Existing two-layer physical and
periodic boundary, wall-motion, temporal interpolation and failure checks also
pass. No tolerance changed.

Provenance, Python syntax, project contract and whitespace checks passed.
The earlier MPI CTest/launcher handles 2018710/2018713 remained present on this
turn's initial observation; no new MPI test was launched. No MPI success, full
matrix, external parity or physical validation is claimed.

## Commands and artifacts

Logs: `build/amr-wide-ghost-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `test.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_ghost_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_amr_tree_reactive_3d_runtime --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH build/spray-amr-reactive-debug/tests/test_amr_patch_tree_ghost_3d
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
