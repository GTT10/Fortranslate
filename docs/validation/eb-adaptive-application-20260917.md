# Temperature-adaptive EB application — 2026-09-17

Implementation: `4e7740bbe46b7d13126eaad0bbdc9aba0875e7dd`, based on
`edc1880`. Fixed and selected serial two-level applications now connect
exposed-fluid temperature tagging, conservative remapping and adaptive restart.

## Behavior and repair

`temperature_adaptive` selects this mode. Requested evaluation times use
`regrid_count`/`regrid_times`; threshold, coarse tile dimensions, buffer and
patch capacity use `tag_temperature`, `tag_tile_shape`, `tag_buffer` and
`tag_maximum_patches`. Prescribed regrid boxes are rejected in adaptive mode.
The tagger selects from exposed coarse/fine fluid cells at each accepted event
endpoint; no selected cells produces zero fine patches. Capacity overflow
rejects before output, without silently dropping selected regions.

The separate `PELEF_EB_TEMPERATURE_CHECKPOINT_V1`–V4 envelopes bind initial
geometry/topology, evaluation times and complete temperature policy. Completed
selected layouts are stored and checked for count, tile alignment, extent,
canonical order, overlap and coarse-domain bounds. Future layouts must be
empty until evaluated. Restored state/budgets and event index/time are checked
by the existing reader. Subsequent events use the resumed solution.
These checks do not reconstruct historical temperatures or establish
cryptographic checkpoint integrity.

Initial adaptive checkpoint tests exposed zero-filled solid cells in newly
created fine patches. The checkpoint's positive finite storage checks rejected
them. The remapper now copies parent conserved state/temperature into every
new child before prolonging fluid state; same-resolution overlap still copies
old payload exactly. This gives solid storage meaningful inherited values.
No checkpoint check or tolerance was loosened, and no active-cell numerical
update kernel or mechanism data changed.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`. Fixed elementary H2/O2 uses explicit chemistry;
selected full H2/O2 uses implicit chemistry and the existing selected H2/O2/N2
composition. Both cases use reaction/transport/hydro on an 8-cubed finest
lattice with a 4-cubed coarse grid and planar wall z=0.31.

At evaluation times 1.5e-8 and 2.5e-8 s, a 1000 K threshold selects eight
2-cubed coarse tiles. A 1e9 K threshold selects none and removes all fine
patches. Each scenario runs fixed common stepping and automatic subcycling.

| Runtime / mode | 1000 K threshold residual | No-hot residual |
| --- | --- | --- |
| Fixed / common | 1.3235535718548270e-16 | 4.2411110160574008e-16 |
| Fixed / subcycling | 2.3867651089517283e-16 | 7.7725141358634883e-16 |
| Selected / common | 1.2052353039586972e-15 | 4.3017347320122123e-16 |
| Selected / subcycling | 1.1045492192755039e-15 | 5.5789922041286307e-16 |

The checker independently sums physical-volume-weighted CSV fields and compares
all component totals with budget output. Restart before and after the first
evaluation reproduces final solution and budget CSV bytes exactly, including
the second event. Changing the buffer rejects without overwriting prior
outputs/checkpoint. Capacity one rejects the hot case without replacing output.
No-hot output contains 48 active level-0 cells.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_two_level_3d test_selected_full_reactive_eb_two_level_3d test_eb_remap_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(regression_(selected_)?adaptive_eb_two_level_3d|regression_scheduled_eb_two_level_3d|unit_eb_remap_3d)$' --output-on-failure --output-log build/eb-adaptive-app-checks/final-ctest.log
python3 -m py_compile tools/check_adaptive_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **4/4 passed**, total 68.81 s, zero failures/skips: remap unit 0.19 s,
prescribed schedule regression 8.35 s, fixed adaptive 19.23 s, selected adaptive
41.01 s. Initial fixed/selected adaptive runs failed at checkpoint serialization
of zero-filled solid storage; their log remains in
`build/eb-adaptive-app-checks/ctest.log`. After the parent-payload repair,
all final tests passed. Provenance, Python syntax, contract and whitespace
checks passed.

Build/final logs: `build/eb-adaptive-app-checks/`. Per-case artifacts:
`build/selected-tree-mpi-debug/tests/adaptive_eb_two_level_3d/` and
`selected_adaptive_eb_two_level_3d/`.

This is serial fixed-ratio, two-level temperature adaptation at prescribed
evaluation times. General multilevel/MPI EB, nonplanar adaptation qualification,
external PeleC parity and physical/convergence validation remain open. No full
configuration matrix or installed-binary audit was rerun.
