# Automatic tree application evolution and continuation — 2026-09-17

## Source and behavior

Base: `6353fde810750bb355daa83f653bbb46ec127431`, local `main`.
The fixed/selected tree driver now accepts successful-step temperature adaptation
controls, runs the topology builder and conservative remapper, and checkpoints
after due updates. Identical proposals skip remapping. Restart restores generated
geometry and accumulated diagnostics without retagging the saved temperature;
context, update count and permitted geometry bounds are checked.

Verified SHA-256:

- Driver `src/driver/amr_patch_tree_application_3d_mod.F90`:
  `88391fbaf7986c97c0410cbf2cc51d06ba3fda30b487da7db7b8b0a8566f259b`
- Regression `tools/check_amr_tree_adaptation_3d.py`:
  `9fa5dcd7bb4900370775f4a64bc328a003577615faca1f71779e96eb123edc4d`

## Completed checks

GNU Fortran 13.3.0 Debug with bounds checks/backtraces/FPE traps; MPI disabled.
Existing build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^regression_amr_tree_(adaptation|restart)_3d$' --output-on-failure
ctest --test-dir build/selected-amr-tree-debug -R '^regression_amr_tree_adaptation_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

- Both frontends built successfully.
- Existing prescribed-schedule restart regression: **1/1 passed**, 37.70 s.
- Automatic adaptation regression, after fixture correction: **1/1 passed**,
  14.48 s. No skips.
- Provenance, project contract, Python `ast.parse` syntax and whitespace passed.

The automatic fixture uses a corner hotspot, hydro with chemistry/transport off,
root 4³, two refinement depths, ratios 2/2 and one-cell tiles. It advances seven
steps to 1e-9 s, generates ten nodes over three levels, and retains the full domain
volume. Both frontends stop after the first generated geometry and resume to
byte-identical fields/budgets versus continuous runs. Fixed and selected results
also match each other exactly. Changed thresholds/cadence and insufficient node
budget are rejected while existing output remains unchanged. Both optional
namelist group orders are exercised. The prescribed regression separately retains
its reacting/transport/wall fixture and three scheduled restart positions.

Final automatic field SHA-256:
`c53685c780c706f4a394118b4160597eaab499cac568f01ad959d7e72df4da82`.
Budget SHA-256:
`adb424363c98202b174ed27050bc4b4b2a2ca082467fbae15e2c830cf68c283e`.

Logs: `build/amr-tree-adaptation-checks/`. Per-case inputs, logs, checkpoint and
outputs: `build/selected-amr-tree-debug/tests/amr_tree_adaptation_3d/{fixed,selected}/`.

## Failure and scope limits

The first automatic fixture was rejected at input validation because transport
was disabled while inherited moving-wall velocities remained nonzero. Setting
those fixture velocities to zero corrected the input; no production validation
or tolerance was relaxed. The first combined CTest invocation therefore reported
one pass and one failure; the corrected automatic test was rerun separately.

This is serial implementation evidence, not MPI/scaling qualification, a full
compiler/configuration matrix, external PeleC parity or physical validation.
Automatic and prescribed schedules are separate modes. The application permits
up to 16 refinement depths; its node budget fails rather than truncating tags.
Restart checks structural compatibility and counters, not historical tag replay.
Coarse/fine donors remain PCM and all levels share a timestep. PLM donors,
subcycling and tree MPI remain open.
