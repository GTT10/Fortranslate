# Regular 3D tree checkpoint storage API — 2026-09-17

Source parent: `3b9462c793e291452546de83c6b81d2ba596378b` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
New module SHA256: `744947c10a6c9cb6c46087f3ee85cbee9d5518b50834ada0a8b10cc0cfb49d4e`.
Test SHA256: `52f42b499da7678d6df59f7f3a6d67a43b81e1c77c266ac6822ae905e2bb331c`.

## Implemented scope

The versioned text snapshot stores arbitrary valid tree topology, state,
temperature, time/step/regrid counters, initial and boundary inventories, and
maximum reflux/minimum limiter/maximum remap diagnostics. Writing validates
input before replacing a file; read failure preserves the entire destination.
Synchronization and the existing EOS-consistency check are applied to all nodes.
Physical write failures may leave partial files; there is no atomic-file claim.

Context generation captures all flow configuration except output path and
continuation time/step limits, ordered NASA7/transport records, and complete
reaction/composition/integrator/tolerance context using the existing reaction
serializer. Caller-supplied mechanism and geometry/schedule identities are also
captured. Context lines compare exactly; cross-compiler formatting compatibility
is not established. The old serializer's implementation is unchanged; its write
entry point is now public for this reuse.

These APIs are included in fixed and selected tree libraries. Application
checkpoint controls, stopping/resuming, continuation-limit checks and schedule
position validation are not yet connected. Historical budget correctness and
the meaning of the caller-supplied schedule identity are application concerns.

## Completed focused checks

Reused GNU Fortran 13.3 Debug builds: fixed unit target in
`build/spray-mpi-migration-gnu-debug` (OpenMPI 4.1.6, MPI/spray enabled), selected
library in `build/selected-amr-tree-debug` (MPI/spray disabled). Both use bounds
and invalid/zero/overflow checks. No compiler/backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_checkpoint_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree_checkpoint_3d$' --output-on-failure
cmake --build build/selected-amr-tree-debug --target pelef_selected_amr_tree_reactive_3d_runtime -j 4
python3 tools/check_project_contract.py --source .
git diff --check
```

Final builds passed and **1/1 focused CTest passed, 0 failures/skips, 0.08 s**.
The test round-trips a three-level tree with ratios 2/3, requiring exact fields,
temperatures, clocks, counters and budgets/diagnostics. It checks flow/transport
context mismatch rejection, permitted continuation-limit/output-path changes,
invalid-temperature write rejection preserving the old file, and a missing
footer after the complete payload preserving the prior destination.
No tolerance/golden change or failed test occurred. The final short reruns cover
added diagnostic persistence and late-read rejection checks.

Provenance, project contract and whitespace passed; both checker scripts passed
Python AST parsing (2/2). Logs: `build/amr-tree-checkpoint-checks/`.
Checkpoint artifacts:
`build/spray-mpi-migration-gnu-debug/tests/amr_tree_checkpoint_3d/`.
