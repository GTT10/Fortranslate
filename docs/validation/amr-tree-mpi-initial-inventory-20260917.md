# Distributed MPI tree initial inventory — 2026-09-17

## Source and change

Base `b7d56b87adaebf88168395309094efd21aa3f4e5`, local `main`.
The MPI application now computes the initial composite inventory using its
existing serial-order compensated distributed traversal. Before checkpoint read,
the same owner-local traversal with absolute cell values computes the reference
scale for cancellation-sensitive momentum ledger checks. The global initialized
state gather, its root copy and subsequent deallocation are removed.

Only sum/correction vectors cross ranks in tree/cell order. Per-node volumes,
covered-cell exclusions, arithmetic order, momentum scaling and restart tolerances
remain unchanged. Root checkpoint reading/validation and save gathers remain.
No mechanism input, numerical solver, golden hash or tolerance changed.

SHA-256 `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
`6d9c3f5207314cf79f762694062cccaba03cdd7b294ef6808706483dccb4402b`.

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_(restart|weighted_adaptation)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed/selected builds passed. CTest **2/2 passed**, 97.86 s total: static
restart 55.56 s, weighted adaptation 42.26 s; no failures/skips. Provenance,
project contract and whitespace checks passed. Existing compiler warnings remain.

Existing tests cover fixed/selected common-step reacting/transport restart and
recursive PLM hydro-only automatic adaptation with weighted owners. They compare
continuous/restarted CSV and budget bytes, serial/MPI checkpoint bytes, and
2-to-3-rank continuation. The static check includes one-rank continuation;
weighted adaptation also switches partition policy on restart. Initial/context,
boundary-ledger and corrupted-checkpoint rejection remain enforced.
No new tests or Python sources were added.

Logs: `build/amr-tree-mpi-initial-inventory-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_{restart,weighted_adaptation}_3d/`.

## Limits

Checkpoint save/load still holds root field payloads. The serial-order vector
traversal is deterministic but synchronizes once per node; no production scaling
benchmark was run. These focused checks are not a full matrix, installed audit,
external solver parity or physical validation.
