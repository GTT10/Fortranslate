# Serial EB tree SDC checkpoint envelope — 2026-09-27

## Source and format

Implementation `a57e9c4a`, based on `489f1428` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The serial tree record now supports `PELEF_EB_TREE_CHECKPOINT_V3` for SDC.
It stores iteration count, chemistry-source presence, external force and the
cumulative composite forcing vector after the existing metadata. Optional
per-patch species/energy reaction estimates follow all gas records under
`SDC_REACTION_SOURCES`. The final envelope marker is unchanged. Embedded gas
records retain their legacy formats; the outer envelope binds SDC controls and
sources. Split readers/writers retain V1 and indexed V2 support.

`eb_tree_checkpoint_3d` adds forcing and packed reaction-source fields. The
composite inventory validation subtracts the independent forcing ledger.
Reader state and metadata publish only after every record, source value,
composite inventory and parent restriction check passes.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_eb_tree_checkpoint_3d` built successfully.
- CTest `unit_eb_tree_checkpoint_3d`: 1/1 passed, 0.51 s test / 0.56 s total.
- New V3 checks use a three-level/five-patch tree and synthetic finite source
  records. Forcing round trips exactly with chemistry enabled and disabled;
  signed reaction fields round trip exactly when present.
- Replacing the final reaction row with a NaN rejects the read after the gas
  records have parsed. Caller gas and source metadata remain unchanged.
- A mismatched external forcing configuration rejects the envelope.
- Existing legacy tree checkpoint/restart and corruption checks remain in the
  same executable.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_checkpoint_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_checkpoint_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-checkpoint-checks/`.
This is a serial record API check; SDC application routing and a coupled SDC
restart trajectory remain pending. MPI indexed SDC records, changed-rank SDC
restart, selected builds, full matrices, installed-binary audits and physical
validation were not performed. Caller-owned stream publication semantics are
unchanged.
