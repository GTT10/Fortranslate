# Embedded EB checkpoint records — 2026-09-17

Implementation `1069a2a5866912224768b3ee95db008e65bc278f`, based on
`caf4b10` on local `main`.

The existing EB checkpoint reader/writer accept optional `checkpoint_unit`.
When present, they borrow an already open formatted unit of the correct action,
ignore the path, and leave positioning and closure to the caller. A successful
embedded read stops after `END_CHECKPOINT`, without consuming a following record.
Caller-owned units stay open on validation or I/O failure as well. Failed reads
retain the existing per-record state/geometry/metadata rollback. The caller must
manage the stream position after a failed read.

All existing thermodynamic, reaction, transport, configuration, geometry, state
and selected-context validation remains in place. The record format and schema
are unchanged. Without the optional unit, path opening/closing and strict
end-of-file validation retain their original behavior. This provides reusable
validated records for a multi-grid container; it does not yet implement the
container transaction or connect two-level application restart. A composite
reader must stage each record privately and validate its own final marker/EOF
before publishing any hierarchy state.

## Independently completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching
`build/selected-tree-mpi-debug`, MPI on, spray off. Fixed core and selected EB
runtime compiled. No MPI or application matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_reactive_eb_3d_checkpoint pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_reactive_eb_3d_checkpoint$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 PASS, 0 failures, 0 skips**, test 0.10 s, total 0.13 s.
The expanded existing checkpoint unit retains its standalone/fixed/selected
roundtrip, compatibility and corruption tests. New checks place two validated
records between caller-owned prefix/suffix records in one formatted scratch
unit; both states and distinct times/steps roundtrip exactly. Unit ownership is
preserved after successful reads/writes and failed reads. A non-checkpoint suffix
rejects without mutating targets and can be recovered by caller positioning.
Closed and unformatted units are rejected without closing an open caller unit.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-checkpoint-record-checks/{build,ctest}.log`.
This verifies the record interface, not two-level checkpoint/restart completion.
