# Owner-independent computational split IDs — 2026-09-18

Implementation `4426a36c37d391f9c1e563ebea6ef67c9684f544`, based on `a2f26e8`.

The shared MPI computational splitter previously reserved children in rank
order, allowing a changed owner layout to change child IDs. Children now order
by global parent ID and then deterministic sibling order. The existing global
reservation still validates the watermark and overflow before assignment.
Each rank broadcasts its compact parent/count metadata in turn; a merge over
sorted parent IDs computes offsets. Only one remote owner's metadata is held
at a time. No particle/state payload gather is introduced. Duplicate splitting
parent IDs across owners reject collectively before publication. Allocation or
assignment failure preserves caller parcels, component masses and watermark.

This applies to computational splitting across regular, AMR and EB callers of
the shared API. Injection, KH-release and splash ID reservation remain separate
rank-ordered paths. Saved IDs remain readable; future computational children
can have different IDs than in older builds. This is not cross-build trajectory
identity or a scaling qualification; metadata communication visits every rank.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `test_spray_mpi` and `pelef_mpi_amr_spray_3d`, and selected H2O2
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.
CTest **3/3 passed**, 21.53 seconds, zero failures/skips: AMR split application
(4.51 s), MPI migration at 2 ranks (8.20 s) and 4 ranks (8.78 s).
Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, full matrix, external solver parity or physical validation.

The MPI test uses the same four parent IDs with unequal child counts 1–4 under
round-robin, reversed and concentrated ownership. Both 2- and 4-rank executions
check the exact parent/sibling ID formula and final watermark 14, including
empty owners. Existing component inheritance, multiplicity/KH-inventory
conservation and capacity rollback checks remain. The prior rank-order assertion
was replaced by the new explicit parent-order contract; no tolerance changed.
The AMR split application continuation test checks gas/particle outputs through
2-to-4 rank restart using its existing tolerances.

## Commands and artifacts

Logs: `build/spray-split-id-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(migration_[24]|amr_split_application)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
