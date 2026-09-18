# MPI AMR coupled computational splitting — 2026-09-18

Implementation `3a14bfb1f75d9385c19576251b74b9b80f239589`, based on `22d2cd1`.

The regular MPI AMR coupled spray operator now honors the existing
`options%maximum_multiplicity`. Shared spray options are validated before
candidate work, a positive limit requires a global ID watermark, and the limit
agrees collectively across ranks. The existing MPI component-aware splitter
runs after gas, particle and optional film updates but before publication.
Failure retains all caller state, component masses, IDs and inout diagnostics;
`split_created` returns zero. Success reports the global number of IDs reserved
by this computational split, excluding earlier physical children. The app's
existing options input/checkpoint identity already carries the limit.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `test_spray_mpi_amr_reactive` and selected H2O2 tests-disabled
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.
CTest **1/1 passed**, `spray_mpi_amr_split_2`, 62.27 seconds (62.31 s total),
zero failures/skips. Provenance, Python syntax, project contract and whitespace
checks passed. No failed builds/tests, full matrix or external solver parity.

The nested three-node mixed-liquid fixture first exercises the existing
default-disabled common/recursive, reacting and physical-boundary paths. The
new periodic case uses ETAB with high surface tension and three parcels of
multiplicity 100. A limit of 50 creates three children, six total parcels and
ID watermark 6. Gas-plus-particle mass, momentum, energy and species inventories
are conserved at the existing `3e-12` relative tolerance. A capacity of three
rejects after candidate evolution while retaining original state/IDs; rank-wise
limit disagreement rejects collectively. No tolerance changes.
This is a focused coupled-library test, not application restart with splitting,
KH/film-plus-splitting qualification, scaling or physical breakup validation.

## Commands and artifacts

Logs: `build/spray-amr-split-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_split_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
