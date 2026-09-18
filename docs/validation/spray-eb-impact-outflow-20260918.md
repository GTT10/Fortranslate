# EB impact/splash exterior outflow — 2026-09-18

Implementation `6ff500f`, base `228c00d`.

`advance_mixture_level_set_impacts` accepts optional `outflow_change(5+nc)`.
When present, each remaining-flight segment is clipped at its first root-domain
crossing. The original level-set impact helper handles an earlier EB hit;
its segment-relative fraction is converted back to the original remaining time
before parent and child continuation. If the exterior crossing comes first,
the parcel's mass, momentum, internal plus kinetic energy and component totals
are accumulated as negative system change, and its record is removed at commit.

The wall receipt remains separate. Deposited films stay in the domain. Shared
IDs are never recycled; exterior events consume the same lineage depth budget
as EB contacts. Parcels, components, films, IDs, hit count and both budgets
publish only after complete success. Without the optional output, the previous
root-contained behavior remains. Source/reactive/application wiring is not yet
connected for the combined impact/outflow mode.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed test was serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, test 0.01 s, total 0.04 s. Existing oblique deposition,
rebound, splash, repeated contacts and rollback checks remain. The new case uses
a 200 m/s mixed parcel striking the oblique EB and follows its splash for 0.1 s.
All four children escape: one film remains, watermark is 5, and hit count is 5
(one impact plus four outflows). Film inventory plus wall receipt minus exterior
change equals the initial mass/momentum/energy within relative 1e-12, floor
1e-20. Component closure uses absolute 1e-23 kg.

A depth-one attempt fails after private splash creation and preserves the
original particle/component population, no films, incoming watermark and both
ledger sentinels. Fixed/selected library builds, provenance, Python syntax,
project contract and whitespace checks passed. No failures or skips. Selected
runtime, broad matrix, installed audit, external solver parity and physical
impact/outflow validation were not run.

Artifacts: `build/spray-eb-impact-outflow-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
