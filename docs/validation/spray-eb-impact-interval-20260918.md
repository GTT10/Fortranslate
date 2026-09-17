# EB impact remaining-flight interval — 2026-09-18

Implementation `e6c1075`, base `a388e22`.

The mixed-liquid EB impact module now advances a bounded population through
remaining flight, appends splash children and deposited films, removes consumed
parents, and commits all inventories and diagnostics together. Child IDs are
consecutive above the supplied watermark. Film thickness uses the existing
planar normal/proximity rule. Contact depth follows descendants, and capacities
include intermediate births and deposits. Properties and films remain frozen
apart from impact deposition; this does not yet connect gas or film evolution,
MPI allocation, or an application. ID-order processing is deterministic but
does not globally sort contacts by physical time.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled and
selected H2O2 tests-disabled build directories were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, total **0.04 s**, no test failures or skips.
The existing oblique synthetic mixture fixture now checks remaining rebound
and splash flight, deposited-parent removal, consecutive child IDs, component
conservation (absolute 1e-23 kg), and mass/momentum/energy with wall receipt
(relative 1e-12, component-wise scale floor 1e-20). Particle/film capacity failures
preserve the initial state. A fluid slab between two EB walls exercises two
rebounds in one interval and rolls back the first contact when limited to one.

Selected runtime library compilation, provenance, Python syntax and project
contract passed. Whitespace checking found one trailing test-source space; it
was removed and the check passed. No tolerances changed. No selected execution,
full matrix, external solver parity or physical validation was run.

Artifacts: `build/spray-eb-impact-interval-checks/{build.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
