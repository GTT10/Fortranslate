# Reactive EB impacts and evolving films — 2026-09-18

Implementation `a5eb0f36c89e9a09536c5685d05b21b7835f4f43`, base `39a39e5`.

The new allocated reactive API composes particle/film half-intervals around
the existing reactive EB gas update. Optional wall settings, film arrays and
an ID watermark are supplied together. All state and diagnostic publication
is transactional across both source halves and the gas step. The original
fixed-population specular interface now delegates to this shared core.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja. Existing fixed tests-enabled and
selected H2O2 tests-disabled configurations were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_reactive test_spray_eb_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(impact_reactive|reactive)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, total **1.66 s** (new test 1.48 s, original specular
test 0.15 s), no failures or skips. Selected runtime library compiled;
provenance, Python syntax, contract and whitespace checks passed. Selected
execution was not run. Existing application compiler warnings remain in the
build log; this change does not alter that source.

The new test uses synthetic water/peroxide parcels, a planar cut cell, implicit
H2O2 chemistry and all molecular transport switches enabled. Incident speeds
20 and 200 m/s exercise deposition and four-child splash with film evolution
enabled. First-half deposits enter the second source half; splash IDs are 2--5.
Gas, parcels and films close mass, momentum, energy and represented species
inventories against exterior/wall/chemistry ledgers within relative 3e-12,
scale floor 1e-30. The chemistry ledger is nonzero and covered cells are unchanged.

An invalid gas solver after the first source half preserves gas, parcels,
component arrays, films, watermark and diagnostics for both impact modes.
A second failure fixture starts 150 nm from the wall at 20 m/s with a 10 ns
step: the first 5 ns half precedes contact, while the second half reaches a
wall with zero film capacity. Its failure restores the already updated gas
and leaves all caller state/diagnostics intact. The source substep limit is
eight in this fixture to bound the rejection path.

Artifacts: `build/spray-eb-impact-reactive-checks/{build.log,build-final.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
This is serial library composition. MPI impact/film ownership, application
inputs/output/restart and broader physical qualification remain open. No full
matrix, external parity, installed audit or physical validation was run.
