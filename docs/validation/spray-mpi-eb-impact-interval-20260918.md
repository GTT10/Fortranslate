# MPI EB remaining-flight transaction — 2026-09-18

Implementation `84a10ee`, base `382d668`.

The new collective interval repeatedly resolves EB contacts and continues
remaining flight. Children inherit local parent pressure, time and contact
depth; binary lookup uses the resolver's strictly increasing local IDs. Every
rank participates until no rank has remaining live motion. The interval privately
retains all intermediate parcels and patches and publishes compacted live
parcels, films, component inventories, watermark, wall receipt and hit count
only after success. Model rejection remains a collective terminal diagnostic.

This is force-free drift with frozen parcel pressure and wall properties.
Films are queried before each batch, including deposits from earlier batches;
events are not globally time-ordered. No gas feedback, film evolution or final
owner migration is performed here. Those remain required for source/application
integration.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled
configurations were reused. Test runtime checks and floating-point traps remain
enabled.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, total **0.99 s** (two ranks 0.46 s; four ranks
0.49 s). Selected library compilation, provenance, Python syntax, contract and
whitespace passed. No failures or skips. Selected execution was not run.

The existing remote-film/contact test now also advances one parcel held only
on the final rank across an oblique wall. Its four splash children continue
remaining flight and the consumed parent is removed. Global mass/momentum/energy
plus wall receipt retains relative 3e-12, scale floor 1e-20. A zero-contact limit
rejects after private impact construction and preserves the original position,
film count, ID watermark, wall ledger and hit-count sentinels.

A fluid slab between two EB planes exercises two collective rebounds in one
0.4 s synthetic force-free interval. A one-contact limit rolls back the first
rebound; a two-contact limit succeeds with a fluid final position. Empty ranks
participate in every batch. Existing remote-film, ID, geometry disagreement and
invalid-correlation assertions remain in the same test.

Artifacts: `build/spray-mpi-eb-impact-interval-checks/{build.log,build-final.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
No broad matrix, scaling, installed audit, external solver parity, convergence
or physical validation was run.
