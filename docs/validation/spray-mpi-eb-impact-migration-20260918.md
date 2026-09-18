# Atomic MPI EB impact and ownership migration — 2026-09-18

Implementation `55a2e68`, base `e9bd65f`.

The new `advance_mpi_eb_tree_impacts` wrapper validates/migrates existing films
on private state, advances collective remaining flight, and migrates live
particles plus stationary films to the current finest owners. It retains input
parcel-pressure alignment until the final migration. All arrays, the global ID
watermark, wall receipt and hit count publish only after both migrations succeed.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled
configurations were reused. The existing impact test retains runtime checks and
invalid/zero/overflow floating-point traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **0.97 s** (two ranks 0.46 s, four ranks 0.47 s).
Selected library compilation, provenance, Python syntax, contract and whitespace
passed. No failures or skips. Selected execution was not run.

The added two-level fixture uses a 2-cube root, left-half ratio-two refinement
and the original 4-cube nodal field phi=x+y-0.5. One parcel on each rank crosses
the oblique surface, generating four children and a film. All final live parcels
and films belong to the finest patch on the final rank. Global mass, momentum
and energy plus wall receipt close within relative 3e-12, scale floor 1e-20.

With particle capacity five, local impact construction succeeds but final
particle reception rejects eight/sixteen arriving children. With particle
capacity 100 and film capacity one, particle migration succeeds privately but
final film reception rejects two/four films. Both cases retain original parcel
positions/masses, empty films, watermark and diagnostic sentinels. Raising both
capacities to 100 succeeds with the expected populations and IDs. Earlier
single-contact and repeated-contact checks remain in the same test.

Artifacts: `build/spray-mpi-eb-impact-migration-checks/{build.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
Gas feedback, evolving films during flight, reactive/source composition and MPI
application persistence remain separate integration work. No broad matrix,
scaling, installed audit, external parity or physical validation was run.
