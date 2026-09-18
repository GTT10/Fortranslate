# MPI EB impact and gas source coupling — 2026-09-18

Implementation `45b44d9`, base `8dd3240`.

The existing MPI AMR mixture source API now accepts film state and a shared ID
watermark in nodal EB mode. Each adaptive source substep composes gas exchange
with collective impacts and remaining child flight. Initial film validation and
final particle/film ownership migration are part of the same transaction.
The global wall receipt is subtracted once from the system boundary inventory.
An optional collective model-rejection diagnostic exits without timestep retries.

## Focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja; reused fixed tests-enabled and selected H2O2 tests-disabled builds.
The new test enables bounds/runtime checking and invalid/zero/overflow traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_step test_spray_mpi_eb_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(impact_step|step)_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_step test_spray_mpi_amr_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(eb_impact_step|amr_step)_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_step_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The initial impact/specular run passed **4/4**, 2.22 s. After adding the negative
model fixture, regular AMR passed **2/2**, while both new negative tests failed:
the initial film height 1e-4 selected a different impact condition instead of the
intended inadmissible splash. Changing that fixture height to 2e-6 (no changes
to production correlations or tolerances) made the final impact run pass **2/2**,
1.18 s (two ranks 0.55 s, four ranks 0.58 s). Thus six distinct tests have passing
results; the two failed intermediate executions remain recorded. No skips.
Selected library compilation, provenance, Python syntax, contract and whitespace
checks passed. No selected-mechanism execution was run.

The fixture has a four-cube root with left-half ratio-two refinement and
phi=y-0.00125 on the original eight-cube nodal lattice. One water/peroxide parcel
starts on the final rank. At 20 m/s normal velocity it deposits into one fine-owner
film; at 200 m/s it creates four outgoing children with IDs 2–5, retains their
lineage and advances their remaining flight. Both cases heat/evaporate before
impact and finish on the correct fluid-side owners. Composite gas, particle and
film mass, momentum, energy and species plus the wall budget close to relative
3e-12 with scale floor 1e-30.

A one-substep limit rejects the whole interval and preserves gas, temperature,
particle position/composition, empty films, IDs and ledger. A preexisting thick
film causes collective model rejection with unchanged gas, particles, film mass/
height, watermark and ledger. The legacy specular and regular AMR source tests
exercise the unaffected optional-argument paths.

Artifacts: `build/spray-mpi-eb-impact-source-checks/` contains `build.log`,
`build-final.log`, `build-model.log`, `build-selected.log`, `ctest.log`,
`ctest-final.log` (including the failed fixture), `ctest-model.log`,
`provenance.log` and `contract.log`.

Stationary-film thermal evolution remains a separate operator. Impact-enabled
reactive composition and MPI application/restart persistence remain open.
No broad configuration matrix, installed audit, external solver parity, scaling
or physical validation was run.
