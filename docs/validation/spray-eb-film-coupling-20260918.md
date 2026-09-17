# Stationary EB film heat and evaporation — 2026-09-18

Implementation `6baf70a`, base `862253e`.

The new serial film source interval reuses the established mixed-film rates and
exchange kernel, including trace depletion. It validates film locations/normals
against the original nodal EB surface, samples just inside fluid, and deposits
species/energy through the EB volume-weighted stencil. Film substeps limit
component loss and temperature change, retry rejected updates and publish the
whole interval transactionally. The wall ledger is system energy gain.

The allocated particle coupling optionally advances existing films after
particle gas deposition and before impact. Newly created films begin exchange
on the next particle substep. This is a first-order source ordering, not exact
contact-time film evolution. Conductivity is a shared constant; lateral film
transport, temperature/port fits, reactive composition, MPI and application
wiring are not included.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja. Existing fixed tests-enabled and
selected H2O2 tests-disabled configurations were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_coupling test_spray_eb_step test_spray_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(impact_coupling|step|impact)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **3/3 passed**, total **0.13 s**, no failures or skips. Selected
runtime library compilation, provenance, Python syntax, contract and whitespace
checks passed. Selected execution was not run.

The coupled synthetic water/peroxide fixture first deposits a film, then uses
a 400 K wall with no surviving droplets. Its film heats and loses mass while
the wall supplies positive energy. Gas plus film inventories close against
the wall ledger within relative 3e-12, scale floor 1e-30. Reversing the film
normal rejects the interval and preserves gas, film mass and ledger. The
200 m/s splash fixture now also enables film evolution while children continue
gas exchange; its gas/particle/film/wall inventory check passes with the same
tolerance. The original specular and impact-interval regressions also passed.

Artifacts: `build/spray-eb-film-coupling-checks/{build.log,build-final.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
No full matrix, external solver comparison, installed audit, convergence or
physical validation was run.
