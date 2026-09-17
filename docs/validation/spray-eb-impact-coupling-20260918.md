# Serial EB impact and gas exchange — 2026-09-18

Implementation `9a0c718`, base `543b8e0`.

An allocated-population source core now shares the existing EB gas-exchange
algorithm with the fixed-population specular wrapper. Optional wall settings,
film arrays and an ID watermark activate the impact interval. Each parcel uses
its sampled gas pressure for boiling/impact correlations; its children inherit
that frozen pressure for remaining flight. Source deposition retains pre-step
positions. Population-size changes resize sampling buffers on the next substep.
The entire interval, including film inventories and watermark, is transactional.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja. Existing fixed tests-enabled and
selected H2O2 tests-disabled configurations were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_step test_spray_eb_impact test_spray_eb_impact_coupling --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(impact_coupling|step|impact)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **3/3 passed**, total **0.11 s**: original specular source regression,
impact/remaining-flight regression, and new coupled impact test. No failures or
skips. Selected runtime library compiled; selected execution was not run.
Provenance, Python syntax, contract and whitespace checks passed.

The new synthetic water/peroxide fixture uses a planar cut cell, hot nitrogen,
drag/heating/evaporation and a 300 K wall. Deposition removes the parcel and
retains one heated film with less mass than the incident population. A 200 m/s
case produces four children with IDs 2--5, which enter subsequent source
substeps. Gas, droplets and stationary films close total mass, momentum, energy
and represented species inventories against the wall ledger within relative
3e-12, scale floor 1e-30. A one-substep limit verifies whole-interval gas,
particle, component, film, ledger and watermark rollback. Covered cells remain
unchanged in the deposition fixture.

Artifacts: `build/spray-eb-impact-coupling-checks/{build.log,build-final.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
This is a serial source interval. Films do not yet exchange heat/mass with gas;
reaction/hydro composition, application wiring and MPI impacts remain open.
No full matrix, installed audit, external parity or physical validation was run.
