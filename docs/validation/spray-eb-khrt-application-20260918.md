# Serial EB KH reactive/application connection — 2026-09-18

Implementation `3967afb`, base `3e1d9a9`.

The allocated reactive wrapper requires shared IDs for KH children and passes a
private watermark through both source halves, including the specular path.
Gas-step failure after private release preserves the incoming population,
components, gas, temperature, IDs and diagnostics. The application accepts KH
mode and validates restart IDs as injection plus computational split counts
plus four-child KH/splash groups. Existing checkpoint fields retain pending KH
state and lineage; no additional format is introduced for KH mode alone.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed tests were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_reactive pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_impact_reactive|spray_eb_khrt_application|spray_eb_split_application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_khrt_application$' --output-on-failure
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Reactive and existing split-application tests passed in the initial run (1.91 s
and 1.55 s). The new KH application fixture initially failed its child-population
assertion: at gas density 0.31 kg/m3, 1000 m/s and 20 ns, representative diameter
changed without a four-child release. Raising dose from 1e-16 to 1e-14 kg alone
also did not trigger release. The final fixture uses gas density 1 kg/m3 and
pressure 435483.87096774194 Pa, preserving the initial pressure/density ratio.
No model equation or tolerance changed. Final KH application CTest **1/1 passed**,
test 1.30 s, total 1.34 s, including the added altered-ID rejection. Thus all
three targeted tests have passing final evidence; there were no skips.

The reactive fixture checks gas failure after first-half KH release and success
with conserved composite inventory to relative 3e-12, scale floor 1e-30. It uses
a retired watermark of 100. The application starts two 10-micrometre water
parcels at 1000 m/s with total dose 1e-14 kg; drag, heating and evaporation are
disabled while gas chemistry, hydro and transport remain enabled. At 10 ns,
the checkpoint contains 18 records and watermark 18; at 20 ns it contains 34.
Continuous and resumed gas/particle/film CSVs and printed cumulative budgets
are byte-identical. Increasing the saved watermark by one rejects restore with
`Invalid restart child IDs`.

The final application driver passed nine launches, including expected failures
for immutable publication, changed liquid/geometry context, truncated archive,
altered IDs and input/output aliasing. Fixed/selected executable builds,
provenance, Python syntax, project contract and whitespace checks passed.
Selected runtime, broad matrix, installed audit, external solver parity and
physical breakup validation were not run.

Logs: `build/spray-eb-khrt-app-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, initial `ctest.log`, dose-only `ctest-app-final.log`,
`ctest-app-density.log`, final `ctest-app-ids.log`). Final application artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_khrt_application/eb_spray_7mgh6tui`.
