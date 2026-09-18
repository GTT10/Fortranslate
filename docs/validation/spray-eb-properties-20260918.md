# Serial EB component property inputs — 2026-09-18

Implementation `385d8be`, base `d743296`.

The serial EB namelist now connects polynomial component density and heat
capacity, bounded Antoine pressure and existing Clausius vapor-surface caps to
the existing liquid initializer and validators. Unused entries require defaults.
The injection effective liquid is constructed from mixed density, heat capacity
and energy rather than inheriting component one's polynomial model. Existing
parcel/film source APIs receive the configured component table. Normalized
restart context includes all new inputs and rejects older contexts.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled builds. Runtime checks were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_properties_application|spray_eb_restart)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **2.44 s**: new properties application 1.06 s and
legacy constant-property restart 1.34 s. The property driver performs eight
launches across solid/Weibull/mass-weighted and hollow/lognormal/number-weighted
injection: full, stopped, resumed and invalid-density cases for each.

The synthetic H2O/H2 components use mass fractions 0.25/0.75 and densities
1030/865 kg/m3 at 350 K, polynomial heat capacities and water Antoine pressure.
Fresh final-time parcels verify diameters using the harmonic mixed density,
component masses, weighting and event dose. Continuous and resumed gas,
particle/film CSVs and printed cumulative budgets are byte-identical. A negative
component-density fit rejects the injected mixture before CSV output. This is
input/coupling evidence with illustrative coefficients, not physical validation
of an H2-containing liquid mixture.

Fixed/selected executable builds, provenance, Python syntax, contract and
whitespace checks passed. No failures or skips. Selected runtime, surface-cap
runtime, broad matrix, installed audit, external solver parity and physical
property qualification were not run.

Logs: `build/spray-eb-properties-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
Property application artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_properties_application/nozzle_et2mzlhv`.
