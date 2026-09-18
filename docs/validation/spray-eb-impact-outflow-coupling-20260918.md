# Serial EB impact/outflow source and reactive coupling — 2026-09-18

Implementation `336988a`, base `e7ebca2`.

The allocated source interval now accepts impact/film arguments together with
all-six-face outflow. It forwards a private component-order exterior inventory
to the impact interval, then maps its component entries into vapor-species slots
(summing duplicate species). EB wall system change remains separate from exterior
system change. Each accepted substep accumulates both inventories; gas, particles,
films, watermark and diagnostics publish only after the whole interval succeeds.
Other exterior modes with impact films remain rejected.

The allocated reactive wrapper now accepts `particle_outflow` with impacts and
films. Both source halves feed their exterior inventories to the full-step
boundary output; the gas update and final splitting remain inside the existing
transaction. The serial application guard and MPI combined mode are not changed.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja. Reused fixed tests-enabled and
selected H2O2 tests-disabled builds. Tests executed serially; selected build is
compile evidence only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_reactive test_spray_eb_impact_coupling --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact_(coupling|reactive)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, 0 failures/skips, 2.44 s total. Source fixture adds
mixed-liquid splash followed by all four children escaping in later substeps,
one retained film and watermark 104. Gas/particle/film inventories close against
wall plus exterior change at existing relative tolerance 3e-12 with 1e-30 floor;
component exterior entries preserve the 0.25/0.75 mixture within 1e-23 kg. A
one-substep cap rejects the interval and preserves all caller state and ledgers.

Reactive fixture splashes next to a root face, retains a film and removes some
children, with the same composite conservation criterion. An invalid gas solver
after the first source half preserves gas, parcel/component/film arrays, watermark
and diagnostics. Existing deposition, film, KH, splitting and exterior tests are
retained. An earlier source-only run passed 1/1; it is not an additional distinct
test. No physical validation or external-solver parity was run.

Artifacts: `build/spray-eb-impact-outflow-source-checks/` contains
`provenance.log`, `build.log`, `ctest.log`, `build-selected.log`,
`build-reactive.log`, `ctest-reactive.log`, and `build-selected-final.log`.
The final source adjustment after checks only indented existing continuation lines.
