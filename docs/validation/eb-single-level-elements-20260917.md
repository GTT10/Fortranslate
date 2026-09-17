# Single-level EB source element inventories — 2026-09-17

Implementation: `68d943a504efb4e9247c9daf40d03dc7423f1517`, based on
`8fdab92` on local `main`.

The selected single-level frontend now passes generated source element names
and atom counts to the shared planar/extruded/level-set application. Generic
inventories replace the species-name whitelist and bind schema-3 checkpoints.
Legacy library callers without metadata retain the original H/O/N path and
schema. The selected initial inventory is populated even when chemistry is off.
Existing relative conservation thresholds are unchanged; the generic path
subtracts a separately bounded machine-roundoff contribution before assessing
relative elemental error, including initially absent elements. The bound uses
64 epsilon times accepted steps, initial species mass and maximum atoms per
unit mass, as in the two-level frontend.

The former unsupported-label rejection fixture now runs A/B chemistry using
its source X-element table. This is synthetic metadata coverage. The selected
restart checker verifies the schema-3 H/N/O table against independent atom
counts, while preserving its existing fixed/selected state comparisons.
The long-interface compiler test's manual dependency list was also missing
the existing polyhedron clipping and EB reflux modules; both were added.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing `build/selected-tree-mpi-debug`.
The optional methanol fixture uses the previously configured
`PELEF_TEST_EB_METHANOL_BUNDLE=build/port-methanol-reference/ffcm1_red.json`
(absolute path in the cache). All checks here are serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target pelef_reactive_eb_3d test_selected_full_reactive_eb_3d \
  test_selected_unsupported_element_reactive_eb_3d test_selected_methanol_eb_3d \
  test_selected_elementary_reactive_eb_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(regression_selected_reactive_eb_3d_full_run|regression_selected_reactive_eb_3d_generic_element_run|regression_selected_methanol_eb_3d|regression_reactive_eb_3d_state_run)$' --output-on-failure
# Rerun the two corrected positive inputs:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(regression_selected_methanol_eb_3d|regression_selected_reactive_eb_3d_generic_element_run)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(regression_reactive_eb_restart_3d_reference|regression_selected_reactive_eb_3d_restart_(reference|checkpoint_stop|restart|composition_mismatch_rejected|check))$' --output-on-failure
python3 tests/python/test_generate_mechanism_bundle.py \
  GenerateMechanismBundleTest.test_long_interfaces_compile_with_probe_and_reactor
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json
python3 -m py_compile tools/check_selected_carbon_eb_3d.py tools/check_selected_reactive_eb_3d_restart.py tests/python/test_generate_mechanism_bundle.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Completed passes: fixed state run 0.11 s; H2/O2 selected run 13.57 s;
methanol/generic-label rerun **2/2 passed** in 2.09 s (2.03/0.02 s);
restart group **6/6 passed** in 9.97 s; long-interface Python/compiler check
**1/1 passed** in 23.808 s. Final passing groups had zero failures/skips.
Provenance, Python syntax, project contract and whitespace passed.

The uniform-density 2-cubed methanol fixture uses implicit chemistry and
transport, completing five steps to 1e-11 s. All eight CSV cells have finite,
positive density/temperature and species closure. Reported invariant and
roundoff-adjusted element errors are zero; minimum transport limiter is 1.
Checkpoint continuation produces byte-identical CSV output. Changing the
checkpoint's carbon label rejects without overwriting that output.

Failures retained as limitations: the initial density-contrast methanol input
changed total mass from 2.4800000000000000e-13 to 2.4799999998720304e-13 kg,
relative difference 5.1600664166805123e-11, exceeding the existing 5e-11 gate.
That input is retained in `build/eb-single-element-checks/density-contrast-methanol.nml`.
It was replaced in the connection fixture by uniform regular/cut density;
no invariant tolerance was relaxed. The application still lacks a general
external-boundary exchange ledger, and the density-contrast result is not
qualified. The former startup-rejection fixture also reached its default step
limit, so its new positive input explicitly uses 1e-9 s. A template line-length
failure and missing compiler-test dependencies were corrected before the
compiler test passed; a diagnostic line was split after a Fortran line-length
build failure.

Methanol source SHA256:
`ddbba8c97d929530cfb213da45485cca4f979dccf97c1063721bd855d7b3f259`.
H2/O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-single-element-checks/`; carbon artifacts:
`build/selected-tree-mpi-debug/tests/selected_methanol_eb_3d/`.
No new extruded/nodal carbon case, broad matrix, installed audit, external
solver parity or physical validation was run.
