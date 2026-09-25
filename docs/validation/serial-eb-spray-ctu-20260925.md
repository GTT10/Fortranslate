# Serial EB spray CTU — 2026-09-25

Implementation `1114ef5`, parent `ac4c6ce`. PeleC pin remains
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Both fixed-array and allocated serial spray APIs now accept an optional extended
geometry and forward it into the coupled EB gas step. Fixed/selected applications
accept `hydro_method='plm_ctu'` with the `pelec` solver, build six-layer configured
geometry once, and forward it in both impact/no-impact dispatch branches.
Checkpoint context adds `hydro-v1 plm_ctu`; PCM context stays unchanged.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; existing builds
`build/spray-amr-reactive-debug` and `build/spray-amr-reactive-selected-debug`,
the latter selecting `mechanisms/h2o2_full.json`. No toolchain changes.

- `spray_eb_reactive`: PASS, 0.29 s.
- `spray_eb_reactive_ctu`: PASS, 1.00 s.
- CTest total: 2/2 passed, 1.33 s, no skips/interrupted runs.
- Fixed and selected application checker invocations with `--ctu --restart`:
  both PASS, including exact gas/particle/film CSV continuation and retained
  existing conservation, reflection/transfer and output-protection checks.
  Wrong hydro-method continuation rejects and preserves all existing outputs.

The unit fixture uses a faceted plane, two liquid components, reactive/transport
coupling and StateRedist. It checks gas/particle/wall inventory, heating,
evaporation, reflection, external acceleration, covered storage and transactional
failure. Existing tolerances were not changed. Application liquid films are empty
in these reflection cases; this does not qualify CTU with active film impacts.
No new external-solver comparison or physical nozzle/LES qualification is claimed.

Initial compilation found one overlong new call line; corrected before testing.
The first selected target name was incorrect; the actual configured target is
`pelef_reactive_eb_3d_selected_spray`. Its build and runtime check passed. Compiler
warnings remain. The new `spray_eb_ctu_restart` CTest registration was checked
with `ctest -N`; equivalent checker commands were run directly.

Commands (repository root; builds used LD_LIBRARY_PATH/LIBRARY_PATH unset):

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-amr-reactive-debug --target test_spray_eb_reactive pelef_eb_spray_3d --parallel 2
cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_reactive(_ctu)?$' --output-on-failure
python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/serial-eb-spray-ctu-checks/fixed --ctu --restart
python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/serial-eb-spray-ctu-checks/selected --selected --ctu --restart
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts under `build/serial-eb-spray-ctu-checks/`: build logs, `LastTest.log`,
`fixed/eb_spray_syr_n1xq/`, `selected/eb_spray_b_330wxg/`. Application directories
retain run/restart/rejection logs, inputs, snapshots and final CSVs.
