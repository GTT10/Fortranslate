# Serial patch-tree pressure outlet application — 2026-09-18

Implementation `339b3e0`, base `136286f`.

Fixed and configure-time selected serial tree frontends enable pressure-outlet
input. The shared application forwards targets to the full split update with
common or recursive timesteps. Existing canonical checkpoint context already
serializes the entire reactive configuration, including pressure targets, so no
checkpoint format change is introduced here. Tree remap uses conserved volume
overlap and EOS recovery without physical ghost sampling.

MPI tree input/stepping remains unconnected. This evidence does not establish
physical validation, nonreflection, MPI parity or installation qualification.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6; reused Debug/Ninja fixed tests-enabled and
selected H2O2 tests-disabled configurations. No MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_amr_tree_reactive_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^regression_amr_tree_pressure_(common|subcycle)_restart_3d$' --output-on-failure
python3 -m py_compile tools/check_amr_tree_restart_3d.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The existing scheduled-regrid hotspot fixture adds x-high pressure 102325 Pa
with moving no-slip y walls, active chemistry and molecular transport. Both
fixed and generated selected elementary mechanisms run with common and recursive
timesteps. Each case stops after step 1 and resumes across the remaining scheduled
changes, reaching step 3/regrid 3. Final state and budget CSVs match uninterrupted
execution byte-for-byte within each mechanism/timestepping configuration.
Changed pressure (103325 Pa), removed pressure (zero), changed physics/schedule,
truncated checkpoint, inconsistent schedule position/ledger and malformed or
aliased checkpoint requests are rejected. Recursive cases also reject changed
timestepping mode. Existing tolerances remain unchanged.

Logs: `build/pressure-outlet-tree-app-checks/{provenance,build,build-selected,ctest}.log`.
Runtime fixtures and per-invocation logs:
`build/spray-amr-reactive-debug/tests/amr_tree_pressure_{common,subcycle}_restart_3d/{fixed,selected}/`.

Final CTest **2/2 passed**, zero skips, 126.91 s total (common 20.85 s,
recursive 106.02 s). Fixed/selected builds, provenance, Python syntax, project
contract and diff checks passed. No failed or interrupted test runs.
