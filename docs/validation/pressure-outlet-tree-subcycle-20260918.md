# Patch-tree subcycle pressure outlet — 2026-09-18

Implementation `c5eef74`, base `4585304`.

`advance_amr_patch_tree_subcycle_hydro_3d` forwards optional
`outlet_pressure(6)` to every recursive level update. Existing time-interpolated
coarse endpoints, flux synchronization, restriction and transactional state/
diagnostic publication remain unchanged. Transport/full/MPI tree wrappers and
application input are still pending.

## Completed focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, reused Debug/Ninja fixed tests-enabled and
selected H2O2 tests-disabled configurations; no MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_subcycle_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_subcycle_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, no skips, 7.65 s total. Fixed/selected builds,
provenance and static checks passed. Initial fixed build failed because the
extended signature exceeded the Fortran line limit; wrapping it resolved the
failure before tests were run. No tolerance changes.

The added PCM/characteristic-PLM cases set x-high/z-high targets to 102325 Pa.
Root-only results retain the existing common-step comparisons; the multilevel
case checks 1/2/6 executed steps, composite conservation, synchronized parent
states and species closure. A rejected large step preserves all fields,
boundary/reflux diagnostics and step counters. This is focused numerical
evidence, not physical validation.

Artifacts: `build/pressure-outlet-tree-subcycle-checks/` contains `provenance.log`,
`build.log`, `build-final.log`, `build-selected.log`, and `ctest.log`.
