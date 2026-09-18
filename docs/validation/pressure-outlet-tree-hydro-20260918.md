# Patch-tree pressure outlet hydro — 2026-09-18

Implementation `26836e8`, base `2b57492`.

`apply_amr_tree_pressure_ghosts_3d` applies the existing pressure kernel to a
private packed conserved-state/temperature ghost buffer after ordinary tree
sampling. Integer tree origin/domain layouts determine which node faces touch
the physical root boundary; internal peer/coarse-fine faces remain unchanged.

Common-step and level-specific tree hydro APIs now accept optional pressure
and apply it after ghost sampling. Common-step RK2 retains its existing flux
accounting, composite reflux and transactional publication. Level-specific
forwarding is compiled but not independently exercised with pressure yet.
Subcycling wrappers, transport, full reactive wrappers, MPI tree callers and
application pressure input remain pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, existing Debug/Ninja fixed tests-enabled
and selected H2O2 tests-disabled configurations. No MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_hydro_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_hydro_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Added cases use 102325 Pa x-high/z-high targets with PCM/characteristic PLM,
a root-only hierarchy and a six-node multilevel tree. Existing regular-grid
state, temperature and surface-ledger comparisons retain their tolerances;
multilevel composite conservation and nonphysical-step rollback are checked
with pressure enabled. This is numerical evidence, not physical validation.

Artifacts: `build/pressure-outlet-tree-hydro-checks/` contains `provenance.log`,
`build.log`, `build-selected.log`, and `ctest.log`.

Final CTest **1/1 passed**, zero skips, 12.14 s total. Fixed/selected builds,
provenance, Python syntax, project contract and diff checks passed. No failed
or interrupted runs and no tolerance changes.
