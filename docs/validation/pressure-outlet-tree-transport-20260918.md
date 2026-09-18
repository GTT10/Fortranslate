# Patch-tree pressure outlet transport — 2026-09-18

Implementation `d33489f`, base `59bd72b`.

Common-step and per-level tree transport APIs accept `outlet_pressure(6)`
and apply the shared tree pressure helper after viscous ghost sampling.
Only physical faces receive targets; peer/coarse-fine sampling retains its
existing behavior. Flux limiting, composite reflux and private state/diagnostic
publication remain unchanged. Level-specific forwarding is compiled; transport
subcycling/full wrappers, MPI tree paths and application inputs remain pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, reused Debug/Ninja fixed tests-enabled and
selected H2O2 tests-disabled configurations; no MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_transport_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Added root-only and six-node multilevel cases combine x-high/z-high targets of
102325 Pa with all-outflow or moving isothermal y walls, linear coarse donors
and viscosity/conduction/species diffusion/barodiffusion. Existing regular-grid
comparisons, composite conservation, species closure and large-step rollback
retain their tolerances. Pressure cases are separate from the unchanged
zero-target PCM-versus-linear donor comparison. Per-level pressure updates
are not independently exercised here. No physical-validation claim.

Artifacts: `build/pressure-outlet-tree-transport-checks/` contains `provenance.log`,
`build.log`, `build-selected.log`, and `ctest.log`.

Final CTest **1/1 passed**, zero skips, 26.81 s total. Fixed/selected builds,
provenance, Python syntax, project contract and diff checks passed. No failed
or interrupted runs and no tolerance changes.
