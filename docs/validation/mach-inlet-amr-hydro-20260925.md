# Serial AMR hydro Mach inlet integration — 2026-09-25

Implementation `6181988`, based on `b1d0a68`.

## Change

Common-step and per-level SSPRK2 hydro forward optional Mach/pressure controls
and inlet temperature to the shared tree ghost sampler at both stages. Reservoir
tables are required only for ordinary inlet faces. The recursive hydro driver
forwards the controls through every level interval, preserving coarse endpoint
interpolation, integrated fluxes, restriction and reflux. Invalid ghost generation
returns before candidate state or diagnostics are published.

MPI hydro, AMR transport/LES, full coupled dispatch and application controls are
not connected by this change. The existing interior-state CFL selectors are
unchanged; inlet-aware application timestep selection remains pending.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed unit targets and selected H2O2 spray-tree
runtime dependencies built successfully. Selected execution was not run.
**2/2 CTests passed**, no failures or skips, total 28.10 s:

- `unit_amr_patch_tree_hydro_3d`: 16.63 s. Added Mach 0.2/700 K inlet cases,
  PCM with extrapolated pressure and characteristic PLM with 120000 Pa inlet
  pressure. Root-only stepping matches regular hydro at existing tolerances;
  refined-tree composite conservation meets the existing 2e-12 scaled bound.
  Positive temperatures, species closure and nonphysical-step state/diagnostic
  rollback pass alongside the existing periodic, wall and outlet cases.
- `unit_amr_patch_tree_subcycle_3d`: 11.43 s. The same Mach cases exercise
  root/common-step comparison and nested ratios two and three, with step counts
  1/2/6. Composite conservation meets the existing 2e-10 scaled bound; parent
  synchronization, positivity, species closure and transactional rollback pass.

Mechanism provenance, Python syntax, project contract and whitespace passed.
No tolerances were changed. These are focused numerical checks, not nozzle
validation or full application qualification.

## Commands and artifacts

Logs: `build/mach-inlet-amr-hydro-checks/{build-fixed.log,build-selected.log,ctest.log}`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_hydro_3d test_amr_patch_tree_subcycle_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_(hydro|subcycle)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
