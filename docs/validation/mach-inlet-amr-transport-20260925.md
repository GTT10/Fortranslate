# Serial AMR Mach inlet molecular transport — 2026-09-25

Implementation `5e836e8`, based on `cf50a9a`.

Common-step and per-level molecular transport now forward optional Mach, inlet
temperature and pressure controls to the viscous tree ghost sampler at each
SSPRK2 stage. Recursive transport forwards the controls through every interval.
Ordinary reservoir faces still require their packed tables; Mach-only faces do
not. Existing conservative limiting, wall flux treatment, restriction, reflux
and transactional publication are retained. MPI transport, LES, coupled dispatch
and application timestep/input integration remain pending.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed unit targets and selected H2O2 spray-tree
runtime dependencies built. Selected execution was not run. **2/2 CTests passed**,
no failures or skips, total 71.70 s:

- `unit_amr_patch_tree_transport_3d`: 32.62 s. Added Mach 0.2/700 K with
  extrapolated or 120000 Pa inlet pressure and moving no-slip/isothermal walls.
  All four molecular terms are enabled. Root-only results agree with regular
  transport at existing tolerances; composite conservation meets the existing
  2e-12 scaled bound. Species closure, temperature positivity and large-step
  state/diagnostic rollback pass.
- `unit_amr_patch_tree_transport_subcycle_3d`: 39.03 s. Both Mach pressure modes
  pass root/common-step comparison, nested 1/4/36 step counts, synchronized parent
  cells, species closure, positive temperatures and rollback. Composite
  conservation meets the existing 2e-10 scaled bound.

Existing non-Mach cases remain included. No tolerance changed. Mechanism
provenance, Python syntax, project contract and whitespace passed. This is focused
numerical evidence, not nozzle validation.

## Commands and artifacts

Logs: `build/mach-inlet-amr-transport-checks/`, containing `build-fixed.log`,
`build-selected.log` and copied `ctest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_transport_3d test_amr_patch_tree_transport_subcycle_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_transport(_subcycle)?_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
