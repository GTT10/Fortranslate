# AMR Mach inlet ghost sampling — 2026-09-25

Implementation `b4d453a`, based on `3695c621204e9f31b608d14396fe356c7525b0b1`.

## Change and scope

The serial tree ghost wrapper and shared cell-source sampler accept optional
species, Mach, temperature and pressure controls. Active Mach faces reconstruct
the packed conserved state and temperature using the existing inlet kernel after
same-level or coarse donor selection, temporal blending and optional spatial PLM.
Reflection and moving-wall momentum/energy adjustments follow reconstruction.
Ordinary reservoir faces retain their existing behavior. Candidate publication
preserves the caller's output on failure. MPI request planning still bypasses
physical state reconstruction.

This is the sampler foundation. AMR gas/LES stepping, MPI wrapper forwarding and
application input are not connected by this change.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug; existing fixed and selected H2O2 build trees.
Fixed unit targets and selected spray-tree runtime built successfully. Selected
execution was not run. Existing compiler warnings remain in the selected build.

**2/2 CTests passed**, no failures or skips, total 1.44 s:

- `unit_amr_patch_tree_ghost_3d`: 1.38 s. Six faces, seven ghost layers,
  extrapolated/prescribed pressure and moving-wall corners match the regular
  Mach ghost builder exactly. An interior refined patch uses temporally blended
  coarse pressure/composition with the MC limiter. Invalid inlet temperature
  preserves the previous output. Existing peer, coarse and wide-wall checks pass.
- `unit_amr_patch_tree_plm_ghost_3d`: 0.01 s. Existing PLM checks pass.

Mechanism provenance, Python syntax, project contract and whitespace passed.
No nozzle or physical validation is claimed.

## Commands and artifacts

Build logs: `build/mach-inlet-amr-ghost-checks/build-fixed.log` and
`build/mach-inlet-amr-ghost-checks/build-selected.log`.
CTest log: `build/spray-amr-reactive-debug/Testing/Temporary/LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_ghost_3d test_amr_patch_tree_plm_ghost_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_(plm_)?ghost_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
