# Recursive AMR PPM — 2026-09-25

Source: base `4bbf47b`, implementation `ca00866`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-amr-recursive-checks/`.

`advance_amr_patch_tree_hydro_level_3d` accepts `characteristic_ppm`.
It fills four ghost layers using limited coarse interpolation, samples coarse
time endpoints at the substep's starting fraction, and performs one CTU
prediction/update. The returned integrated flux is dt times the time-centered
flux. PCM/PLM retain two stages and half-dt flux weights. Existing recursive
scheduling accumulates child fluxes, refluxes and restricts after synchronization.
All candidate fields and integrated fluxes remain private until success.

The starting-fraction sampling follows the initial-state CTU predictor;
it does not use the second RK endpoint as another flux stage. This change
establishes dispatch and conservation, not coarse/fine convergence order
or exact reference pressure evolution.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_hydro_level_3d test_amr_patch_tree_subcycle_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_ppm_(level|subcycle)_3d$' --output-on-failure`:
  **2/2 passed**, level 1.32 s, recursive 9.74 s, total 11.11 s.
- Level fixture covers periodic and physical single-root comparison, peer
  flux matching, integrated flux divergence, inactive-level preservation,
  coarse-time donor influence, invalid time ordering and failed-step rollback.
- Recursive fixture covers the root and three refined patches with nested
  refinement ratios 2 and 3; update counts are exactly 1/2/6. Five boundary
  combinations cover periodic, outflow/slip wall, pressure outlets and Mach
  inlets with extrapolated/prescribed pressure. Composite conservation remains
  below the existing 2e-10 scaled tolerance; synchronized parent states,
  species closure, positive temperatures and fields/diagnostics/count rollback
  pass. No tolerances changed.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, and `git diff --check`: passed.

No failed tests occurred. Compiler warnings remain. The default PCM/PLM
variants, full matrix, selected/MPI builds, reaction/transport splitting,
AMR application selection, convergence and external/physical validation were
not run in this increment. MPI and application PPM integration remain open.
