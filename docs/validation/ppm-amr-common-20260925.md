# Common-step AMR tree PPM — 2026-09-25

Source: base `615b4e8`, implementation `090aec2`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-amr-common-checks/`.

`advance_amr_patch_tree_hydro_3d` now accepts `characteristic_ppm`.
It samples four ghost layers, retains the existing conservative limited
coarse-donor interpolation, applies physical/Mach/pressure boundary data,
and invokes the supplied-halo CTU kernel once per patch. The existing
integrated flux containers feed restriction/reflux and the composite boundary
ledger. The time-centered CTU result is committed after one stage; PCM/PLM
retain their two-stage SSPRK2 update. Failure preserves state, temperature,
boundary budget and maximum correction.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_hydro_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_ppm_3d$' --output-on-failure`:
  **1/1 passed**, 8.94 s test / 8.98 s total.
- The PPM variant reuses the existing single-root and five-refined-patch
  multilevel fixtures. Five boundary combinations per geometry cover periodic,
  outflow/slip wall, prescribed outlet pressures, extrapolated-pressure Mach
  inlet and prescribed-pressure Mach inlet. Single-root U agrees with regular
  boundary PPM within 2e-13 scaled error, T within 1e-9 K and boundary ledger
  within 1e-14. Composite volume change minus boundary budget is below 2e-12
  scaled error; positive T, species closure, and oversized-step rollback of
  fields and diagnostics pass. Existing tolerances were retained.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No failed tests occurred. Compiler warnings remain. The existing PCM/PLM test
variant, full configuration matrix, selected builds, recursive level stepping,
MPI, application input and physical/convergence validation were not run.
Those PPM integrations remain open. Conservative EOS transverse correction
still differs from the reference primitive-pressure equation.
