# Pressure-traced CTU application dispatch — 2026-09-25

Source: base `5d37d02`, implementation `97097af`. GNU Fortran 13.3.0,
OpenMPI 4.1.6, Ninja Debug, Python 3.12.3. Reused fixed tests-enabled
`build/spray-amr-reactive-debug` and selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug`. Artifacts:
`build/ppm-pressure-dispatch-checks/`.

The common periodic and supplied-halo CTU entries now default to separate
primitive pressure work and multidirectional flattening when solver is `pelec`.
Existing application input `reconstruction="characteristic_ppm"` with
`riemann_solver="pelec"` therefore enables the policy without a new input field.
Other solvers retain conservative/EOS correction. Explicit kernel overrides
remain available for comparisons. This dispatch reaches regular, AMR, MPI and
spray callers of these CTU entries; it does not change PCM/PLM or EB hydro.

Restart context already binds reconstruction and solver. The application check
now additionally verifies that changing `pelec` to `hllc` is rejected. The test
establishes continuation with this implementation, not bitwise continuation of
old executables whose `pelec` PPM policy differed.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d test_reactive_ppm_halo_3d pelef_amr_reactive_3d pelef_mpi_amr_reactive_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_reactive_3d_runtime --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_step_3d|unit_reactive_ppm_halo_3d|mpi_amr_two_level_ppm_restart_3d)$' --output-on-failure`:
  **3/3 passed**, 50.22 s total; no failures/skips.
  - Step: 3.06 s; four-step conserved inventory/species closure, rollback,
    zero-step identity and exact periodic/physical transport-split dispatch
    for both `hllc` and `pelec`.
  - Halo: 5.09 s; default periodic/halo policy exactly matches explicit pressure
    tracing, with existing pressure-work, shock, subdomain and boundary checks.
  - MPI application: 42.03 s; reacting/transporting pressure-outlet and moving
    isothermal wall, dynamic regrid, exact coarse/fine CSV agreement with serial,
    four-to-two-rank restart, boundary ledger, and mismatched restart contexts.
- Application artifacts:
  `build/spray-amr-reactive-debug/tests/mpi_amr_two_level_ppm_restart_3d/run-vob31m5h/`.
- `python3 -m py_compile tools/check_project_contract.py tools/check_amr_physical_boundary_restart.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

The separate serial application CTest was updated to select `pelec` but not
rerun; serial reference execution is included in the MPI application check.
Selected runtime was built, not executed. Spray and patch-tree runtime tests,
full matrix, compiled external parity, convergence and physical validation were
not run. Compiler warnings remain. The project is not complete PeleC parity.
