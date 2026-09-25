# Regular PPM inflow PLM flux replacement — 2026-09-25

Source: base `e066558`, implementation `8349866`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-inflow-plm-checks/`.

The supplied-halo CTU API accepts an optional six-face inflow mask. After all
CTU solves, pressure-traced `pelec` faces marked as inflow are recomputed from
fourth-order limited primitive PLM slopes and the pinned wave-tracing equations,
without transverse corrections. The existing cell-wide flattening coefficient
is reused. Primitive face states feed the PeleC Riemann entry directly. Final
fluxes are published only when every replacement succeeds.

The regular physical-boundary hydro dispatcher forwards `faces=='inflow'`.
This includes regular gas and serial-spray callers of that dispatcher. Other
solvers retain their existing behavior. AMR, distributed and patch-tree callers
have not yet been given local physical-boundary masks. This increment does not
claim complete boundary parity across those paths.

The new primitive predictor returns only variables consumed by `pc_cmpflx`:
density, three velocities, pressure and species fractions. It does not return
an independently traced internal energy because that field is not read by
this Riemann interface. Fourth-order slopes match the pinned default
`PeleC::plm_iorder=4`; alternative slope orders are not exposed here.

## Reference

PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, inspected via `git show`
in `/mnt/sn850/Projects/gtt2/external/PeleC`:
`Godunov.cpp::pc_low_order_boundary`, `PLM.H::plm_slope`, `pc_plm_d`,
`pc_plm_d_passive`, and `Params/param_includes/pelec_defaults.H`.
PLM.H SHA256 `cbdbee7d0848604b7fbba73ee9d6d05c82cebd530ae421cfaae2fe8388da04fb`;
its inspected contents are retained with artifacts. Actual `dt` is passed to
the predictor, following the reference call rather than its contradictory
comment about half-time prediction.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure`:
  **1/1 passed**, 5.75 s test / 5.79 s total; no failures/skips.
- Added linear PLM spatial edges, fully flattened prediction and CFL rollback;
  all six marked boundary planes change on the nonuniform fixture while every
  interior face flux remains exactly unchanged. Existing six-face Mach inlet
  inventory and failed-step rollback now run for both `hllc` and `pelec`.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No selected build, MPI/application run, full matrix, compiled PLM oracle,
convergence or physical validation was performed. Compiler warnings remain.
