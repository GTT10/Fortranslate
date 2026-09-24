# Primitive-state PeleC Riemann entry — 2026-09-25

Source: base `364dadf`, implementation `f5b91cd`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3. Reused fixed tests-enabled
`build/spray-amr-reactive-debug` and selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug`. Artifacts:
`build/ppm-primitive-riemann-checks/`.

`reactive_pelec_primitive_flux_x` accepts density, three velocities, pressure,
mass fractions and the adjacent cell sound-speed mean directly. It validates
thermodynamics from pressure/density without recovering pressure from a
conserved normal predictor. A shared private solver retains the existing
conserved-entry policy. The primitive entry uses the pinned reference's
stationary threshold, wave-speed regularization and star-pressure floor, and
keeps the pressure difference grouped before addition to the velocity numerator.
It publishes fluxes only on success and exposes the interface density, normal
velocity and pressure needed for subsequent transverse pressure work.

This is the Riemann prerequisite for reference pressure tracing. The existing
CTU kernels still use their conservative transverse correction; the new primitive
entry has not yet been connected there. The source-aligned arithmetic does not
establish full external-solver parity. PeleF thermodynamic admissibility limits
and species-flux closure remain in force.

## Reference evidence

Inspected pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` with
`git -C /mnt/sn850/Projects/gtt2/external/PeleC show <commit>:Source/<file>`:

- `Godunov.H::pc_cmpflx` forwards reconstructed pressure directly to `riemann`;
  `pc_transdo`/`pc_transdd` evolve primitive pressure separately from total energy.
- `Riemann.H`: SHA256 `073ae14b2f7925de2858ef312d02fc50f033c2a344dccc5b24ecde1f77ed8253`.
- `Constants.H`: SHA256 `0b2757e83a7e608b3782049ab5215ae3e439230b9fc0c0624cdc096207f57d48`.

The last two headers are retained in the artifact directory. Constants used by
the primitive entry are smallu `1e-12`, small_num `1e-8`, and very_small_num
`epsilon*1e-100`. The conserved entry retains its earlier roundoff policy.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_pelec_riemann --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_reactive_1d_runtime --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_pelec_riemann$' --output-on-failure`:
  **1/1 passed**, 0.00 s reported test / 0.05 s total; no final failures/skips.
- Existing conserved contact/shock/directional tests remain. Added primitive
  equal-state physical flux, `1e-8` m/s contact-speed retention, exactly retained
  prescribed pressure, unequal-state general-EOS agreement, and invalid pressure/
  sound-speed rejection with zero outputs. No tolerance was loosened.
- Initial compilation rejected an accidentally inserted statement before the
  new routine's declarations; this was removed. The first primitive test then
  measured velocity `1.0000014068850616e-8` instead of `1e-8`: the old velocity
  numerator added and subtracted equal large pressures separately. Grouping
  their difference as in the pinned source fixed this cancellation. Initial
  build/test and diagnostic logs are retained.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected execution, MPI, full matrix, compiled external oracle, convergence and
physical validation were not run. Compiler warnings remain.
