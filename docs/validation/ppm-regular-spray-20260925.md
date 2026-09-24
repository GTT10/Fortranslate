# Regular serial spray PPM — 2026-09-25

Source: base `9eb0421`, implementation `6ceeb2e`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3. Existing fixed tests-enabled
`build/spray-amr-reactive-debug` and selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug` were reused separately.
Artifacts: `build/ppm-regular-spray-checks/`; fixture:
`build/spray-amr-reactive-debug/tests/spray_ppm_mach_restart/spray-boundaries-_w6j3ra6`.

The shared serial spray application accepts
`hydro_reconstruction='characteristic_ppm'`. Its existing gas/particle split
already forwards this selector to `advance_reactive_full_3d`; no extra gas
stage is introduced. Existing checkpoint context binds the hydro selector.
The fixed and selected front ends share this application implementation.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^spray_ppm_mach_restart$' --output-on-failure`:
  **1/1 passed**, 12.16 s test / 12.21 s total.
- Fixed fixture uses Mach inflow/outflow with periodic transverse faces,
  molecular transport, Smagorinsky coefficient .17, parcel drag/heating/
  evaporation and chemistry disabled. Uninterrupted and resumed gas/parcel
  CSVs are byte-identical; final history ledger rows agree exactly.
  Existing relative mass/energy/momentum budget gates pass. Mach inflow
  retains interior composition rather than prescribed reservoir H2.
- Changing PPM to PLM on restart rejects before history output creation;
  the existing inlet temperature/Mach/pressure rejection checks also pass.
- `python3 -m py_compile tools/check_spray_boundaries.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No failed tests occurred; compiler warnings remain. Selected runtime,
regular MPI PPM, full configuration matrix, ignition, convergence and
physical spray validation were not run. This is a coupled numerical fixture,
not nozzle/diesel qualification or exact PeleC primitive-pressure parity.
