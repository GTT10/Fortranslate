# Physical-boundary PPM application input — 2026-09-25

Source: base `7d19d1a`, implementation `7dde6146ad5c074a318978b7cc5bbe71fdc303dc`.
GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3. Existing fixed tests-enabled
`build/spray-amr-reactive-debug` and selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug` were reused separately.

The shared configuration validator now permits physical-boundary PPM.
The reader's per-application capability flag remains required, so only the
regular fixed/selected gas front ends enable this selector. They already
forward boundary data and per-step inlet timestep limits to the physical
PPM library implementation. AMR/MPI application selection remains disabled.

`cases/reactive_hotspot_3d/ppm_mach.nml` supplies a 4^3 nonuniform-temperature
H2/O2/N2 fixture, low-x Mach inlet, high-x outflow and periodic y/z. Chemistry
is disabled; molecular viscosity, heat conduction, species diffusion and
barodiffusion are enabled. It does not exercise ignition or restart.

## Completed checks

Artifacts, generated fixed/selected inputs and exact checker:
`build/ppm-boundary-input-checks/`.

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2`: passed.
- `python3 build/ppm-boundary-input-checks/check.py`: both application runs
  passed, each completed **4 steps to 1e-6 s**. All 64 parsed output rows agree
  exactly. Fields are finite and density/temperature positive; mass changes
  under the open boundaries. Each reports maximum normalized all-component
  boundary-accounted residual **1.0733601507606494e-16**.
- `python3 -m py_compile build/ppm-boundary-input-checks/check.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, and `git diff --check`: passed.

No test failures occurred. Existing compiler warnings remain. No full matrix,
MPI/AMR PPM, pressure-target outlet application, restart, convergence or
external/physical validation was run in this increment. The conservative
transverse EOS update still differs from the pinned PeleC primitive-pressure
update; exact application agreement here is between PeleF mechanism paths.
