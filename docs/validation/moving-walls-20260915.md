# Tangential moving walls — 2026-09-15

Source: parent `6996930` plus this commit's regular 3D configuration, ghost fill,
transport/driver plumbing, case, tests and documentation. GNU Fortran 13.3 Debug,
`build/port-spray-debug`; focused local checks only.

`wall_velocity(3,6)` specifies a constant tangential velocity for adiabatic or
isothermal no-slip faces. Transport ghost momenta reflect about wall velocity;
total energy changes by the kinetic-energy difference, preserving internal
energy and temperature. Wall energy flux is viscous traction dotted with wall
velocity plus Fourier heat transfer. The existing surface ledger includes this
work through both RK stages and split transport half-steps. Euler walls remain
impermeable. Normal wall motion, nonfinite velocity and velocity on unsupported
face types are rejected; nonzero wall velocity requires molecular viscosity.

## Checks

- `unit_reactive_boundary_transport_3d` passed. For initially uniform tangential
  speed U=20 m/s and a top wall at V=50 m/s, the expected net x-momentum exchange
  is `2 mu (V-2U) A dt/dy`, and wall work is
  `2 mu (V-U) V A dt/dy`. Numerical surface exchanges match within relative
  1e-12; the volume energy change agrees within 1e-8. The same work check passes
  with the top wall isothermal at the initial gas temperature. Normal motion
  fails without changing the conserved state or temperature.
- `reactive_boundary_3d_no_slip`, `reactive_boundary_3d_isothermal` and the new
  `reactive_boundary_3d_moving_wall` passed (3/3). After the final guard and
  isothermal unit-test addition, the moving application and transport unit test
  were rerun: 2/2 passed in 4.44 seconds.
- In the moving-wall application, integrated energy exchange is
  `2.5189904286993581e-8` J. The maximum normalized all-component conservation
  residual after boundary accounting is `2.9096786113157094e-16`.
- Project contract, pinned FFCM1_Red provenance and `git diff --check` passed.
  No Python code changed. The initial build caught a misplaced local `speed`
  declaration; it was corrected before the successful builds and tests.

Commands:

```sh
cmake --build build/port-spray-debug --target test_reactive_boundary_transport_3d pelef_reactive_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^reactive_boundary_3d_(no_slip|isothermal|moving_wall)$' --output-on-failure
ctest --test-dir build/port-spray-debug -R '^(unit_reactive_boundary_transport_3d|reactive_boundary_3d_moving_wall)$' --output-on-failure
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
```

Application artifact:
`build/port-spray-debug/tests/boundary_3d_moving_wall/moving_wall.csv`.
The most recent logs are in `build/port-spray-debug/Testing/Temporary/LastTest.log`.
These tests establish flux/work accounting for a small regular grid, not
moving geometry, AMR/EB/particle coupling, an experimental lid-driven cavity
benchmark or full PeleC equivalence.
