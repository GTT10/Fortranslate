# Regular 3D physical-boundary increment, 2026-09-15

Parent: `dda251f`. This increment adds per-face extrapolated outflow and
inviscid slip walls to regular fixed/selected 3D hydro/chemistry. The supported
inputs and remaining boundary work are in the
[case guide](../../cases/reactive_boundaries_3d/README.md).

The existing PCM/characteristic-PLM fluxes consume two freshly filled ghost
layers at each RK stage. Surface flux integration includes both RK stages;
application conservation checks account for this exchange. A failed candidate
leaves both state and boundary ledger uncommitted. Existing periodic dispatch
is retained. Physical-boundary molecular transport is rejected explicitly.

## Focused evidence

GNU Fortran 13.3.0, Debug in `build/port-spray-debug`:

```sh
cmake --build build/port-spray-debug --parallel 8 --target \
  test_reactive_boundaries_3d test_reactive_euler_3d_dimensional_reduction pelef_reactive_3d
ctest --test-dir build/port-spray-debug --output-on-failure \
  -R '^unit_reactive_(boundaries_3d|euler_3d_dimensional_reduction)$'
ctest --test-dir build/port-spray-debug --output-on-failure -R '^reactive_boundary_3d_'
```

- Two unit tests passed in 2.42 s: periodic dimensional reduction and the new
  boundary test (PCM/PLM, Rusanov/HLLC/PeleC, mixed faces, surface balance,
  wall mass/energy flux, pressure force, invalid configuration and rollback).
- Two application cases passed in 89.86 s. Normalized boundary-accounted
  Euler residuals were `1.083e-16` (open case) and `2.798e-17` (closed/reacting).
  The reacting case is expensive in the native Debug chemistry path; no
  additional broad matrix was run.
- The selected FFCM1_Red Release target built in `build/port-spray-cvode` and
  ran the mixed open case successfully. Input/output/log are retained in
  `boundary-selected-run/` there: `thermo_model='selected'`, three composition
  names CH3OH/O2/N2 with mole fractions 0.05/0.20/0.75; remaining inputs follow
  the committed `outflow.nml`. Chemistry is disabled in this selected smoke.
- Mechanism provenance, project contract and `git diff --check` passed.

Logs: `build/port-spray-debug/boundary-{build,tests,app-build,app-tests}.log`
and `build/port-spray-cvode/boundary-selected-build.log`. Individual application
stdout is in the Debug build's `Testing/Temporary/LastTest.log`.

These tests establish the implemented finite-volume boundary contracts, not
external PeleC field parity or complete project qualification. General AMR/EB,
physical-boundary transport, inflow/no-slip and particle coupling remain open.
