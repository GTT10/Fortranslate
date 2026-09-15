# Prescribed-temperature stationary walls in regular 3D

Parent: `f0a8ccc`. `isothermal_wall` imposes a stationary no-slip wall with
`wall_temperature(6)`. Molecular heat conduction is required. Wall conductivity
uses the specified temperature and adjacent-cell composition/pressure;
Fourier heat input is `k_wall*(T_wall-T_cell)/(dx_normal/2)`.

The physical face energy flux is set directly. There is no extrapolated
negative ghost temperature. Hydro reflection and viscous wall shear retain
the no-slip implementation; wall heat joins the split-step boundary ledger.
The thermal timestep sums wall conductances divided by adjacent-cell heat
capacity at corners, with half the interior diffusion allowance reserved.
Mechanism temperature bounds and candidate EOS validity are checked.

## Focused evidence

GNU Fortran 13.3.0 Debug, `build/port-spray-debug`:

```sh
cmake --build build/port-spray-debug --parallel 8 --target \
  pelef_reactive_3d test_reactive_boundary_transport_3d test_reactive_boundaries_3d
ctest --test-dir build/port-spray-debug --output-on-failure \
  -R '^(unit_reactive_(boundaries_3d|boundary_transport_3d)|reactive_boundary_3d_isothermal)$'
```

Three tests passed in 8.40 s: analytical Fourier heat flux and thermal timestep
on all six faces with hot/cold walls, energy uptake, invalid-wall rollback,
configuration requirements, existing boundary/transport behavior and the
coupled wall example. Its boundary-accounted Euler residual was `1.062e-16`
after gas energy decreased by `3.918e-7 J`.

The existing FFCM1_Red Release configuration built and ran the selected front
end. Input/output/stdout are in `build/port-spray-cvode/isothermal-run/`:
CH3OH/O2/N2 mole fractions 0.05/0.20/0.75 and otherwise the committed
`isothermal.nml`. Its Euler residual was `3.402e-17`, energy change `-1.800e-7 J`.
Chemistry is disabled in these application examples.

Logs use prefix `isothermal-` in both build directories. Mechanism provenance,
project contract and whitespace checks passed. No Python source changed;
no broad matrix or external heat-transfer benchmark was run.

Moving walls, characteristic boundaries, general 3D AMR/EB and particle
boundary coupling remain open. This is an implementation increment, not
project completion; see the [boundary guide](../../cases/reactive_boundaries_3d/README.md).
