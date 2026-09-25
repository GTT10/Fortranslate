# Single-level reactive EB examples

`plm_ctu.nml` runs a small plane-cut H2/O2 case with fourth-order PLM slope
selection, unsplit CTU hydro, StateRedist, chemistry and molecular transport.
From the repository root:

```sh
./build/spray-amr-reactive-debug/pelef_reactive_eb_3d cases/reactive_eb_3d/plm_ctu.nml
```

The input selects `hydro_method='plm_ctu'` and `riemann_solver='pelec'`.
Omitting `hydro_method` retains `pcm`. CTU uses multidirectional shock
flattening and the existing `difmag` coefficient (default 0.1). Physical
boundaries default to outflow hydro with closed transport faces. Extended geometry is built once
from plane, extruded-section or nodal input by nonperiodic domain-face
extension; state halos are rebuilt at each hydro step.

The selected executable accepts the same method with `thermo_model='selected'`
and explicit composition, as in the selected EB examples. Checkpoint restart
requires the same hydro method. The example is a numerical integration check, not a physically
validated combustion configuration. See the [focused evidence](../../docs/validation/eb-ctu-application-20260925.md).

## Prescribed inflow

`plm_ctu_inflow.nml` enables a prescribed inlet on the lower y face. Both the fixed
and selected single-level gas executables support these inputs:

- `boundary_faces(6)`: order x-low, x-high, y-low, y-high, z-low, z-high. Specify
  all six as supported physical labels for physical boundary transport. The default
  six `legacy` labels retain the previous boundary behavior; mixing modes rejects.
- `inflow_density(6)`, `inflow_pressure(6)`, `inflow_velocity(3,6)`: prescribed
  inlet density in kg/m³, pressure in Pa and Cartesian velocity in m/s. The inlet
  uses the initial mixture composition, including selected-mechanism composition.

Physical boundary inputs require `plm_ctu` and `pelec`. Ghost donors are rebuilt
at every transport stage and hydro update. Open-face donor signal speeds and
transport diffusivity limit the time step. Restart requires the same boundary
labels and inlet data. Periodic transport remains unsupported;
these inputs are not yet enabled in the spray or tree front ends.

See the [inflow/restart checks](../../docs/validation/eb-inflow-application-20260925.md).

## Fixed domain walls

`plm_ctu_wall.nml` uses six stationary isothermal domain faces. The supported wall
labels are `slip_wall`, `no_slip_wall` (adiabatic), and `isothermal_wall` (no-slip).
Set `wall_temperature(6)` in kelvin for isothermal faces; both transport and
thermal conduction must be enabled. The wall thermal rate also limits the timestep.
Restart requires the same wall temperatures. The embedded cut surface retains
its existing model; these controls apply to the six Cartesian domain faces.

Extended plane, extruded and nodal geometry reflects at fixed walls, including
repeated reflection when the six-cell halo exceeds the physical grid width.
Outflow/inflow faces retain extrapolated geometry. Moving walls and periodic
transport remain outside this input path. These are numerical integration cases,
not physical wall-flow validation. See the
[wall application checks](../../docs/validation/eb-wall-application-20260925.md).
