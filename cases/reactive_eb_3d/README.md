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
  all six as `inflow` or `outflow` for physical boundary transport. The default
  six `legacy` labels retain the previous boundary behavior; mixing modes rejects.
- `inflow_density(6)`, `inflow_pressure(6)`, `inflow_velocity(3,6)`: prescribed
  inlet density in kg/m³, pressure in Pa and Cartesian velocity in m/s. The inlet
  uses the initial mixture composition, including selected-mechanism composition.

Physical boundary inputs require `plm_ctu` and `pelec`. Ghost donors are rebuilt
at every transport stage and hydro update. Open-face donor signal speeds and
transport diffusivity limit the time step. Restart requires the same boundary
labels and inlet data. Periodic and wall transport policies remain unsupported;
these inputs are not yet enabled in the spray or tree front ends.

See the [inflow/restart checks](../../docs/validation/eb-inflow-application-20260925.md).
