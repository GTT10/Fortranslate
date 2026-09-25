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
boundaries in this application are outflow. Extended geometry is built once
from plane, extruded-section or nodal input by nonperiodic domain-face
extension; state halos are rebuilt at each hydro step.

The selected executable accepts the same method with `thermo_model='selected'`
and explicit composition, as in the selected EB examples. Checkpoint restart
requires the same hydro method. EB tree/MPI/spray applications do not yet accept
this selector. The example is a numerical integration check, not a physically
validated combustion configuration. See the [focused evidence](../../docs/validation/eb-ctu-application-20260925.md).
