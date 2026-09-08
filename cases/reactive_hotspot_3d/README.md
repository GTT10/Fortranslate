# Reactive 3D hotspot

This periodic 6^3 case exercises elementary H2/O2 chemistry around the
general-EOS three-dimensional hydro step. The initial pressure is uniform,
the velocity is zero, and a Gaussian temperature hotspot changes the local
density through the mixture EOS. Each timestep uses the transactional
`R(dt/2)-H(dt)-R(dt/2)` split.

```sh
pelef_reactive_3d /path/to/cases/reactive_hotspot_3d/hotspot.nml
```

The case is a coupled numerical regression, not an external ignition
validation. Its checker enforces physical fields, EOS and species identities,
Euler conservation, H/O/N elemental conservation, and nonzero chemical
evolution.
