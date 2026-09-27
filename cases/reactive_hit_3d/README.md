# Regular 3D HIT velocity input

`pelef_reactive_3d` and its selected-mechanism counterpart accept an optional
`reactive_hit` namelist after the usual uniform-reactor input. This supplies
the initial velocity field for periodic turbulence calculations. LES requires
`PELEF_ENABLE_SPRAY=ON`; velocity-file initialization itself does not.

The example expects `hit_ic_4_32.dat` from PeleC commit
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, under `Exec/RegTests/HIT/`.
Its SHA-256 is
`983ebe8d4db570e56a634346871993e346ba6b3bb7eebf02ca8f39a6e5194776`.
Retrieve that blob from a local PeleC checkout containing the pinned commit:

```sh
git -C /path/to/PeleC show bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8:Exec/RegTests/HIT/hit_ic_4_32.dat > hit_ic_4_32.dat
/path/to/pelef_reactive_3d /path/to/pelec/cases/reactive_hit_3d/hit.nml
```

`input_file` is relative to the process working directory. ASCII files have
one `#` header followed by `inres**3` records of `x,y,z,u,v,w`, with x varying
fastest. `binfmt=.true.` instead reads the same records as native-endian
binary64 values, with no header or Fortran record markers. The coordinate
cube must use the same strictly increasing center array on each axis.

The reference's lookup, periodic repetition and trilinear interpolation are
retained, including extrapolation below the first center and wrapped upper
neighbors. Input velocities are scaled by
`mach_t0 * sound_speed / (sqrt(3) * uin_norm)` and added to the configured
mean velocity. Density, temperature and species retain their uniform initial
values; total energy includes the updated kinetic energy. Input coordinates
and application coordinates must use the same length units (PeleF uses SI).

This example uses NASA7 nitrogen and molecular transport. It does not reproduce
the reference's gamma-law EOS, constant viscosity derived from Taylor Reynolds
number, forced turbulence, or long-time statistics. The 16³ example is a
usable input template; focused checks used smaller grids. See
[implementation evidence](../../docs/validation/hit-input-20260927.md).

## Prescribed constant transport

For periodic SDC with `transport_enabled=.true.`, add:

```fortran
&reactive_constant_transport
 enabled=.true., viscosity=1.0e-5, conductivity=0.02
/
```

These are dynamic viscosity in Pa s and thermal conductivity in W/(m K).
Both must be finite and nonnegative. The same values reach old/new SDC fluxes
and the transport timestep bound, with or without LES. Omit the namelist to
retain the molecular mixture model. Species diffusion retains its mixture
model; for the HIT reference's zero species diffusion, set
`species_diffusion_enabled=.false.,barodiffusion_enabled=.false.` in
`reactive_3d`. Nonperiodic boundaries, split stepping, or disabled transport
reject enabled constant transport. Reynolds-based coefficient construction is
not automatic; supply the converted SI values explicitly.
