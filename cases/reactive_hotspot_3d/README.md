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

## PPM

The fixed and selected regular gas executables accept
`reconstruction = "characteristic_ppm"` with periodic or supported physical
boundaries. This selects
the conservative twelve-solve CTU step with flattening enabled and contact
steepening disabled. The `pelec` solver uses separate pressure tracing and
artificial viscosity; the regular gas application has no checkpoint/restart interface.
See the [focused input checks](../../docs/validation/ppm-input-3d-20260925.md).

`ppm_mach.nml` is a small nonreacting H2/O2/N2 example with molecular transport,
Mach inflow on low x, outflow on high x, and periodic y/z:

```sh
pelef_reactive_3d cases/reactive_hotspot_3d/ppm_mach.nml
```

For a selected H2/O2 build, set `thermo_model='selected'` and supply
`composition_count=3`, `composition_species='H2','O2','N2'`, and
`composition_mole_fractions=0.3,0.15,0.55`. This is an application coupling
fixture, not a nozzle validation. See the
[physical-boundary input checks](../../docs/validation/ppm-boundary-input-20260925.md).

### Artificial viscosity

`difmag=-1` (default) selects the reference value `0.1` for PCM, PLM and PPM,
independent of solver. Set `difmag=0` to disable it, or a positive value to
choose its strength. Fixed and selected regular gas executables print the
resolved coefficient. Older zero-policy runs/checkpoints require explicit zero.
See [default-policy checks](../../docs/validation/regular-gas-default-difmag-20260925.md).

## Optional regular-grid LES with SDC

Build with `PELEF_ENABLE_SPRAY=ON`. Fixed and selected regular 3D executables
accept this additional namelist together with `time_integrator='sdc'` in
`reactive_3d`:

```fortran
&reactive_les
 enabled=.true., model='smagorinsky', cs=.17, prt=.7, sct=.7, ci=0
/
```

Models are `smagorinsky`, `wale`, `vreman`, and `dynamic`. Dynamic uses
`test_model='box'`, `test_ratio=2`, `grid_ratio=1`, and
`harmonic_mean=.false.` by default; it requires `ci=0`. Nondefault filter
controls require Dynamic. The namelist defaults to disabled. Enabled LES
requires SDC and the optional runtime; unsupported builds reject it.
`sdc_reactor='cvode'` additionally requires SUNDIALS and enabled chemistry.
Boundary flux inventories and the LES diffusion timestep bound are included.
See [focused application evidence](../../docs/validation/les-regular-app-20260927.md)
for the tested scope. These are experimental LES capabilities.
