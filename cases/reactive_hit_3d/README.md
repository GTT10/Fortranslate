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
values; total energy includes the updated kinetic energy. `coordinate_scale`
multiplies the file's coordinates before interpolation (default 1). Set it to
`.01` for centimetre input and supply the application domain in metres.
It must be finite and positive; it does not multiply velocity samples.

This example uses NASA7 nitrogen and molecular transport. It does not reproduce
the reference's gamma-law EOS, forced turbulence, or long-time statistics.
Constant viscosity from Taylor Reynolds number is available through the input
below. The 16³ example is a
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
reject enabled constant transport. Alternatively, derive them from the initialized HIT state:

```fortran
&reactive_constant_transport
 enabled=.true., mode='hit', lambda0=.5, reynolds_lambda0=100, prandtl=.71
/
```

Use one constant-transport namelist. `lambda0` is in metres. Density, cp and
`urms=mach_t0*sound_speed/sqrt(3)` come from successful HIT initialization;
`viscosity=rho*urms*lambda0/reynolds_lambda0` and
`conductivity=viscosity*cp/prandtl`. This RMS is independent of `uin_norm`.
Positive length/Re/Pr and positive initialized scales are required. Explicit
nonzero coefficients cannot be mixed with `mode='hit'`. The default mode is
`specified`. See [coefficient checks](../../docs/validation/hit-transport-scales-20260927.md).

## Calorically perfect reference gas

Use [gamma_law.nml](gamma_law.nml) with `pelef_reactive_3d` for a one-component
pseudo-air model. `thermo_model='gamma_law'` selects a constant cp/cv ratio;
`&gamma_law gamma=1.4 /` sets it explicitly. Molecular weight is 28.97 kg/kmol,
formation energy is zero, and the supported temperature interval is 1–1e6 K.
This model requires nonreacting SDC and no named composition. Transport, when
enabled, must use explicit constant transport with species diffusion disabled.
The selected-mechanism executable retains its own mechanism.

The example converts the reference's 2π cm cube through `coordinate_scale=.01`,
uses a 0.06283185307179586 m domain and converts the 0.5 cm Taylor length to
0.005 m. Time remains in seconds. The 16³ template has not been qualified as
a reference trajectory. See [gamma-law checks](../../docs/validation/gamma-law-app-20260927.md)
and [coordinate conversion checks](../../docs/validation/hit-coordinate-20260927.md).

## Full reference executable preparation

Export the pinned source and recursive dependencies from a populated local
PeleC checkout. The destination must be new. This leaves that checkout intact:

```sh
python3 tools/prepare_pelec_reference.py --reference /path/to/PeleC --destination build/pelec-pinned-hit-source
env -u LD_LIBRARY_PATH -u LIBRARY_PATH GIT_CEILING_DIRECTORIES="$PWD/build" cmake -S build/pelec-pinned-hit-source -B build/pelec-pinned-hit -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DPELE_ENABLE_MPI=OFF -DPELE_ENABLE_PARTICLES=OFF -DPELE_ENABLE_FCOMPARE=OFF -DBUILD_TESTING=OFF
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/pelec-pinned-hit --target PeleC-HIT --parallel 2
```

The exported `pelef-reference-provenance.json` records source identities.
`GIT_CEILING_DIRECTORIES` keeps upstream CMake from accidentally reporting the
enclosing PeleF repository's commit as the reference identity; its embedded
Git labels are `unknown` for these archive exports.

[pelec_reference.inp](pelec_reference.inp) prepares a 4³, one-step, hydro-only
reference run with PPM and two SDC iterations. Supply the pinned velocity file
in its working directory. It uses CGS units and `1e-8` seconds. This preparation
does not establish evolved-field parity; the current completed comparison is
the [initialization kernel check](../../docs/validation/hit-initialization-reference-20260927.md).
