# Reactive 3D entropy wave

This case exercises the regular-grid, periodic, general-EOS multispecies Euler
path. A constant-pressure H2/O2/N2 mixture carries a diagonal density and
temperature entropy wave through all three coordinate directions. Chemistry
and molecular transport are intentionally disabled so the 3D conservative
hydrodynamic boundary can be qualified independently.

Run from a writable directory:

```sh
pelef_reactive_3d /path/to/cases/reactive_entropy_wave_3d/entropy_wave.nml
```

The deterministic CSV is x-fastest and contains primitive variables, all
species mass fractions, and all species conserved densities.
