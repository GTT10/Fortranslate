# Reacting hotspot in 2D

A Gaussian temperature hotspot is initialized in a stoichiometric H2/O2/N2
mixture with a small radical seed. Cell-local constant-volume chemistry is
Strang split around the unsplit two-dimensional CTU hydro update.

- `hotspot.nml` uses characteristic PLM.
- `hotspot_characteristic_ppm.nml` uses time-traced characteristic PPM with
  bounded contact steepening and shock flattening enabled.

Both inputs are structural reactive-flow regressions rather than ignition-delay
validation cases for a complete mechanism.
