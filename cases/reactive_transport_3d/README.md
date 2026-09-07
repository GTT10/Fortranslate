# Reactive 3D molecular-transport hotspot

This paired periodic case isolates regular-grid molecular transport from the
three-dimensional general-EOS hydro update. Both inputs use the same inert
full-H2O2 hotspot. `control.nml` disables transport; `transport.nml` enables
Newtonian viscosity, Fourier conduction, mixture-averaged species diffusion,
barodiffusion, correction velocity, and species-enthalpy energy transport.

The case is a numerical conservation and coupling regression, not an external
physical validation. The unit and regression suites separately qualify active
viscous, thermal, and species-diffusion modes; this pair qualifies the public
driver, configuration, full transport database, timestep selection, and CSV
contract.

```sh
pelef_reactive_3d /path/to/cases/reactive_transport_3d/control.nml
pelef_reactive_3d /path/to/cases/reactive_transport_3d/transport.nml
```
