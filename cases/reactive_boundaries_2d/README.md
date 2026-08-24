# Reactive physical-boundary cases

- `couette.nml`: periodic x direction, two no-slip adiabatic walls, moving upper wall.
- `thermal_channel.nml`: periodic x direction, slip isothermal walls at 800 K and 1200 K.
- `inflow_outflow.nml`: fixed-state x inflow, zero-gradient x outflow, periodic y direction.
- `prescribed_species_wall.nml`: lower-wall H2/O2 conversion flux with zero
  net mass and species-enthalpy transport.

These are qualification cases for the 0.18.0 physical-boundary and 0.69.0
prescribed-species wall-flux milestones.
