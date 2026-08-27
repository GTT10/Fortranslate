# Reactive EB wall-transport regression

`reference.nml` advances the plane-EB hotspot without molecular transport.
`transport.nml` enables bulk viscosity and conduction, selects a 1500 K
isothermal no-slip embedded wall, and moves it at 5 m/s tangential to the plane.

The comparison requires bulk conduction to reduce the active temperature span,
the hot wall to increase at least one cut-cell temperature, and viscous wall
traction to create a nonzero tangential cut-cell velocity. Both outputs must
retain finite positive active thermodynamics and species closure.

