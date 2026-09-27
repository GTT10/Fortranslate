# Heated embedded wall with heptane parcels

This bounded application case uses a 4³ grid over a 4 mm cube, a plane at
`x=0.7 mm`, initially 400 K gas, a 450 K embedded wall and four heptane parcels.
It takes two `1e-7 s` SDC steps with chemistry OFF and pinned Simple transport.
It demonstrates coupled wall heating and accounting; it is not the full pinned
Spray-EB example or a physically validated impinging spray.

Prepare the heptane bundle and selected MPI spray build as described in the
[Simple transport integration](../../docs/validation/simple-transport-integration-20260927.md).
Copy `tree.nml`, `spray.nml` and `level_set.dat` into a new run directory, then
run `pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml` there.
Outputs are protected against replacement.

`embedded_wall_temperature=0` (default) retains zero embedded molecular heat
flux. A positive value enables the quadratic heat stencil in the MPI SDC path;
transport and thermal conduction must be enabled. Common and recursive steps
support the setting, with or without particles. The serial application rejects
the enabled setting. The quadratic stencil requires equal cell widths and valid
fluid donors; unresolved stencils are rejected rather than replaced.

The app caps every interval by `transport_cfl / heat_frequency`, including when
`automatic_timestep=.false.`. This frozen-conductivity bound complements existing
state acceptance/retry checks. It is not a general nonlinear stability guarantee.
The embedded wall temperature is stored in checkpoints; changed, newly enabled
or newly disabled settings are rejected during restart. Existing checkpoints
without this setting remain readable when it is disabled.

Wall heat appears in the energy row of `fortran.csv.budget.csv`; particle-wall
impulse and physical domain boundaries retain their separate accounting.
The supplied input enables molecular wall heat only. Add
`embedded_noslip=.true.` to `eb_tree` to enable stationary-wall viscous traction;
SDC, transport and viscosity must be enabled. It works independently of wall
heat: set `embedded_wall_temperature=0` for an adiabatic no-slip wall. Both
common and recursive MPI gas/particle paths support it; serial rejects it.
The interval cap includes the frozen viscous frequency, and checkpoints bind
the no-slip setting. Stationary traction supplies wall momentum without wall
mechanical energy. Full pinned Spray-EB execution and physical validation
remain open.
