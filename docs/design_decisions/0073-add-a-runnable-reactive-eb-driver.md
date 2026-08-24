# Decision 0073: add a runnable reactive EB driver

## Context

The `0.80.0` EB path can complete one first-order reactive hydro update, but it
has no public configuration, time loop, timestep selection, diagnostics, or
output. A verified kernel is not yet a usable calculation workflow.

## Decision

Add a separate serial `pelef_reactive_eb_2d` application. Reuse the existing
reactive 2D thermodynamic configuration and initialization, and read a second
`embedded_boundary` namelist for a plane or circle level set and the weighted
StateRedist target volume fraction. Require the generated geometry to contain
at least one active and one cut cell.

Compute the hyperbolic CFL rate over active cells only, clip the last step to
the requested final time, and repeatedly call the complete first-order EB
hydro transaction. Integrate every conserved component with the cell volume
fraction and write geometry classifications and metrics beside the primitive
state in CSV.

The public path accepts only PCM hydro, zero-gradient outflow domain faces,
stationary slip-wall EB physics, and disabled chemistry and molecular
transport. Both the file reader and direct simulation API reject settings
outside that boundary.

## Consequences

PeleF can now run a two-dimensional multispecies general-EOS embedded-boundary
calculation from a committed input file to a checked output file. The active
CFL is not divided by the minimum volume fraction because weighted StateRedist
provides the small-cell stabilization. Higher-order reconstruction,
face-centroid interpolation, chemistry/transport splitting, other outer
boundaries, AMR, and MPI remain future work.
