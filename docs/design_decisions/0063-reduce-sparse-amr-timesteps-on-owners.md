# Decision 0063: reduce sparse AMR timesteps from owners

## Context

Sparse MPI AMR physics operates on globally single-copy patch payloads. The
serial timestep routine requires a materialized tree, so using it in a future
distributed driver would either gather every field or risk selecting a step
from an incomplete local view.

## Decision

Evaluate the qualified reactive CFL and molecular-transport timestep kernels
only for locally owned patches. Convert every level result to a root-step limit
with the same cumulative refinement-ratio scaling used by the serial patch
tree: `r` for hyperbolic work and `r^2` for parabolic work. Reduce the local
minima once with `MPI_MIN`.

Ranks without a local patch contribute `huge()`. Validate sparse ownership,
mechanism width, every evaluated patch state, and optional transport data
collectively. Publish a positive result only after all ranks accept; every
rejection returns `dt = 0` and `ok = .false.`.

## Consequences

The stable root timestep no longer requires a replicated field tree and is
identical to the serial hierarchy calculation for the same state. The routine
does not choose a stop-time-clipped step, advance the hierarchy, trigger
regrids, or manage output and restart cadence; those remain responsibilities
of a future runnable sparse MPI AMR driver.
