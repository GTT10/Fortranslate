# Decision 0035: synchronize each child set inside recursive tree advancement

## Context

The patch-tree geometry and nested flux registers can represent separated
children at arbitrary depth, but a single deepest-to-root synchronization pass
is not sufficient to determine when each child state may use an updated parent
ghost value. Hyperbolic subcycling requires every parent to retain its state at
both ends of the parent interval while its children advance.

## Decision

Advance the reactive patch tree with a depth-first recursive routine. Each call
advances one parent patch once, stores its start and end states, and accumulates
the coarse boundary flux for every local child. It then fills all child ghosts
from the time-interpolated parent states and recursively advances every child
for each refinement-ratio substep. Returned time-integrated child boundary
fluxes accumulate in the matching parent-owned registers.

After every local child finishes its subcycles, reflux and average down the
complete child set into that parent. A public advance deep-copies the whole
solution first and restores it if any nested operation fails. Initially accept
only strictly interior child patches; this avoids claiming physical-boundary or
periodic-seam behavior before it has a dedicated gate.

## Consequences

Branches may stop at different depths while other branches continue refining,
and every descendant correction reaches the root in one accepted interval.
The cumulative refinement ratios define both the child call counts and the
root-level CFL reduction. The first qualification uses PCM hydro. Chemistry,
transport, dynamic tree rebuilds, same-level exchange, physical-boundary
children, and distributed ownership remain explicit follow-on work.
