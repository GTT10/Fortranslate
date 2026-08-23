# Decision 0031: synchronize every multipatch split operator

## Context

The fixed two-level multipatch hydro path already advances each separated
child with its own flux register. Adding molecular transport and chemistry
creates intermediate states consumed by the next Strang-split operator. Those
states must represent the same composite solution on the parent and every
patch; deferring synchronization until the end would feed stale covered parent
cells into later coarse/fine ghost interpolation.

## Decision

Use the same conservative synchronization boundary after every spatial
operator. Each transport half interval advances the parent once with SSPRK2,
advances every patch for `r^2` substeps using parent-time-interpolated ghosts,
then refluxes one diffusive register per patch and averages down the full set.
Cell-local chemistry advances the parent and every patch over the same physical
half interval, then performs set-wide average-down. Hydro retains `r` fine
substeps and its own per-patch registers.

Compose one accepted interval as
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. The root timestep is the minimum of
the parent limit, `r` times every fine hyperbolic limit, and `r^2` times every
fine parabolic limit. Deep-copy the solution before the composition and restore
it if any database validation, chemistry solve, EOS recovery, flux update,
ghost fill, reflux, or average-down fails.

## Consequences

Every split boundary owns a synchronized composite state, and the existing
single-patch chemistry and transport kernels can be reused without changing
their numerical formulas. This performs more average-down and ghost-fill work
than a deferred synchronization scheme but keeps correctness explicit.
At the `0.39.0` milestone, dynamic patch-set rebuilds, arbitrary-depth
multipatch recursion, adjacent-box exchange, and distributed patch ownership
remained separate integrations. Decision 0032 subsequently connects the
two-level patch-set engine to the tag-driven application.
