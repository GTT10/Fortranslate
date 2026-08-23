# Decision 0026: advance reactive AMR through adjacent-level recursion

## Status

Accepted for PeleF `0.31.0`.

## Context

The arbitrary-depth hierarchy stores a runtime chain of already-qualified
adjacent coarse/fine relations. Reactive hydro, chemistry, and transport were
previously coupled only in the separate two-level application. Adding a fixed
three-level driver would repeat that limit and make flux synchronization order
hard to reason about.

## Decision

Store conserved state, temperature, and ghost cells in an allocatable level
array. Advance one level, recursively advance its child over the same physical
interval, and synchronize the current relation before returning. Use `r`
child calls for hyperbolic hydro and `r^2` calls for explicit parabolic
transport. Return time-integrated outer fluxes from each recursive call so its
parent owns exactly one flux register for that relation.

Advance chemistry on every level and average down deepest first. Compose the
operators as reaction--transport--hydro--transport--reaction and deep-copy the
complete solution before the step for transactional rollback.

## Consequences

- reactive advancement has no compile-time level limit;
- each recursive frame owns only one adjacent interface and flux register;
- synchronization naturally proceeds from deepest level to root;
- the global timestep includes cumulative hyperbolic and parabolic scaling;
- the existing dynamic executable remains the qualified two-level regrid path;
- recursive tag generation, regridding, output, multiple patches,
  physical-boundary refinement, and MPI patch ownership remain future work.
