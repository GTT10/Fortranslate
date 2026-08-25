# Decision 0088: support configured outflow-side EB patches

## Context

The reactive EB AMR hydro kernel required every fine rectangle to remain one
root cell away from each physical boundary. That allowed all fine exterior
states to be sampled from an adjacent root cell, but prevented a configured
static patch from refining a feature that reaches an outflow side. The EB flux
register already distinguishes actual coarse/fine interfaces from a patch side
coincident with the physical domain boundary.

## Decision

Permit the configured single-patch input and geometry helper to include root
index one or the final root index in either direction. Keep the established
coarse-time interpolation on each true coarse/fine side. On a side coincident
with an outflow physical boundary, construct the exterior conserved state and
temperature by copying the current fine boundary cell, which implements the
qualified zero-gradient outflow closure.

Pass the current fine candidate state and temperature explicitly into exterior
construction. Require this payload whenever any patch side is physical and
validate its shape, finiteness, and positive temperature. Continue to omit
flux-register accumulation and reflux where no coarse/fine interface exists.
Keep dynamic tag-planner candidates strictly internal until boundary-aware
buffering, clustering, and topology replacement are qualified separately.

## Consequences

A configured static two-level reactive EB hierarchy may now refine through an
outflow domain side without indexing a nonexistent root neighbor. Its physical
ghost state follows the evolving fine solution, while the other sides preserve
subcycled coarse-time interpolation and conservative synchronization.

This decision does not add inflow, wall, or periodic physical-side refinement;
dynamic physical-boundary regridding; deeper EB levels; molecular transport;
or distributed EB ownership.
