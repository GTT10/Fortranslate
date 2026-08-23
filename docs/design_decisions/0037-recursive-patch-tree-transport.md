# Decision 0037: reuse nested registers with parabolic tree subcycling

## Context

Molecular transport is conservative and needs the same parent-child flux
correction as hydro, but its explicit stability limit scales with `dx^2`.
Applying only refinement-ratio hyperbolic subcycling would make fine transport
updates unstable, while advancing all patches at the finest global interval
would discard the existing recursive ownership model.

## Decision

Advance transport through a depth-first recursion parallel to hydro. Advance a
parent once for its physical interval, then advance each child `r^2` times.
Fill child ghosts from the parent start/end states at each fine interval
midpoint. Accumulate time-integrated coarse and fine diffusive boundary fluxes
in the matching parent-owned child register, then reflux and average down the
complete local child set before returning.

Track transport calls independently from hydro calls. Extend the stable root
timestep with every patch transport limit multiplied by the square of its
cumulative refinement ratio. Require the transport database whenever transport
is enabled, before mutating the solution.

## Consequences

Every separated branch follows its own parabolic recursion and a leaf
correction propagates through all ancestors. The full public operator becomes
`R-T-H-T-R` and retains whole-tree rollback. The number of calls grows quickly
with depth, which is expected for explicit parabolic subcycling and remains a
future optimization target alongside dynamic trees and distributed ownership.
