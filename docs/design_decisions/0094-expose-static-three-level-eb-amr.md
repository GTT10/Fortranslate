# Decision 0094: expose static three-level EB AMR

## Context

The three-level transfer, recursive hydro, EB-cut conservation closure, and
Strang chemistry transactions were directly callable but not reachable from a
namelist-driven executable. Consequently they qualified numerical kernels but
did not yet form a runnable application lifecycle.

## Decision

Extend the existing serial reactive EB AMR configuration with an explicit
`three_level_enabled` mode. Retain the existing root-index rectangle as the
middle patch and add one rectangle expressed in middle indices for the
finest patch. Require the finest rectangle to retain the qualified two-cell
margin from every middle boundary.

Initialize root state normally, then apply PCM prolongation from root to
middle and middle to finest. Select each root interval from the minimum of the
root limit, `r1` times the middle limit, and `r1*r2` times the finest limit.
Advance the complete three-level Strang transaction and write separate
geometry-aware CSV files only after the requested final time is reached.

## Consequences

The public executable now runs the static three-level reactive EB hierarchy,
including a finest coarse/fine interface cut by the embedded boundary. The
input mode is mutually exclusive with multipatch and dynamic regridding, and
it rejects checkpoint/restart because the established schemas encode only
root plus one patch or root plus sibling patches.

Dynamic three-level topology, a new checkpoint schema, molecular transport,
arbitrary-depth EB recursion, and MPI-distributed ownership remain separate
work.
