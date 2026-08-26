# Decision 0120: schedule tagged sparse EB regrids

## Context

The sparse MPI EB hierarchy can advance to a target time and can rebuild an
explicitly supplied child topology, but applications still have to materialize
root fields and coordinate those operations outside the public clock. EB
geometry also does not retain the level-set source needed to construct an
arbitrary new fine rectangle.

## Decision

Add a collective tag-driven sparse regrid entrypoint. Gather root tiles only to
the root physics owner, evaluate the established temperature-gradient planner
there, and broadcast only the compact ordered patch-plan metadata. Require the
caller to supply a geometry-builder callback that reconstructs each planned EB
child from its coarse bounds and refinement ratio. Feed those geometries into
the transactional explicit regrid from Decision 0119.

Extend the public sparse target-time loop with an optional complete regrid
control bundle: accepted-step interval, criteria, refinement ratio, geometry
builder, evaluation count, and topology-change count. On a due step, advance a
private physics candidate and regrid that candidate before publishing either
one. Commit state, topology, clock, counters, limiter, and traffic only when the
combined operation succeeds.

## Consequences

The public two-level sparse EB clock can now create, move, resize, or remove
fine children from solution tags without exposing a replicated numerical field
to its caller. Invalid or rank-inconsistent controls reject before advancement,
and a late geometry or regrid failure discards the otherwise valid step.

Planning is owner-local, but geometry construction intentionally runs on every
rank and the regrid still uses the replicated compatibility window introduced
in `0.127.0`. Replica-free old/new overlap migration, sparse checkpoint/output,
arbitrary-depth dynamic EB topology, and decomposed root physics remain future
work.
