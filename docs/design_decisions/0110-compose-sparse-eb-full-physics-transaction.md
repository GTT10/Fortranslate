# Decision 0110: compose a sparse EB full-physics transaction

## Context

Sparse chemistry and average-down are direct, while the qualified hydro and
molecular-transport operators still require complete temporary hierarchies. A
caller should nevertheless be able to retain owner-only state across a full
reactive split and observe one transaction boundary.

## Decision

Advance the first chemistry half-step directly on a sparse candidate, then
materialize once for the central transport-hydro-transport sequence. Keep that
complete candidate through both transport stages, scatter it back to sparse
owners, and execute the final chemistry half-step directly.

Retain the caller's sparse object until every stage succeeds. Publish chemistry,
hydro, and transport counts plus the limiter minimum only with the final sparse
commit.

## Consequences

The full `R-T-H-T-R` entrypoint has sparse persistent input and output, exact
owner accounting, serial parity, and outer rollback after late failures. Only
one complete hierarchy window is created instead of crossing the boundary for
each operator.

Hydro and transport numerical work still consumes replicated fields. Direct
sparse stencils and flux synchronization, point-to-point traffic, public time
loop control, regridding, checkpointing, and output remain future work.
