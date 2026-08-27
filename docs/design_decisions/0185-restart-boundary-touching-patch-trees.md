# Decision 0185: Restart boundary-touching patch trees

## Context

Arbitrary-depth outflow-boundary children are qualified during a fresh serial
or sparse-MPI run. The checkpoint format stores every child bound and geometry,
but the public restart case still used only interior children. A process or
rank-count boundary therefore had no gate proving that physical-side contact
survives reconstruction and ownership redistribution.

## Decision

Move the existing public patch-tree split-run hotspot to the x-upper physical
boundary. Require the uninterrupted, stopped, and restarted composite outputs
to contain four levels and reach that exact side at every level. Preserve
serial field parity and the two-rank checkpoint to four-/eight-rank restart
with changed ownership weight.

Keep checkpoint schema 4. Stored child bounds and self-describing geometry
already encode physical-side contact; adding a derived boundary flag would
create a second source of truth.

## Consequences

The qualified outflow-boundary topology can cross independent processes and
sparse ownership maps without a schema change. Non-outflow physical children
and periodic-seam topology remain unsupported by the EB configuration.
