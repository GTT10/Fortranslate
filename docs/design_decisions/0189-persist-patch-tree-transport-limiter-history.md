# Decision 0189: Persist patch-tree transport limiter history

## Context

Arbitrary-depth serial and sparse-MPI checkpoints retain the minimum accepted
root timestep, but not the cumulative minimum transport limiter theta. A
restarted application therefore resets theta to the neutral value `1` and
reports only the continuation suffix, even though its numerical state and
other clock metadata represent the complete logical run.

## Decision

Store minimum transport theta beside time, minimum timestep, step count, and
regrid count. Advance the base checkpoint envelope from schema 1 to 2 and the
public fingerprinted envelope from schema 4 to 5. Validate finite `[0,1]`
metadata, require the neutral value `1` for a zero-step checkpoint, and restore
the value before any continuation step.

Sparse I/O includes theta in write-control consensus and in the selected-root
restart metadata broadcast. Rejected reads publish theta `1` with the existing
empty distribution and sparse tree.

## Consequences

The reported limiter minimum now spans a serial or changed-rank sparse restart
boundary. Existing formatted checkpoints are intentionally rejected by the
strict schema check; numerical fields and the physics fingerprint layout do
not change.
