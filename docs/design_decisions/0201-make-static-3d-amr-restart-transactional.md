# Decision 0201: make static 3D AMR restart transactional

## Context

Milestone `0.207.0` qualified a serial, static, synchronized two-level 3D AMR
hydrodynamic hierarchy. Restart is the next prerequisite for longer runs and
distributed ownership, but restoring arrays without proving that their
thermodynamics, mesh, solver, and coarse/fine relationship match the current
run can silently create a different problem. A partial read must not damage a
valid in-memory hierarchy.

## Decision

Write checkpoints only after a completed root step, reflux, average-down, and
temperature recovery. Use an explicit magic record and schema version. Store
the full NASA7 species definitions, solver and mesh fingerprint, both level
states and temperatures, time, cumulative step count, original conservation
baseline, and reflux history.

On restart, compare the complete fingerprint before allocating payload
candidates. Read and validate every payload record privately, recover both
temperature fields independently, require stored/recovered agreement and
covered-parent synchronization, and reject trailing data. Commit arrays and
metadata to the caller only after the closing record and end-of-file are
validated. Permit a later requested final time while requiring the checkpoint
time and step to fit inside the new stopping limits.

Qualify continuation by stopping the public full-H2O2 case at an interior
root-step boundary and requiring both final CSV files to match the
uninterrupted run byte for byte.

## Consequences

The qualified static hierarchy can be interrupted and resumed without
changing its numerical result. Missing, corrupt, physically inconsistent, or
configuration-incompatible files fail before caller state is published. The
schema is deliberately strict and currently same-build, serial, hydro-only,
and single-patch. Distributed layout, topology reconstruction, AMR source
history, EB data, checksummed/atomic production storage, and schema migration
remain separate decisions.
