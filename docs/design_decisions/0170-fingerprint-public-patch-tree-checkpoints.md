# Decision 0170: Fingerprint public patch-tree checkpoints

## Context

Arbitrary-depth EB checkpoints validated schema, species, state width, stored
topology, and fields, but a public restart could silently change physics or
future regridding controls.

## Decision

Add schema 2 with a structured mesh, EB, physics, StateRedist, hierarchy, and
regridding fingerprint. Public serial and MPI applications always supply it.
Reject mismatches before reading numerical payloads. Exclude final time,
maximum steps, paths, checkpoint cadence, MPI ranks, and ownership weighting
so legitimate continuation and redistribution remain supported. Preserve
schema 1 only for explicit low-level compatibility callers.

## Consequences

Public restart now fails transactionally when an incompatible control such as
CFL changes, while rank-neutral restart and longer continuation remain valid.
The formatted checkpoint remains self-describing and human-inspectable.
