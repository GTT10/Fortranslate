# Decision 0095: checkpoint static three-level EB AMR

## Context

The public serial EB AMR application can advance a static root, middle, and
finest hierarchy, but its established checkpoint schemas represent either one
optional child or a set of sibling children. Neither format can encode a
nested child of a child without weakening its compatibility contract.

## Decision

Add a distinct versioned formatted checkpoint for the static three-level
mode. Store the complete conserved state and temperature field on every level,
both nested rectangles and refinement ratios, accepted time and step metadata,
base density, ordered species names, EB geometry, and numerical controls.

On restart, reconstruct all three EB geometries from the configured topology
and read the stream into private candidates. Require exact schema, mechanism,
geometry, topology, and numerical compatibility, validate the end marker, and
recover active temperatures from conserved state through the EOS before
publishing any level. Write scheduled checkpoints only after a committed
three-level Strang interval and retain the existing stop-after-write and final
checkpoint controls.

## Consequences

A static three-level run can stop and resume transactionally while preserving
accepted time, minimum timestep, and step cadence. Final time, maximum steps,
output paths, and checkpoint scheduling remain restart-mutable. The existing
single-patch and patch-set schemas remain byte-for-byte independent.

The format is serial and formatted, and the topology must match the configured
static hierarchy. Dynamic three-level regridding, topology-changing restart,
parallel I/O, arbitrary depth, and MPI-distributed ownership remain separate
work.
