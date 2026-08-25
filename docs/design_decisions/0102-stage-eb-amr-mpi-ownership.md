# Decision 0102: stage two-dimensional EB AMR MPI ownership

## Context

The serial reactive EB AMR hierarchy now supports one patch, separated
siblings, and strict three-level nesting with hydrodynamics, chemistry, and
molecular transport. Direct distributed advancement requires a deterministic
answer to which rank owns each root region and fine patch before fields,
coarse/fine traffic, flux registers, and regridding can become rank-local.

## Decision

Keep the serial EB geometry and patch-set descriptors replicated. Partition
the root into contiguous, nonoverlapping y-tiles, using up to one tile per
rank, and treat every fine sibling as one additional work entity. Assign all
entities deterministically with a greedy 64-bit work schedule. Root work is
its Cartesian cell count. Fine work is its cell count multiplied by the
refinement ratio raised to exponent zero, one, or two for storage,
hyperbolic-subcycle, or parabolic-subcycle weighting.

For this first bridge, field replicas remain allocated on every rank while
only the recorded owner payload is authoritative. Synchronize each root tile
and each child state and temperature from its owner after proving collective
agreement on the root and child topology. Perform all validation before any
communication that depends on the owner map, and publish synchronized outputs
only after communicator-wide acceptance.

## Consequences

Two-dimensional EB AMR now has a qualified MPI ownership boundary, 64-bit
load accounting, collective work-model rejection, and transactional
owner-authoritative synchronization on one, two, four, and eight ranks. The
root tiling gives every tested rank a concrete ownership unit even when the
fine patch count is smaller than the communicator.

This remains a correctness-first replicated bridge. Rank-local sparse field
storage, direct EB hydro/chemistry/transport execution, coarse/fine and halo
point-to-point traffic, distributed flux registers, topology migration,
parallel checkpoint I/O, and rank-count field parity for a complete time loop
remain outside this milestone.
