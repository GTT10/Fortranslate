# Decision 0148: Advance EB patch-tree transport recursively

## Context

The arbitrary-depth numerical EB tree owns CFL selection, recursive
hydrodynamics, and active-cell chemistry, while molecular transport remains
limited to fixed two-/three-level and sibling patch-set APIs. Applying a
fixed-depth operation relation by relation would repeat ancestor Euler stages
and would not provide one SSPRK2 rollback boundary for a branching tree.

## Decision

Advance one transport Euler node per recursive call with the established EB
diffusive-flux, redistribution, exterior-fill, flux-register, reflux, and
temperature-recovery kernels. Retain the parent start and Euler-end fields,
subcycle every ordered child by the relation refinement ratio, and close every
refined subtree against the parent's outer diffusive flux before deepest-first
synchronization.

Build SSPRK2 from two complete recursive Euler candidates. After the second
Euler tree, blend every runtime node with its accepted stage-zero state,
recover active temperatures through the EOS, synchronize deepest first, and
publish the state, minimum limiter theta, and optional per-level node counts
only after final validation.

## Consequences

Molecular transport now follows arbitrary runtime depth and branching without
a fixed level count. A three-level runtime chain retains field and temperature
parity with the established fixed-depth SSPRK2 path. A separate four-level
branching gate verifies exact recursive scheduling, composite conservation,
changed state, thermodynamics, and rollback.

This decision does not compose chemistry, transport, and hydrodynamics into
one `R-T-H-T-R` transaction and does not add a public time loop, dynamic
tagging, checkpoint I/O, or MPI ownership for the numerical tree.
