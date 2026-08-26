# Decision 0145: Reduce EB patch-tree CFL limits over all nodes

## Context

The `0.135.0` reactive EB patch tree owns state and temperature at arbitrary
depth and across branching relations, but it cannot yet select a stable root
interval. Fixed two- and three-level drivers already apply refinement-ratio
scaling to the qualified active-cell EB CFL limit. Making the low-level AMR
state module depend on the runnable driver would invert the module layering.

## Decision

Extract the existing single-node CFL calculation into
`reactive_eb_cfl_2d_mod`. Keep the public driver routine as a compatibility
wrapper and call the shared kernel directly from the patch-tree module.

Traverse every level and patch. Multiply each node-local CFL interval by the
cumulative product of the preceding relation refinement ratios, then take one
minimum in root time. Validate the tree, reactive component count, finite CFL,
every active-node conversion, and cumulative-scale arithmetic before
publishing. A fully covered node imposes no constraint; reject an entirely
inactive tree. Return zero on failure and treat the hierarchy as read-only.

## Consequences

The arbitrary-depth numerical tree can now select the hyperbolically stable
root interval without a fixed level count or duplicate CFL formula. Existing
single-level and fixed-depth APIs retain their names and behavior, with explicit
nonfinite-CFL rejection added to the shared kernel. This decision does not add
tree hydrodynamics, chemistry, transport, public time integration, dynamic
tagging, checkpoint I/O, or MPI ownership.
