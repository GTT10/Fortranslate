# Decision 0091: recurse EB hydro across regular nested interfaces

## Context

Three-level synchronization establishes ownership but does not advance the
hierarchy. Recursive hydrodynamics needs two temporal interpolation contexts,
two flux-register lifetimes, and a deterministic commit order. Directly
reusing the two-level EB register at a finest interface crossed by the embedded
boundary exposed a nonconservative geometric case, so that capability cannot
be admitted implicitly.

## Decision

Advance the root once. For every root-to-middle refinement substep, advance
the middle from a root-time-interpolated exterior state and then advance the
finest through its complete middle-to-finest subcycle. Close the inner flux
register and average down finest to middle before starting the next middle
substep. After all middle substeps, close the outer register and synchronize
the complete hierarchy deepest first.

Require the finest rectangle to remain two middle cells from the middle
boundary. Require every face on the finest coarse/fine interface to have unit
open-area fraction. Reject a fractional interface transactionally before any
hydro work. Root and middle may continue to contain cut and covered EB cells.

## Consequences

The qualified three-level kernel has exact update counts, independent nested
flux histories, conservative mass/energy/species synchronization, EOS-positive
parents, and whole-hierarchy rollback. The interface restriction makes the
current numerical boundary explicit instead of accepting a measured
conservation defect.

Supporting an EB-cut nested interface requires a dedicated geometric
coarse/fine flux and re-reflux treatment. Chemistry splitting, dynamic
three-level ownership, checkpoints, public output, transport, arbitrary depth,
and MPI distribution also remain separate milestones.
