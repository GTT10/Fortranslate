# Decision 0021: conservative dynamic AMR regridding in 1D

## Status

Accepted for PeleF `0.26.0`.

## Context

The static hierarchy established transfer and reflux rules but could not place
refinement from the evolving solution. The next step needs deterministic
tagging and state migration without prematurely introducing general AMReX-style
box clustering or boundary-patch ghost filling.

## Decision

Tag one configured state component using the maximum adjacent jump normalized
by a local magnitude floor. Also require an absolute jump threshold. Buffer the
bounding interval of all tags and expand it to a configured minimum width,
producing at most one strictly interior fine patch.

On regrid, average the old fine patch down before it is discarded, prolong the
new patch conservatively from the synchronized coarse state, and copy the exact
old fine representation wherever old and new patches overlap at the same
refinement ratio. Treat an empty tag set as fine-level removal. Reject a tagged
physical-boundary cell until boundary-patch ghost rules exist.

## Consequences

- creation, movement, resizing, and removal preserve composite integrals;
- retained fine cells do not lose subcell information during patch movement;
- criteria and transfers remain independent of the fluid-state width;
- widely separated tags may create an intentionally broad single patch;
- multiple patches, physical-boundary refinement, multilevel recursion,
  solver scheduling, and MPI ownership remain follow-on work.
