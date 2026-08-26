# Decision 0146: Advance EB patch-tree hydrodynamics recursively

## Context

The `0.153.0` reactive EB numerical tree can select a stable root interval at
arbitrary depth, but only fixed two-/three-level and sibling patch-set types can
advance 2D EB hydrodynamics. Flattening the runtime tree into those APIs would
lose actual parent links, duplicate intermediate fields, and reintroduce a
depth bound.

## Decision

Advance one node per recursive call. Retain its start and uncorrected end states
for child exterior interpolation, create one EB flux register per ordered
child, and invoke each child recursively for the relation refinement ratio.
Accumulate parent and child fluxes with their actual time weights, then reflux
and average down each child into its actual parent.

Close every refined subtree against flux through the parent node's outer EB
boundary. Apply any density, total-energy, and species residual only to active
parent cells not represented by children, enforce density/species consistency,
and recover temperature through the EOS. Finish with deepest-first hierarchy
synchronization and commit only a completely valid candidate. Publish optional
per-level advance counts only with that commit.

## Consequences

Runtime depth and branching now determine the hydro schedule without a fixed
level count. Nested and sibling relations share the qualified level, exterior,
flux-register, reflux, and average-down kernels. Failure at any depth preserves
the accepted hierarchy exactly. This decision does not add chemistry,
molecular transport, a public time loop, dynamic tagging, checkpoint I/O, or
MPI ownership for the new tree.
