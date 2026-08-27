# Decision 0162: Store self-describing EB tree checkpoints

## Context

The arbitrary-depth reactive EB tree can change branching topology from
solution tags, so a fixed two- or three-level checkpoint cannot represent its
accepted state. Rebuilding geometry only from runtime configuration would also
couple restart to one level-set implementation and could change cut-cell
metrics across versions.

## Decision

Define a distinct versioned formatted checkpoint for the complete serial tree.
Store ordered species names, lifecycle metadata, every parent/child rectangle,
and all EB geometry metrics for every node. Store conserved state and the
accepted temperature field cell by cell.

Read into a private candidate. Bound levels, patch counts, and geometry cell
counts before allocation; validate species order, relation ordering, topology,
array extents, finite data, and the terminal marker. Recompute temperature from
the conserved general-EOS state before publishing the candidate.

## Consequences

An arbitrary-depth branching tree restarts independently of the geometry
builder that originally created it, while malformed, incompatible, truncated,
or thermodynamically invalid input publishes only neutral outputs. Formatted
files prioritize inspectability and deterministic round trips over compactness.
Sparse MPI root-only I/O and rank-neutral redistribution remain separate work.
