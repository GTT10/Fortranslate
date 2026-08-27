# Decision 0161: Regrid EB tree tags on owners

## Context

The serial arbitrary-depth EB tree can derive a complete branching plan from
temperature gradients, but the sparse MPI tree stores every numerical field on
one owner only. Materializing that tree to plan or rebuild topology would
reintroduce an all-rank field replica at a normal lifecycle boundary.

## Decision

Evaluate each prospective parent only on its current owner. Reduce compact tag
bounds and counts, construct the same parent-major plan on every rank, and ask
the caller to build identical EB geometry collectively. Use the deterministic
work-weighted distribution for the candidate topology.

Initialize new children by sending PCM parent data directly from old owners to
new owners. Retain geometrically identical same-resolution overlap with direct
rectangle transfers, then synchronize the candidate deepest-first and verify
its composite conserved integral. Commit topology, distribution, fields, and
diagnostics only after all ranks accept the complete candidate.

## Consequences

The sparse MPI tree can create, reshape, deepen, preserve, and collapse its
branching topology from owner-local tags without a complete numerical-tree
materialization. Compact topology and EB geometry remain replicated. A rank
mismatch, invalid geometry, EOS failure, transfer failure, count overflow, or
conservation failure leaves the accepted tree unchanged. Arbitrary-depth
checkpoint/restart and composite output remain later lifecycle boundaries.
