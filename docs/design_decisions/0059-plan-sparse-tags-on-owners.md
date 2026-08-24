# Decision 0059: plan sparse tags on owners

## Context

After explicit-plan regrid became field-sparse, the tag-driven API still
gathered the complete current tree on every rank so the serial planner could
tag, cluster, and prolong candidate levels. The resulting plan was compact, but
deriving it retained the last all-rank field replica in sparse AMR.

## Decision

Create a rank-local planning copy and average it down through the existing
direct child-to-parent path. At each candidate depth, only the owner of a parent
evaluates gradients and clusters tags. Agree on its tagged-cell count and child
bounds with integer reductions, then construct the replicated compact hierarchy
metadata deterministically.

Build candidate states on their owners with the direct parent-to-child
prolongation path. Repeat through the configured depth and pass the completed
plan to the replica-free explicit regrid transaction. Count owner tagging
evaluations and candidate transfers separately from final prolongation and
overlap messages.

Reset cached temperatures before recovering them from newly installed
conserved fields. A previous temperature belongs to a different state and must
not influence the iterative EOS inversion endpoint.

## Consequences

Both topology-change APIs retain globally single-copy field data. Ranks share
compact topology and ownership metadata plus collective acceptance decisions,
not complete state or temperature trees. Candidate construction may repeat
prolongation for each depth, matching the qualified serial planning order while
keeping memory rank-local.
