# Decision 0114: target sparse EB average-down transfers

## Context

Direct sparse reactive average-down computes one compact coarse-footprint
restriction on each child owner, but broadcasts that buffer to every rank.
Only root tile owners whose row bands intersect the child footprint can use the
payload. The broadcast therefore scales traffic and temporary allocation with
the communicator rather than with the actual coarse/fine relationship.

## Decision

For each child, derive a rank set from the root tiles whose y intervals overlap
the child's coarse footprint. Deduplicate owners, retain a local copy when the
child owner is in the set, and send the restriction once to every other owner
with a blocking point-to-point message. Unrelated ranks do not allocate a
restriction buffer.

Keep the existing collective acceptance after each child so EOS recovery on
all intersecting root owners remains one transaction. Count actual remote sends
on the child owner, but publish the count only when the complete average-down
succeeds.

## Consequences

Communication is proportional to distinct coarse owners touched by a child,
including zero messages when all work is local. Conserved restriction, covered
cell behavior, serial parity, and rollback semantics are unchanged.

Root physics synchronization, coarse-time exterior transfer, flux and reflux
traffic, distributed root StateRedist, public time-loop control, regridding,
checkpointing, and output remain future work.
