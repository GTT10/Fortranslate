# Decision 0191: Persist patch-tree operator counters

## Context

The full-physics kernels expose committed per-level chemistry, transport, and
hydro patch advances, but public drivers previously discarded them. A restart
therefore had no run-wide operator-work ledger, and sizing history only to the
current tree would lose deeper-level counts after dynamic coarsening.

## Decision

Maintain three cumulative vectors sized to the configured maximum AMR depth.
Add the common capacity and all three vectors to base checkpoint schema 4 and
fingerprinted schema 7. Populated levels contribute committed deltas; dormant
slots retain their prior values.

Sparse kernels return owner-local patch counts. The public MPI application
uses communicator sums before accumulation. Sparse checkpoint writes require
all ranks to agree on presence, capacity, and values; the selected I/O root
writes them and broadcasts them during restart before owner redistribution.

## Consequences

Serial and changed-rank continuations report the same complete-run per-level
operator counts as uninterrupted execution, including across topology shrink
and regrowth. A configured restart depth smaller than the stored counter
capacity is incompatible. Older strict-schema checkpoints are intentionally
rejected, and failed reads leave optional counter arrays unallocated.
