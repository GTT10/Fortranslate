# Decision 0192: Persist patch-tree regrid history

## Context

Patch-tree checkpoints retain the number of topology changes, but public
drivers discarded the number of scheduled tag evaluations and their cumulative
tagged-cell count. Restart diagnostics therefore described only the process
suffix and could not distinguish an unchanged evaluation from no evaluation.

## Decision

Count each tag/regrid invocation after its full transaction succeeds and sum
its global tagged-cell result. Store both values in base schema 5 and
fingerprinted schema 8. Require `regrids <= evaluations <= steps + 1`,
nonnegative tagged cells, and an all-zero history when no evaluation exists.

Sparse writes require all ranks to agree on presence and value before field
gather. The selected I/O root serializes the values, and restart broadcasts
them with the existing clock metadata before owner reconstruction.

## Consequences

Serial and changed-rank restarts report the same adaptation history as an
uninterrupted run. Evaluations that keep the current topology remain visible;
failed transactions do not alter history. Older strict-schema checkpoints are
intentionally incompatible.
