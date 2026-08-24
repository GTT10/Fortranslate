# Decision 0053: send sparse owner migrations point to point

## Context

Same-hierarchy sparse migration retained only the new owner's payload but used
one broadcast per patch, sending complete fields through every rank. An owner
map change has exactly one authoritative source and one required destination.

## Decision

Pack a reactive patch's state, temperature, narrow ghosts, and wide ghosts into
one contiguous real payload. When ownership changes, send that payload directly
from the old owner to the new owner with one ordered MPI message. When ownership
is unchanged, assign locally without communication; unrelated ranks allocate no
payload buffer.

Keep the deterministic level/patch order and place a collective acceptance
boundary after each transfer. Count a committed transfer on the sending rank so
the communicator sum must equal the number of changed-owner patches.

## Consequences

Same-hierarchy rebalancing traffic is no longer replicated across the
communicator, and every reactive patch array is reconstructed exactly on its new
owner. One-copy sparse storage and existing gather parity are preserved.

The ordered blocking schedule favors simple failure boundaries over overlap.
Physics-stage streaming, halo exchange, and topology-changing regrid still use
collective communication and remain separate point-to-point work.
