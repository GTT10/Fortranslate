# Decision 0045: compose owner physics in one MPI transaction

## Context

The MPI patch-tree bridge has qualified owner-only chemistry, recursive hydro,
and recursive molecular transport entry points. Calling them independently
does not guarantee that a later failure restores changes committed by earlier
operators in the same physical timestep. The serial patch-tree API already
defines one transactional `R-T-H-T-R` interval.

## Decision

Add an outer MPI reactive step with the serial operator order
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. Synchronize every owner patch and make
the root patch owner authoritative for time, step, level counters, transport
counters, and regrid statistics before copying the complete outer backup.

Retain collective acceptance and rollback inside each component operator. The
outer wrapper additionally accepts each completed operator collectively. If
any stage fails, restore the outer backup on every rank and report zero
chemistry, hydro, and transport calls. Reject a missing transport database
before the first operator runs.

## Consequences

One public MPI call now matches the serial patch-tree full-physics transaction,
including state, temperature, ghosts, bookkeeping, exact owner call counts,
and rollback after a later-stage failure. Component entry points remain useful
for verification and split orchestration.

The transaction still operates on owner-authoritative replicas and broadcasts
complete patches and flux fields. Sparse rank-local storage, ownership
migration after regrid, and scalable point-to-point communication remain
separate work.
