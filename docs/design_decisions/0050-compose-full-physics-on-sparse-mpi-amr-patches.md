# Decision 0050: compose full physics on sparse MPI AMR patches

## Context

Chemistry, recursive hydro, and recursive molecular transport each have a
qualified direct sparse entry point. Calling the older replicated outer
transaction would still reconstruct a complete reactive tree on every rank
for a normal full-physics interval.

## Decision

Compose the direct sparse component transactions as
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. Collectively validate the interval,
sparse hierarchy, and optional transport database before mutation. Take one
outer sparse backup and collectively accept each component result, the final
ghost refresh, and final sparse validity check.

Accumulate local operator calls only after their stages succeed. If any later
stage is rejected, restore the outer backup on every rank and report zero
committed chemistry, hydro, and transport calls.

## Consequences

The complete full-physics interval now keeps persistent patch payloads
rank-local and matches the serial operator order, result, bookkeeping, and
transaction boundary. Missing transport data is rejected without mutation,
and a hydro failure after valid chemistry and transport prefixes restores the
entire starting sparse solution exactly.

Component routines retain nested backups and correctness-first collective
streaming. Reducing temporary copies, replacing broadcasts with point-to-point
schedules, and topology-changing distributed regrid remain separate work.
