# Decision 0047: advance chemistry directly on sparse MPI AMR patches

## Context

The `0.54.0` sparse container removes persistent non-owner payloads, but all
qualified physics still transitions back to a complete replica. Chemistry is
the simplest operator with which to establish direct sparse execution because
the reaction solve is cell-local; only post-reaction AMR average-down and
ghost consistency require communication.

## Decision

Advance each patch only on the rank that stores it. After all reactions,
traverse relations deepest-to-root. Stream each child interior to its parent
owner, which performs conservative average-down and temperature recovery.
Then traverse root-to-deepest: fill root physical ghosts on its owner, stream
each parent to child owners for coarse/fine ghost fill, and stream adjacent
siblings one at a time for normal and PPM-wide ghost replacement.

Retain a rank-local sparse backup. Apply communicator-wide acceptance after
every patch reaction and synchronization boundary. Any rejection restores all
local payloads and reports zero accepted patch calls.

## Consequences

Chemistry no longer needs a complete reactive patch tree on any rank. The
same sparse path covers separated multilevel branches and adjacent PPM child
sets with exact serial parity and rollback.

The current synchronization is intentionally broadcast-based and may expose
one transient patch or one parent's child set on non-owners. Hydro, transport,
full-physics composition, point-to-point schedules, and topology-changing
distributed regrid remain separate work.
