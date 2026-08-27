# Decision 0154: Route arbitrary-depth EB chemistry on owners

## Context

Sparse owners can select the global timestep, but chemistry still requires the
complete serial numerical tree. Reaction itself is node-local; the distributed
boundary is the deepest-first average-down that must publish each reacted child
into its parent before the next shallower relation is synchronized.

## Decision

Advance every node on its owner inside one private sparse candidate. Recover
temperature locally and collectively accept after each node. Then traverse the
topology relations deepest-first and children in deterministic order. Apply
average-down locally when parent and child share an owner. Otherwise send the
complete child conserved state directly to the parent owner and apply the same
serial EB restriction kernel there.

Collectively accept after every child restriction and commit only after the
complete sparse candidate validates. Count owner-local chemistry calls by
level and count each distinct-owner child send exactly once.

## Consequences

Chemistry and hierarchy restriction no longer require materializing numerical
fields on nonowners. The result remains exactly equal to the serial transaction,
including parent temperature recovery and ordered multiple-child updates.

Recursive hydrodynamics and transport still require time-dependent parent
context, interface fluxes, reflux, and conservation closure. They remain
separate owner-local routing milestones.
