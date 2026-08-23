# Decision 0034: mirror patch-tree ownership in the flux-register hierarchy

## Context

The arbitrary-depth patch tree can prolong and average fields, but reactive
time integration also requires the fine and coarse time-integrated fluxes at
every child boundary. A flat register array cannot identify which parent field
receives a correction when several parent patches exist at the same level.

## Decision

Store one register array for every parent-owned child set and one flux register
for each local child patch. Keep the register topology identical to the
relation and parent topology of the patch tree, including allocated zero-length
arrays for parents without children.

Synchronize from the deepest relation toward the root. For every parent, apply
all child reflux corrections and then average down the complete local child
set. Deep-copy both fields and registers before synchronization and restore both
if any validation, reflux, or restriction operation fails.

## Consequences

Every future recursive hydro or transport node has an unambiguous register in
which to accumulate its boundary fluxes. Synchronization propagates a deepest
correction through every ancestor without special cases for branching. The
reactive tree driver must still implement recursive subcycling, ghost fill, and
flux accumulation; dynamic ownership changes must rebuild both fields and
register topology together.
