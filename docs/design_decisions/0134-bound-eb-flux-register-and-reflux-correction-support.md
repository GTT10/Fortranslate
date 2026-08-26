# Decision 0134: bound EB flux-register and reflux correction support

## Context

Each EB flux register previously allocated a correction over the complete
coarse level even though coarse/fine flux mismatch is accumulated only in the
one-cell exterior ring of its fine patch. Sparse transport also moved a complete
corrected root to and from every remote child for each Euler stage. Those
round trips preserve deterministic cumulative reflux, but most payload values
cannot be changed by the selected child.

For a regular interface, a mismatch changes only its exterior coarse cell. For
a cut exterior cell, the established conservative re-redistribution can also
change a connected cardinal or diagonal neighbor, including a fine-covered
recipient. No correction can propagate farther during one child reflux.

## Decision

Store the flux-register correction on the patch rectangle expanded by one
coarse cell and preserve absolute coarse indices as the allocatable array's
lower and upper bounds. Validate those bounds with the register and iterate
only that support during reflux.

For sparse transport, transfer the current corrected root between the root
physics owner and a remote child only on the patch rectangle expanded by two
coarse cells. A child initializes its compatibility workspace from the
uncorrected root candidate, overwrites that rectangle with the latest
cumulative values, performs the existing reactive reflux transaction, and
returns the same rectangle. The root owner merges returned rectangles in the
established child order.

## Consequences

Flux-register storage scales with patch support instead of root-level area,
and repeated child correction messages carry only the region that reflux can
change. Message counts, correction order, EOS recovery, child ownership, and
atomic rollback remain unchanged.

Each distinct child owner still receives a complete root start/end/temperature
and x/y-flux bundle once per Euler stage. Replacing that compatibility input
with compact exterior and interface-flux data is subsequent work.
