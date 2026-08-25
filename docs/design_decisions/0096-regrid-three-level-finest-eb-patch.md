# Decision 0096: regrid the three-level finest EB patch

## Context

The public three-level lifecycle previously required both nested rectangles to
remain fixed. The qualified two-level regrid transaction already conserves the
parent/child composite state while moving a child, but applying its unrestricted
planner to the finest level could violate the two-middle-cell stencil margin
required by recursive EB advancement.

## Decision

When `three_level_enabled` and `dynamic_regridding` are both selected, keep the
configured root-to-middle patch fixed and plan only the middle-to-finest patch.
Evaluate temperature-gradient tags on active middle cells, remove the outer
two-cell band from the planning domain, and apply the existing buffer and
minimum-size controls within that interior.

Reuse the transactional two-level topology replacement on the middle/finest
pair: average down the old finest state, prolong the replacement from the
synchronized middle candidate, copy geometrically consistent overlap cells,
recover active temperatures through the EOS, and publish both levels only
after complete validation. Apply the operation at initialization when
requested and after accepted root steps at `regrid_interval` cadence.

## Consequences

The public three-level application can move and resize the finest rectangle
without changing recursive subcycling, reflux, chemistry, CFL reduction, or
output. Tags outside the qualified interior are intentionally ignored, and an
empty interior plan retains the current finest patch.

This mode does not remove the finest level, change the root-to-middle patch,
create sibling finest patches, or checkpoint a dynamic topology. Those
capabilities, molecular transport, arbitrary depth, and MPI ownership remain
separate work.
