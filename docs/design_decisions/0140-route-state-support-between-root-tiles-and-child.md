# Decision 0140: route state support between root tiles and child

## Context

Root tile owners already compute the coarse transport start/end fields and
retain the interface fluxes needed by each child. The root physics owner still
extracted every child exterior context, supplied cumulative reflux support,
merged the returned correction, and scattered the final root rows. Those
steps kept state traffic centralized after the tile computation was complete.

## Decision

Retain three state/temperature views on every transport tile owner: the stage
start, the uncorrected stage end used for coarse-time exterior interpolation,
and the current corrected end used for cumulative reflux support. For each
child, route only intersecting patch-plus-two row fragments to its owner. Build
the exterior context there through the globally indexed support API.

After child-local reflux, route corrected support fragments directly back to
the intersecting tile owners in deterministic child order. Commit those tile
states without a root-owner scatter. Keep the complete temporary root bundle
for compatibility validation and cut-boundary flux closure in this milestone.

## Consequences

Child state context and reflux correction no longer pass through the root
physics owner. Overlapping child corrections remain ordered because each child
finishes its tile updates before the next child assembles support. The
uncorrected end view preserves the established exterior interpolation even
when an earlier child has changed cumulative corrected support. Message counts
are derived from remote tile/child intersections for both SSPRK2 Euler stages.
