# Decision 0178: Prolong cut EB parents conservatively

## Context

Limited-linear EB AMR prolongation is active across every public lifecycle,
but its first qualified kernel uses PCM for every cut parent. That avoids
topology-sensitive interpolation, yet introduces a first-order band beside the
embedded boundary and discards smooth coarse variation during regridding.

## Decision

For an EB-cut parent, compute directional conserved-state differences between
active coarse fluid-volume centroids. Apply MC limiting with neighbors on both
sides and use the available one-sided derivative when covered geometry removes
one side. Bound every reconstructed conserved component by the active 3-by-3
coarse-neighbor envelope.

Evaluate offsets at active fine fluid centroids. Subtract the
fine-volume-fraction-weighted mean offset within each parent before applying
the slopes. This makes the reconstructed offset contribution sum to zero under
the same weights used by EB average-down.

Retain PCM for covered parents and regular parents whose child topology is not
fully regular. Recover every active child temperature through the EOS and
retain the existing parent-local PCM retry if any linear child is inadmissible.

## Consequences

Configured linear prolongation no longer inserts an automatic PCM band at cut
parents, and average-down still recovers their conserved state to roundoff.
The method is monotone component by component over the local active envelope.
It is not an exact reproduction of AMReX interpolation, a multidimensional
least-squares gradient, or a higher-order cut-cell predictor.
