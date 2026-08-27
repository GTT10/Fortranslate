# Decision 0181: Grow rank-deficient cut-parent stencils

## Context

The connected 3-by-3 cut-parent least-squares fit is compact and exactly
reproduces a supported affine field. Near a narrow fluid passage, however, all
immediate connected centroids can lie along one direction even when the fluid
path turns just outside that box. Accepting the rank-one fit there discards a
smooth gradient component that a slightly larger local neighborhood resolves.

## Decision

Build the 3-by-3 normal system first. If it is rank deficient, rebuild the
system and component envelope over a parent-centered 5-by-5 box. Include a
candidate only when an open-face path from the parent reaches it without
leaving the box. Solve the grown full-rank system normally; use the existing
minimum-norm rank-one or zero-gradient fallback only if the grown system is
still deficient.

## Consequences

Compact full-rank stencils are unchanged. Narrow or turning fluid paths can
recover both components of a smooth affine gradient without crossing covered
geometry or a disconnected fluid component. The larger stencil is used only
for deficient cut parents, so its extra bounded search does not affect regular
parents or already-resolved cut parents. Existing fine-child limiting,
fine-volume-weighted conservation, EOS recovery, and PCM retry remain the
acceptance boundary. This is not quadratic reconstruction or exact AMReX
interpolation.
