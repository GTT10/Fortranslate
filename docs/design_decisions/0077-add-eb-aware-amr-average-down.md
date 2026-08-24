# Decision 0077: add EB-aware AMR average-down

## Context

The regular one-dimensional AMR hierarchy already restricts fine conserved
averages to a covered parent region. Embedded-boundary cells change the measure:
a fine cell contributes only its fluid volume, so an unweighted child average
does not preserve a coarse/fine composite integral through a cut interface.

AMReX's EB average-down kernel weights each fine value by its volume fraction,
divides by the total fine volume fraction, and uses the first child if that
total is zero. Correct transfer also assumes that the coarse and fine geometry
describe the same physical fluid measure.

## Decision

Introduce `amr_eb_hierarchy_2d_mod` for one static, aligned, rectangular fine
patch. Store its inclusive coarse bounds and integer refinement ratio. Accept a
patch only when the level dimensions, physical bounds, spacing, and restricted
parent/child volume fractions agree.

For each parent in the patch, restrict every state component with fine
volume-fraction weights. Retain the AMReX first-child fallback for a completely
covered child block. Leave the rest of the coarse array unchanged. Define a
matching composite integral that excludes the geometrically covered coarse
rectangle and includes the fine patch once.

Wrap generic restriction for reactive states. Recover temperature from the
restricted conserved state on every active parent and commit state plus
temperature only after all conversions succeed. Covered parents retain their
input reactive values rather than accepting arbitrary solid-cell fallback data.

## Consequences

Constant states, affine states sampled at exact fluid centroids, and the
coarse/fine composite integral survive restriction to roundoff. Misaligned or
measure-inconsistent hierarchies are rejected before transfer, and EOS failure
cannot leave a partially restricted reactive hierarchy.

This module does not provide fine initialization, ghost fill, time advancement,
subcycling, flux-register reflux, multiple patches, dynamic regridding,
multilevel redistribution, or MPI ownership.

Primary references:

- [AMReX EB MultiFab average-down dispatch](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EBMultiFabUtil.cpp)
- [AMReX two-dimensional EB average-down kernel](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EBMultiFabUtil_2D_C.H)
