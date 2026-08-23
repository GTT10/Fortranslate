# Decision 0028: retain aligned fine cells across multilevel regrids

## Status

Accepted for PeleF `0.33.0`.

## Context

The first dynamic multilevel rebuild preserved composite integrals by averaging
the old hierarchy to the root and prolonging a new hierarchy. It discarded old
fine-scale information whenever any patch bound changed, even where an old and
new patch represented the same physical cells at the same resolution.

## Decision

After conservative hierarchy construction, intersect each common old/new fine
level in physical coordinates. Directly copy conserved state and temperature
only when the level spacings agree and both overlap edges align with cell
boundaries. Count transferred cells for diagnostics. After every level copy,
average down the rebuilt hierarchy from deepest level to root.

If spacing differs, retain the already constructed conservative prolongation.
Do not interpolate between unequal fine resolutions in this slice.

## Consequences

- moving and resized patches retain all aligned fine-scale overlap exactly;
- copied deepest states become authoritative in every covered parent;
- composite conservation from synchronized-root reconstruction is retained;
- changed refinement ratios use a safe conservative fallback;
- multiple-patch intersections and cross-resolution remapping remain future
  work.
