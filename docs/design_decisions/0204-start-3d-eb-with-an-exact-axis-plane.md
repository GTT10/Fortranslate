# Decision 0204: start 3D EB with an exact axis plane

## Context

The qualified 3D path had regular-grid and static-AMR hydrodynamics, but no
embedded-boundary data model. Porting general triangulated geometry,
small-cell stabilization, AMR cut-face registers, regridding, and distributed
storage in one increment would mix geometric errors with EOS, flux-balance,
and hierarchy errors. The existing 2D EB implementation supplies useful
semantics, but its nodal polygon construction and redistribution algorithms
cannot establish that a new 3D metric convention is correct by themselves.

## Decision

Introduce a deliberately narrow analytic geometry: one static x-, y-, or
z-normal plane, with fluid on its coordinate-positive side. Store Cartesian
cell volume fractions, normalized fluid centroids, all three face apertures
and tangential face centroids, EB area and physical centroid, and the
solid-to-fluid unit normal and area integral. Require the type validator to
check allocation bounds, finiteness, ranges, cell classification, interface
metrics, and the cellwise identity between Cartesian aperture divergence and
the EB normal integral.

Treat a plane at or beyond the domain boundary as a fully regular or covered
domain. Reject an interior plane exactly on a grid face: this first ownership
model requires a positive-volume cut cell to store the EB metric.

Recover wall pressure through the existing NASA7 mixture EOS. Apply the
impermeable slip-wall flux using the outward normal of the fluid control
volume, combine it with aperture-weighted x/y/z Cartesian fluxes, and divide
by the cut-cell fluid volume. Construct PCM Riemann fluxes on open faces with
zero-gradient physical-domain exterior states. Publish a forward-Euler state
and recovered-temperature candidate only after every face, wall, divergence,
and EOS operation succeeds.

## Consequences

The exact metric identity makes uniform pressure and uniform flow tangent to
the plane invariant for all three coordinate orientations. Internal face
fluxes cancel in fluid-volume integrals, and failure paths preserve the input
state. The increment establishes the sign, area, indexing, and transaction
conventions needed by later 3D EB work without claiming general geometry.

The update has no small-cell redistribution, so its stable timestep remains
limited by the smallest positive volume fraction. It is PCM, serial,
single-level, inviscid, nonreacting during the step, and has no restart or
field-output application. Oblique and curved surfaces, grid-face-aligned
interior ownership, high-order reconstruction, chemistry/transport,
AMR/reflux/regridding, and MPI remain separate increments.
