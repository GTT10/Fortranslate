# Decision 0205: stabilize 3D EB with conservative flux redistribution

## Context

The first 3D EB Euler kernel divided the Cartesian and wall-flux balance by
the cut-cell fluid volume. That is the correct conservative residual, but its
explicit stability interval shrinks with the smallest positive volume
fraction. The axis-plane metric contract was already independently qualified,
so the next increment could isolate stabilization from general geometry,
higher-order reconstruction, and AMR synchronization.

The mature 2D path contains both FluxRedist and weighted StateRedist. Porting
the complete higher-order StateRedist neighborhood and fallback machinery at
the same time as the first 3D stabilization would make failures ambiguous.

## Decision

Add a first-order FluxRedist operator over the six face-connected Cartesian
neighbors. A cut cell forms a volume-fraction-weighted neighborhood average
using only active cells reached through positive-area faces. Its published
right-hand side blends the conservative cut residual with that average. The
operator then returns the removed excess to the same neighbors in proportion
to their total fluid volume.

Keep regular-cell residuals unchanged before receiving excess and force
covered-cell residuals to zero. Reject a cut cell with no active open-face
neighbor. Validate the complete geometry, array extents, and finiteness before
publishing the redistributed field.

Provide two composition boundaries: a low-level general-EOS update from a
caller-supplied conservative residual, and an end-to-end PCM EB Euler update
that builds face fluxes and wall divergence before redistribution. Both
recover active-cell temperature and leave caller state and temperature
unchanged on any failure. Retain the raw Euler routine as a diagnostic
reference.

## Consequences

The volume-fraction-weighted integral of every redistributed component equals
the original integral. A uniform active residual remains uniform. All x/y/z
normal plane orientations and the four transverse neighbors of an x-normal
cut sheet use the same six-face rule.

For the qualified `kappa=0.05` case, a synthetic residual is reduced to
`9.5238%` of its raw cut-cell magnitude. In the end-to-end hydro gate, a
CFL-0.5 full-grid timestep makes the raw cut density negative, while the
redistributed update remains physical and conservative.

This is deliberately FluxRedist, not a claim of current PeleC/AMReX weighted
StateRedist parity. General geometry, higher-order StateRedist, an EB-aware
timestep API, characteristic cut-face reconstruction, chemistry/transport,
AMR/reflux/regridding, restart, MPI, and a public 3D EB application remain
separate work.
