# Decision 0075: stage EB PLM fluxes at face centroids

## Context

The runnable EB path has valid open-area fractions and conservative
StateRedist, but its Cartesian Riemann states are piecewise constant and its
fluxes are evaluated only at Cartesian face centers. PeleC's EB hydro driver
obtains AMReX face centroids, computes Godunov fluxes into temporary
face-centered arrays, and supplies those geometry fields when forming the EB
divergence. PeleC's two-dimensional face stencil linearly combines the local
face-center flux with one tangential neighbor according to the signed
normalized face-centroid offset.

## Decision

Store one normalized tangential centroid offset with every x- and y-face.
Derive it from the same linearly clipped positive segment as the open fraction,
so shared geometry remains deterministic and full or closed faces have zero
offset.

Add a separate EB reconstruction module. Keep PCM as the default public
behavior and add selectable frozen-composition characteristic PLM. Reuse the
qualified regular reactive characteristic projection and normal trace, but
form a nonzero directional slope only when both normal stencil neighbors are
active. A covered neighbor or outer boundary therefore produces a local
zero-slope fallback without reading a covered primitive state.

First compute unweighted Riemann fluxes at Cartesian face centers. Then
interpolate each partial-face flux toward the signed tangential neighbor using
the normalized centroid offset. If that neighbor is outside the domain or its
face is closed, retain the local face-center flux. Apply the open fraction only
later in the existing conservative divergence.

## Consequences

The EB application can run either PCM or time-traced characteristic PLM and
uses open-face centroid geometry in both modes. Linear moving-contact fluxes
are exact on active regular stencils, and geometry-adjacent stencils never
require a covered cell. Existing callers that omit reconstruction arguments
retain PCM behavior.

This is not full PeleC EB Godunov parity. It excludes unsplit transverse
prediction, PPM, higher-order StateRedist slopes, periodic EB ghost
neighborhoods, and multilevel flux registers.

Primary references:

- [PeleC `pc_umdrv_eb`](https://github.com/AMReX-Combustion/PeleC/blob/development/Source/Hydro.cpp)
- [PeleC face-centroid interpolation stencil](https://github.com/AMReX-Combustion/PeleC/blob/development/Source/EB.cpp)
- [AMReX embedded-boundary documentation](https://amrex-codes.github.io/amrex/docs_html/EB.html)
