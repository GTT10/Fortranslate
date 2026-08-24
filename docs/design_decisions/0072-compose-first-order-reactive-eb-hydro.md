# Decision 0072: compose a first-order reactive EB hydro update

## Context

The geometry, slip-wall flux, open-area divergence, and weighted StateRedist
kernels were independently verified, but callers still had to manufacture all
Cartesian face fluxes. PeleC's EB hydro path constructs face fluxes, forms the
pre-redistribution divergence, and passes `U_old + dt*dUdt` through the selected
redistribution method. A smaller complete path is needed before porting the
high-order cut-cell predictor and face-centroid interpolation.

## Decision

Add one two-dimensional first-order EB hydro module. For every
positive-aperture Cartesian face, call the selected general-EOS reactive Riemann solver
with the adjacent cell-centered conserved states. Require both cells on an
interior open face to be active. Leave every closed face exactly zero. At a
nonperiodic domain face, use the adjacent fluid cell on both sides to implement
a piecewise-constant zero-gradient boundary.

Feed those face fluxes to the existing open-area and integrated slip-wall
divergence. Apply the resulting right-hand side through the default weighted
StateRedist provisional-state update, and commit state and temperature only
after every active cell passes EOS recovery. Preserve the face-flux builder as
a separate public verification seam.

## Consequences

PeleF now has an end-to-end first-order reactive Euler update for a stationary
embedded boundary. Uniform stationary pressure is testable through the entire
pipeline, not only through hand-supplied face fluxes. This milestone does not
claim PeleC's high-order EB reconstruction, transverse prediction,
face-centroid flux interpolation, periodic or inflow boundary treatment,
moving walls, chemistry/transport splitting, AMR coupling, or MPI ownership.

The operator sequence follows PeleC's current
[`pc_umdrv_eb`](https://github.com/AMReX-Combustion/PeleC/blob/development/Source/Hydro.cpp)
path, which constructs EB face fluxes and a pre-redistribution divergence
before calling AMReX-Hydro redistribution.
