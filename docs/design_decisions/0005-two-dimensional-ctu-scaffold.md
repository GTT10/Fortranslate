# Decision 0005: add a separate periodic 2D CTU scaffold

## Status

Accepted for the uniform-grid single-species milestone.

## Context

The verified one-dimensional path already contains characteristic normal prediction and selectable Riemann solvers. Extending the existing 1D array and driver directly to multiple dimensions would mix indexing, boundary, reconstruction, and update changes in one step and make regressions difficult to localize.

PeleC's multidimensional Godunov path also applies transverse flux corrections after normal prediction. A two-dimensional analytical problem is needed before AMR, chemistry, or spray can be layered onto multidimensional flow.

## Decision

Introduce a separate `pelef2d` executable and `state(variable,i,j)` storage. The 2D operator will:

- initially support periodic uniform Cartesian meshes only;
- reuse the existing primitive conversion, limiter, characteristic tracing, and Riemann solvers;
- obtain y-direction fluxes through explicit x/y momentum rotation;
- build time-centered normal predictor states in both directions;
- evaluate provisional x and y fluxes;
- correct each face state with half of the transverse provisional flux divergence from its originating cell;
- evaluate final face fluxes and apply one unsplit conservative update;
- scale only the transverse correction if it would violate density or pressure positivity;
- retain an option to disable transverse correction for differential verification.

## Consequences

Benefits:

- the existing 1D implementation remains unchanged;
- direction rotation and transverse correction have focused unit tests;
- dimensional reduction can be checked directly against the 1D solver;
- periodic conservation follows from shared face fluxes;
- the isentropic vortex provides an analytical multidimensional convergence gate;
- the effect of transverse correction is measurable rather than nominal.

Limitations:

- current 2D boundaries are periodic only;
- componentwise MC/minmod slopes feed the normal characteristic tracing;
- there are no source terms, PPM/WENO, embedded boundaries, AMR, or 3D double-transverse corrections;
- the state remains single-species and constant-`gamma`.

The next milestone will add multispecies conserved-state infrastructure while preserving the five-variable Euler layout as a specialization.
