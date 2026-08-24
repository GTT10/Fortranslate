# Decision 0076: add second-order weighted StateRedist

## Context

The serial reactive EB path already matches the PeleC-default weighted
StateRedist configuration with `eb_srd_max_order = 0`. That construction is
conservative and stabilizes small cells, but every merge neighborhood scatters
one constant `Qhat`; it therefore diffuses an affine state even when geometry
and the provisional solution are locally smooth.

AMReX locates each `Qhat` at a weighted neighborhood centroid, fits slopes from
neighboring `Qhat` values, limits predictions at those centroids, and evaluates
the reconstruction at every cell receiving part of the neighborhood. Its
`max_order=2` setting denotes a linear reconstruction and second-order method,
not a quadratic polynomial.

## Decision

Store normalized x/y offsets of each fluid-volume centroid alongside volume
fraction. Compute those moments from the same two clipped positive polygons
used for cell volume, leaving regular and covered offsets exactly zero.

Retain `state_redist_max_order=0` as the API default and support the AMReX
second-order value `2`. Form `cent_hat` with exactly the self, neighbor,
volume-fraction, and overlapping-neighborhood weights used by `Qhat`. Fit a
linear slope on connected active 3-by-3 neighborhoods, growing to an active
5-by-5 stencil when the normal matrix is rank deficient. Apply the AMReX
pairwise centroid limiter. Then apply one additional common slope scale so all
self and merge-recipient evaluations remain inside the active component range.

Scatter reconstructed values through the existing partition. Since `cent_hat`
is the weighted centroid of those same recipients, the weighted sum of every
linear displacement is zero. The higher-order correction therefore preserves
the established componentwise conservation contract. The runnable EB cases
select order 2 explicitly.

## Consequences

Affine fields at fluid centroids survive overlapping StateRedist neighborhoods,
uniform-state behavior and order-0 callers remain unchanged, and discontinuous
recipient values obey a maximum principle. The extra recipient limiter is a
serial safety layer beyond the upstream centroid limiter; it uses a global
active-cell component range and will need a reduction or localized equivalent
when this path becomes distributed.

This does not implement AMReX order 4, boundary or periodic ghost
neighborhoods, multilevel redistribution, or EB coarse/fine synchronization.

Primary references:

- [AMReX `StateRedistribute`](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EB_StateRedistribute.cpp)
- [AMReX `MakeStateRedistUtils`](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EB_StateRedistUtils.cpp)
- [AMReX StateRedist centroid limiter](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EB_StateRedistSlopeLimiter_K.H)
- [AMReX two-dimensional EB slopes](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EB_Slopes_2D_K.H)
