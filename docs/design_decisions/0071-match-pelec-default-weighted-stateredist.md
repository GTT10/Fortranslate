# Decision 0071: match PeleC's default weighted StateRedist

## Context

PeleC selects AMReX-Hydro `StateRedist` by default and sets
`eb_srd_max_order = 0`. Unlike FluxRedist, this method redistributes the
provisional conserved state. Its conservation depends on accounting for cells
that belong to more than one merge neighborhood; independently averaging each
small cell would double-count those shared volumes.

## Decision

Implement the two-dimensional, zeroth-order algorithm used by AMReX. A cell
with volume fraction below the default target `0.5` chooses its first neighbor
from the dominant component of the face-aperture-difference normal. An equal
normal-component magnitude or insufficient accumulated volume adds the
orthogonal neighbor and their connecting diagonal, for at most three
neighbors.

Count how many neighborhoods contain every active cell (`nrs`). For merging
cell `i`, set its neighbor weight to
`beta_i = (0.5-kappa_i)/sum_(j in M_i) kappa_j`. Reduce each recipient's self
weight by `beta_i/nrs_j`, form the resulting volume-weighted neighborhood state
`Qhat_i`, distribute each `Qhat_i` back through the same partition, and divide
every recipient by `nrs_j`.

Apply this operator to `U_old + dt*R`. Recover the general-EOS primitive state
and temperature for all active cells and commit only if every recovery
succeeds. Keep the earlier FluxRedist API as a distinct supported method.

## Consequences

The partition is conservative for every component even when merge
neighborhoods overlap, and it preserves a uniform state. The implementation
matches PeleC's current default zeroth-order choice on one nonperiodic
Cartesian level. It does not claim higher-order StateRedist reconstruction and
limiting, periodic or physical ghost-cell neighborhoods, multilevel
redistribution, or distributed ownership.

The construction follows the current AMReX
[`MakeITracker`](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EB_StateRedistItracker.cpp),
[`MakeStateRedistUtils`](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EB_StateRedistUtils.cpp),
and
[`StateRedistribute`](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EB_StateRedistribute.cpp)
implementations. The selected default and order follow the
[PeleC runtime parameters](https://github.com/AMReX-Combustion/PeleC/blob/development/Source/Params/_cpp_parameters).
