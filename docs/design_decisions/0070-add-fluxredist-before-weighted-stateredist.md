# Decision 0070: add FluxRedist before weighted StateRedist

## Context

An explicit conservative divergence divides by the cut-cell volume fraction,
so an arbitrarily small cell can impose an arbitrarily small stable timestep.
PeleC currently selects AMReX-Hydro `StateRedist` by default. AMReX also defines
a first-order conservative FluxRedist construction whose conservation and
small-cell scaling can be qualified independently before higher-order weighted
state neighborhoods are added.

## Decision

For every cut cell, use itself and active Cartesian face neighbors connected by
a positive aperture. Form a volume-fraction-weighted neighborhood right-hand
side `Rnc`, replace the cut-cell value by
`kappa*Rc + (1-kappa)*Rnc`, and distribute
`kappa*(1-kappa)*(Rc-Rnc)` equally per unit neighbor fluid volume. Reject a cut
cell with no in-domain connected neighbor.

Apply the redistributed right-hand side to the full reactive conserved state
in temporary storage. Recover the general-EOS temperature in every active cell
and commit state and temperature together only if all cells are physical.

## Consequences

The update is conservative, preserves uniform right-hand sides, and removes the
unscaled inverse-volume contribution from the cut cell. It is a qualified
first-order FluxRedist path, not a claim of PeleC's default weighted
StateRedist, higher-order slopes, periodic/physical-boundary neighborhoods,
AMR redistribution, or MPI distribution.

The formulas and default-method distinction follow the
[AMReX embedded-boundary documentation](https://amrex-codes.github.io/amrex/docs_html/EB.html)
and [PeleC runtime parameters](https://github.com/AMReX-Combustion/PeleC/blob/development/Source/Params/_cpp_parameters).
