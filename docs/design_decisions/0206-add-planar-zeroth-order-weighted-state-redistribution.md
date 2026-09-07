# Decision 0206: add planar zeroth-order weighted state redistribution

## Context

The qualified 3D EB kernel already exposed both the raw cut-cell Euler update
and first-order FluxRedist. FluxRedist redistributes the conservative residual,
whereas PeleC's principal small-cell route redistributes a provisional
conserved state through overlapping weighted neighborhoods. Conflating those
algorithms would hide their different diffusion and admissibility behavior.

The only qualified 3D geometry is an exact x-, y-, or z-normal plane. Its
small cut cell has one full regular cell in the aperture-difference normal
direction. Generalizing the neighborhood builder before general 3D geometry
exists would create an untestable interface.

## Decision

Add a distinct zeroth-order weighted StateRedist operator. For every active
cell, let `nrs` count its own neighborhood and every small-cell neighborhood
that names it. A cut cell with volume fraction `kappa < target` selects the
single regular receiver along the nonzero component of its aperture-
difference normal and assigns

```text
alpha_neighbor = (target - kappa) / kappa_receiver.
```

The receiver's self partition is reduced by
`alpha_neighbor / nrs_receiver`. The same partition first forms each
volume-weighted neighborhood state `Qhat`, then scatters `Qhat` back through
the self and receiver weights and divides every recipient by its `nrs`.
Using the identical partition in both directions preserves
`sum(kappa U)` component by component. Covered output is zero in the
standalone operator, while the reactive update keeps covered caller state
unchanged.

Reject a non-axis aperture-difference normal, a missing or nonregular receiver,
an invalid target outside `(0,1]`, a negative self partition, invalid geometry
or extents, and nonfinite input. Do not silently broaden this operator to an
oblique or curved geometry.

Provide a low-level transactional update that forms `U* = U + dt R`, applies
weighted StateRedist, and recovers every active NASA7 state. Add a matching
end-to-end Euler route. Refactor the hydro module so raw, FluxRedist, and
StateRedist all consume one common PCM face-flux and EB-divergence builder.

## Consequences

At `kappa=0.05`, `target=0.5`, and a full receiver, the receiver weight is
`0.45`. If the provisional cut state is `-U` and the receiver is `U`, the
zeroth-order neighborhood gives

```text
Qhat_cut = 0.6363636363636364 U
U_receiver,new = 0.9181818181818182 U.
```

The x/y/z analytical results, uniform-state identity, component conservation,
covered-cell contract, invalid-input transactions, and EOS rollback are
independently testable. The full-grid CFL control that makes the raw cut
density negative remains physical through both FluxRedist and StateRedist,
but their resulting states intentionally differ.

This decision does not claim second- or higher-order StateRedist, general
multi-neighbor geometry, direct AMReX/PeleC field parity, characteristic EB
reconstruction, a public EB timestep/application, chemistry or transport at
EB cells, AMR/reflux/regridding, restart, MPI ownership, or scalable I/O.
