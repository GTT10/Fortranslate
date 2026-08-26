# Decision 0131: blend sparse SSPRK2 root on tile owners

## Context

Each sparse EB transport Euler stage gathers root tiles to the root physics
owner because diffusive face construction and StateRedist depend on neighboring
cells. After the second stage, the SSPRK2 root blend separately gathers the
interval-start and second-Euler candidates to that owner, averages them,
recovers temperature, and scatters the result back to tile owners.

The final conserved-state average and EOS temperature recovery are cell-local.
They do not use a face stencil, a redistribution neighborhood, or any other
root tile's numerical values. The extra two gathers and one scatter therefore
cross the sparse ownership boundary without a numerical dependency.

## Decision

Retain the interval-start and second-Euler candidates in their existing sparse
tile allocations. On every root tile owner, average the two local conserved
states, extract the tile's exact EB geometry row band, and recover temperature
through the existing transport EOS helper. Accept this phase collectively
before blending children and performing the established sparse average-down.

Count only the root transfers performed by the two Euler stages. Do not gather
either blend candidate or scatter the blended root. Preserve the complete
outer transaction so any local EOS failure rejects without publishing state,
Euler counts, limiter values, or transfers.

## Consequences

Each remote root tile now sends four payloads per SSPRK2 transport call instead
of seven: one gather and one scatter in each Euler stage. The final blend uses
no root numerical-field communication or complete-root allocation. Root and
child storage remains globally single-copy between stages.

The Euler stages remain level-wide on the root physics owner. Removing those
gathers requires finite-halo diffusive flux construction, StateRedist support,
and deterministic shared-face result ownership, which is separate work.
