# Decision 0115: target sparse EB hydro root traffic

## Context

Direct sparse hydro keeps fine children owner-local but materializes every root
tile on every rank. It then broadcasts root state, temperature, fluxes, and
each cumulative reflux correction across the communicator. Only the root
physics owner and actual child owners need complete root fields; root tile
owners need only their final row bands.

## Decision

Pack and send each nonlocal root tile once to the root physics owner. Advance
the existing level-wide root kernel there, then pack the root start, updated
state, temperature, and fluxes into one bundle sent once to each distinct
remote child owner.

Keep reflux ordering unchanged. Before a remote child advances, send it the
current corrected root; after reflux, return the corrected root to the physics
owner. When every child is accepted, scatter only the final row band to each
remote root tile owner and invoke targeted sparse average-down.

Count each packed MPI message as one payload transfer and publish the count
only with the final sparse commit. Retain the caller object until the complete
hydro transaction succeeds.

## Consequences

Unrelated ranks no longer allocate complete root fields or receive hydro
numerical broadcasts. Communication scales with nonlocal root tiles, distinct
child owners, and remote children while preserving serial ordering, parity,
and rollback.

The root kernel itself remains a single level-wide physics operation. Root
reconstruction/StateRedist decomposition, targeted transport traffic, public
time-loop control, regridding, checkpointing, and output remain future work.
