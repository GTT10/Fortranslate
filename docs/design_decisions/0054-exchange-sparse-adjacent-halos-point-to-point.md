# Decision 0054: exchange sparse adjacent halos point to point

## Context

Sparse chemistry, hydro, and molecular transport retain only owner-local patch
arrays, but adjacent same-level ghost refresh broadcast every complete sibling
state and temperature array to all ranks. Only the two owners sharing an
adjacent face require that data, and only the reconstruction-width boundary
layers are consumed.

## Decision

Iterate adjacent sibling pairs in deterministic hierarchy order. Pack each
owner's state and temperature boundary into one payload per direction, using
one layer for narrow ghosts and the configured four layers for PPM. When the
owners differ, exchange both payloads with one `MPI_Sendrecv`; when they match,
copy the boundaries locally. Ranks unrelated to the pair skip allocation and
communication.

Update the narrow state/temperature ghost from payload layer one and update the
wide PPM ghost arrays from every transmitted layer. Count a completed
cross-owner pair once on its left owner so a communicator sum can be compared
with the independently derived adjacent-face count.

## Consequences

Adjacent sparse halo traffic scales with interface payload size and reaches
only the two participating owners. The blocking deterministic order retains
the existing collective failure boundaries and exact serial parity.

Parent interval states, child-to-parent synchronization, average-down, shared
flux reconciliation, and topology-changing overlap transfer still use the
collective correctness schedule and remain later point-to-point work.
