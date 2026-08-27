# Decision 0173: Expose embedded-wall controls first in the single-level app

## Context

The `0.180.0` wall kernel accepted nondefault values through the shared
boundary-set API, but the public namelist always constructed its adiabatic slip
default. Adding inputs to checkpoint-capable AMR applications without extending
every restart compatibility format would allow a changed wall to resume an
incompatible solution silently.

## Decision

Add embedded-wall kind, thermal mode, temperature, and velocity to the existing
`&embedded_boundary` namelist. Validate that isothermal mode has thermal
transport, no-slip has viscous transport, and a moving wall is no-slip. Apply
the values transactionally in the public single-level EB driver, which has no
checkpoint/restart state.

Treat an isothermal or no-slip selection as an active nondefault wall. Reject
that selection in the public AMR driver preflight until its checkpoint and
fingerprint formats carry the same controls. Preserve the low-level AMR and
MPI boundary-set interfaces from `0.180.0`.

## Consequences

The installed single-level application can now exercise embedded heat and
momentum transfer from input alone, and invalid transport combinations fail
before state allocation or clock advancement. AMR applications do not ignore
parsed wall values. Public AMR activation and restart compatibility remain one
explicit follow-up milestone instead of an implicit format change.

