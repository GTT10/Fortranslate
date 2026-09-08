# 0216: Dispatch selected mechanisms in regular 1D flow

## Status

Accepted for `0.223.0`.

## Context

The configure-time mechanism boundary could generate, compile, probe, and run
an arbitrary normalized schema-1 bundle in a constant-volume reactor, but all
CFD applications still linked the committed elementary or full-H2/O2 modules.
Directly linking a selected bundle into `pelef_core` is unsafe because a valid
bundle may intentionally declare the same Fortran module and symbol names as
the committed generated full-H2/O2 source.

The regular 1D solver also obtained the transport derived type through a fixed
database module that imports both committed mechanism loaders. That dependency
would recreate the collision even if the selected application never called a
fixed loader. Initial composition was hard-coded by chemistry model rather
than supplied in bundle order.

## Decision

Build a separate `pelef_selected_reactive_1d_runtime` whenever tests or a user
selected bundle require it. Compile only the regular 1D state layout,
configuration, reconstruction, transport, and evolution sources into a private
Fortran module directory. Link that runtime to
`pelef_selected_runtime_core`, not to `pelef_core`.

Move the mechanism-independent `gas_transport_species` type, record validity,
compatibility checks, and array-to-database builder into `gas_transport_mod`.
Let `transport_database_mod` reexport that interface while retaining the fixed
elementary/full-H2/O2 loaders. The selected application builds the same
transport database directly from its generated primitive arrays.

Extend the regular 1D configuration with a bounded exact-name composition and
an explicit `allow_selected` parser switch. The fixed application leaves that
switch false and rejects `chemistry_model = "selected"`. The configured
selected application enables it, maps nonduplicate names into bundle order,
and supplies the normalized mole fractions through optional arguments on the
regular 1D initialization/evolution interfaces. Existing fixed callers remain
source-compatible.

Configure `pelef_reactive_1d_selected` from the bundle's declared loaders and
embed the normalized JSON SHA-256. Validate reaction structure and
molecular-mass balance, transport records, the common NASA7 interval, and the
requested temperature range before output creation. Reject AMR and
checkpoint/restart controls because this executable is a regular-grid serial
application.

## Consequences

A configure-time selected mechanism can now use the existing conservative 1D
hydrodynamics, native chemistry splitting, and molecular transport without a
generated-module collision. The pinned selected full-H2/O2 case must remain
byte-identical to the fixed application, while an independently ordered
two-species fixture proves that the path is not just a fixed-name alias.

There is intentional code recompilation in the small selected 1D runtime. This
is preferable to renaming generated modules or adding fragile link-order
assumptions, and its parity is guarded in every test configuration.

This decision does not introduce runtime YAML/JSON loading, CVODE integration
inside CFD, selected 2D/3D/AMR/EB/MPI applications, mechanisms above 32
species, thread safety, performance qualification, detailed-fuel validation,
or external PeleC/experiment parity.
