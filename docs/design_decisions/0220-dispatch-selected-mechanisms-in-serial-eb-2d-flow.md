# 0220: Dispatch selected mechanisms in serial EB 2D flow

## Status

Accepted for `0.227.0`.

## Context

The serial single-level reactive 2D embedded-boundary application already
supports plane and circle geometry, cut-cell hydro stabilization, masked
chemistry, molecular transport, and configurable embedded walls. Its program
unit nevertheless loaded only the committed elementary or full-H2/O2 tables
and directly owned geometry construction, boundary construction, simulation,
output, and diagnostics.

Linking a selected generated mechanism into `pelef_core` would permit module
and symbol collisions with committed generated mechanisms. Copying the EB time
loop into a selected program would create a second implementation. Selected
composition must also reach configured physical-boundary states, not only the
regular initial field, and selected chemistry policy must reach both masked
Strang half-steps without changing fixed callers.

## Decision

Extract `reactive_eb_2d_application_mod`. Fixed and selected front ends load
and validate their own thermo, reaction, transport, and composition data, then
call this driver for geometry and boundary construction, simulation,
deterministic EB CSV publication, volume-weighted invariants, and diagnostics.
Thread optional bundle-order mole fractions through the regular initializer
and reactive-boundary builder. Thread an optional chemistry-integrator policy
through the EB time loop and both active-cell chemistry stages. Absence
preserves fixed public behavior.

Build `pelef_selected_reactive_eb_2d_runtime` on the selected regular-2D
runtime in a private Fortran module directory. Compile only the generic EB
geometry, reconstruction, hydro, transport, driver, and application modules.
Configure and install `pelef_reactive_eb_2d_selected` from each normalized
bundle's declared interfaces and SHA-256.

The selected front end requires `chemistry_model = "selected"`, validates the
bundle and exact-name composition, checks every initial hotspot temperature
and any isothermal embedded-wall temperature against the common NASA7 range,
and rejects lexical input/output aliases before creating output. It rejects
AMR and checkpoint/restart controls because this selected application is
single-level and has no restart boundary.

## Consequences

Configure-time bundles with 2--32 species can use the existing serial
single-level 2D EB plane/circle geometry, PCM or characteristic-PLM hydro,
FluxRedist or StateRedist controls, native explicit/implicit chemistry,
molecular transport, and embedded-wall controls. A two-species fixture
requires active chemistry/transport, cut-cell masking, closure, and repeat-byte
identity. Pinned full-H2/O2 fixed/selected pairs require byte-exact output for
plane chemistry and characteristic-PLM circle transport with an isothermal
moving no-slip wall.

The application's four outer-domain boundaries remain outflow-only. Embedded
wall configuration is not a claim of general outer physical-boundary support.

The optional composition and integrator arguments preserve source
compatibility for callers rebuilt from source. They change Fortran module
procedure interfaces, so mixing new modules with old objects or archives is
not a supported binary ABI.

This milestone does not add selected EB AMR/3D, checkpoint/restart, MPI,
runtime mechanism loading, CVODE inside CFD, mechanisms above 32 species,
thread safety, performance qualification, detailed fuels, external PeleC
field parity, or experiment validation.
