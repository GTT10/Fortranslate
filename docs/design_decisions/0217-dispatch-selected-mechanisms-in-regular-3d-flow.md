# 0217: Dispatch selected mechanisms in regular 3D flow

## Status

Accepted for `0.224.0`.

## Context

The `0.223.0` boundary dispatched a configure-time normalized mechanism through
the regular 1D solver, but the periodic regular 3D application still loaded the
committed full-H2/O2 mechanism directly. Copying that application and replacing
only its loader would duplicate the timestep loop, invariant publication, and
diagnostics. Linking a selected bundle into `pelef_core` would also reintroduce
the generated-module collision that the selected runtime was designed to avoid.

The 1D selected front end contained mechanism-structure, mass-balance,
transport-order, common-temperature-range, and exact-name composition checks
that are dimension independent. Repeating those checks in each future selected
application would create divergent startup contracts.

## Decision

Move selected composition resolution into `selected_composition_mod` and
selected thermo/reaction/transport startup checks into
`selected_mechanism_runtime_mod`. The regular 1D and 3D selected front ends use
the same bounded 2--32-species contract. Fixed parsers continue to reject the
`selected` model unless an explicitly configured selected front end enables it.

Extract the regular 3D application orchestration into
`reactive_3d_application_mod`. It receives validated species, reactions,
transport data, and bundle-order base composition, then owns mesh creation,
problem initialization, timestep selection, transactional `R-T-H-T-R`
advancement, invariant checks, deterministic CSV output, and diagnostics. The
fixed and selected front ends differ only in loading, validation, provenance,
and composition resolution.

Build `pelef_selected_reactive_3d_runtime` in a private Fortran module directory.
It reuses the selected 1D thermochemical/regular-face runtime and privately
compiles the multidimensional and 3D regular-grid modules. Generic transport
consumers import `gas_transport_mod`, so neither fixed transport loaders nor
committed generated mechanisms enter the selected link graph. Configure and
install `pelef_reactive_3d_selected` from the bundle's declared loader symbols
and normalized JSON SHA-256.

## Consequences

Configure-time selected mechanisms can now run the existing periodic
regular-grid 3D PCM or characteristic-PLM hydro, native chemistry, and molecular
transport path without symbol collision or a copied application loop. The
pinned selected full-H2/O2 transport case must remain byte-identical to the
fixed front end. An independently ordered two-species fixture proves active
chemistry and exact-name dispatch rather than a fixed-name alias. A short
nonuniform full-H2/O2 characteristic-PLM hotspot additionally keeps the full
selected reaction family and higher-order branch live with exact fixed/selected
output.

The selected 3D runtime intentionally recompiles its bounded dependency graph.
The resulting executable remains serial, periodic, single level, native-reactor
only, and limited to complete schema-1 bundles with 2--32 species. This decision
does not add physical boundaries, selected 2D/AMR/EB/MPI, runtime YAML/JSON
loading, CVODE inside CFD, thread safety, performance qualification,
detailed-fuel validation, or external PeleC/experiment parity.
