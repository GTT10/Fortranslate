# 0212: Compile selected normalized mechanism bundles as isolated targets

## Status

Accepted for `0.219.0`.

## Context

Milestone `0.218.0` established a deterministic YAML-to-JSON-to-Fortran chain
for the pinned full H2/O2 mechanism. The generated module was still fixed in
the core source list, so a second valid normalized bundle could not prove that
its declared interfaces compiled and linked without editing production source.
Source provenance also needed to remain checked if a user supplied the YAML
alongside an already normalized bundle.

## Decision

Add optional `PELEF_MECHANISM_BUNDLE` and `PELEF_MECHANISM_SOURCE` CMake cache
paths. Read the bundle's declared module, loader, kernel, Jacobian, thermo,
transport, and symbol-prefix names during configuration. Regenerate the
Fortran source into the build tree, compile it in an isolated static library
and Fortran module directory, and link an installable probe configured from
those declared names.

When a source path is supplied, require its basename and SHA-256 to equal the
bundle provenance both during configuration and immediately before generated
source publication. Track both inputs as CMake configure dependencies so a
plain rebuild refreshes the interface-name and expected-hash snapshots before
regeneration. The probe loads every declared table, checks species and
transport ordering and admissibility, and executes finite production-rate and
Jacobian evaluations. Kinetic activity is a fixture policy rather than an ABI
requirement because a representable mechanism may legitimately be inactive at
the generic probe state.

Keep the selected target outside `pelef_core`. The ordinary executable and
install sets therefore remain unchanged unless selection is explicit. Compile
a distinct two-species fixture in every test configuration, register a
selected user bundle in CTest when tests are enabled, and qualify the pinned
full H2/O2 bundle in a clean tests-disabled configuration. Also strengthen the
existing full-reactor Jacobian test with an energy-constrained centered
directional finite-difference comparison. Require generated identifiers to be
valid, noncolliding Fortran names and split a legal 128-character reaction
equation across bounded-width continuation lines.

## Consequences

Normalized bundles now have a configure-time compile, link, runtime, source-
pinning, and install boundary independent of the fixed H2/O2 adapter. Generated
module-name collisions do not enter the core target because each selected
module has a separate target and module directory.

This decision does not add runtime YAML/JSON parsing, application-level
mechanism dispatch, CHEMKIN, mechanisms above 32 species, new thermodynamic or
reaction families, detailed-fuel qualification, a CVODE-class production
integrator, or PeleC field parity. If no source path is supplied, the bundle is
accepted as the provenance authority; its embedded source hash is not evidence
that an external source file was rechecked.
