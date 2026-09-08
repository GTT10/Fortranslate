# 0219: Dispatch selected mechanisms in serial AMR 1D flow

## Status

Accepted for `0.226.0`.

## Context

The serial reactive 1D AMR application already supports solution-driven
two-level regridding, arbitrary-depth hierarchies, dynamic multipatch level
one, characteristic PPM, native chemistry, and molecular transport. Its
program unit nevertheless loaded only the committed elementary or full-H2/O2
tables and owned the simulation, output, and diagnostics directly.

Linking a selected generated mechanism into `pelef_core` would permit module
and symbol collisions with the committed generated mechanisms. Copying the AMR
algorithms into a selected program would create two hierarchy implementations.
The selected composition must instead reach the existing root initializer in
all three public AMR modes without changing existing fixed callers.

## Decision

Extract `amr_reactive_1d_application_mod`. Fixed and selected front ends load
and validate their own mechanism and composition, then call this driver for
mode dispatch, simulation, deterministic composite CSV publication,
conservation diagnostics, and stdout reporting. Thread an optional
bundle-order base mole-fraction array through two-level, multilevel, and
multipatch initialization and simulation entry points; absence preserves the
fixed public behavior.

Build `pelef_selected_amr_reactive_1d_runtime` on the selected regular-1D
runtime in a private Fortran module directory. Compile only the generic 1D AMR
hierarchy, regrid, reacting evolution, and application modules. Configure and
install `pelef_amr_reactive_1d_selected` from each normalized bundle's declared
interfaces and SHA-256.

The selected front end requires `chemistry_model = "selected"` and
`amr_enabled`, validates bundle structure, exact-name composition, transport
ordering, and every initial temperature implied by a hotspot or constant-
pressure entropy wave before creating output. It rejects checkpoint/restart
controls because this existing application does not implement them.

Require every complete schema-1 bundle to declare its native flow-chemistry
policy as `chemistry_integrator = "explicit"` or `"implicit"`. Emit that value
as a generated public constant and pass it explicitly through every selected
regular 1D, AMR 1D, 2D, and 3D chemistry stage. Fixed callers omit the new
optional trailing argument and retain their established policy. This removes
species count as the selected-mechanism dispatch key.

## Consequences

Configure-time bundles with 2--32 species can use the existing serial 1D
two-level, arbitrary-depth, or dynamic multipatch AMR algorithms with native
chemistry and molecular transport. A two-species fixture requires active
chemistry, regridding, closure, composite-domain coverage, and repeat-byte
identity. Pinned full-H2/O2 fixed/selected pairs require byte-exact output for
active two-level chemistry/transport, three-level characteristic PPM, and
three-patch molecular transport. A three-level outflow case requires the
finest level to touch the physical lower boundary while characteristic PPM and
hybrid WENO7-Z are active.

The optional composition and integrator arguments preserve source
compatibility for callers rebuilt from source. They change Fortran module
procedure interfaces, so mixing new modules with old objects or archives is
not a supported binary ABI.

This milestone does not add selected AMR checkpoint/restart, EB or MPI
dispatch, runtime mechanism loading, CVODE inside CFD, mechanisms above 32
species, thread safety, performance qualification, detailed fuels, external
PeleC field parity, or experiment validation.
