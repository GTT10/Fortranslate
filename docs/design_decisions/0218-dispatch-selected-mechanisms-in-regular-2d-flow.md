# 0218: Dispatch selected mechanisms in regular 2D flow

## Status

Accepted for `0.225.0`.

## Context

The regular 2D solver already owns the project's broadest single-level
reacting-flow boundary: PCM and characteristic PLM/PPM, optional CTU, native
chemistry, molecular transport, and periodic or physical domain boundaries.
Its public executable nevertheless loaded only the committed elementary or
full-H2/O2 tables and also duplicated simulation, output, and diagnostics in
the program unit.

The selected 3D runtime compiled copies of several 2D modules because its face
and transport implementation depends on them. Adding a second selected 2D
front end on that 3D target would work, but would expose unrelated 3D modules
and preserve duplicate ownership. Linking a generated bundle into
`pelef_core` would permit module-name collisions with the committed generated
mechanisms.

Physical walls introduce two additional selected-mechanism contracts. A
prescribed species flux has bounded input storage whose unused tail must not
be silently ignored, and an isothermal ghost fill can evaluate NASA7 data at
the reflected temperature `max(tiny, 2*Twall - Tinterior)` rather than only at
the wall temperature.

## Decision

Extract `reactive_2d_application_mod`. Fixed and selected front ends load and
validate their own mechanism and composition, then call the same driver for
simulation, deterministic CSV publication, invariant/extrema evaluation, and
diagnostics. Thread an optional bundle-order base composition through regular
2D initialization, the diagonal composition wave, boundary construction, and
the simulation loop while preserving fixed callers when it is absent.

Add opt-in selected fields to `simulation_config_reactive_2d_mod`. The fixed
parser continues to reject `chemistry_model = "selected"`. The configured
selected front end resolves exact species names after bundle loading, validates
all actual-species wall fluxes and their zero tails, and checks initial,
hotspot, overlapping-double-hotspot, thermal-channel, isothermal-wall, and
reflected initial ghost temperatures against the common NASA7 interval before
output creation.

Build `pelef_selected_reactive_2d_runtime` in its own Fortran module directory
on top of the selected 1D runtime. Privately compile only the regular 2D mesh,
configuration, boundary, transport, evolution, and application modules. Make
the selected 3D runtime link this target and remove its duplicate 2D sources.
Configure and install `pelef_reactive_2d_selected` from each normalized
bundle's declared interfaces and SHA-256.

## Consequences

Configure-time selected bundles with 2--32 species can use the existing serial
single-level regular 2D periodic and physical-boundary algorithms without a
copied time loop or a fixed generated-mechanism dependency. The fixture bundle
declares and exercises the generic adaptive explicit chemistry route. The
pinned full-H2/O2 bundle declares and exercises the generic adaptive implicit route,
characteristic PLM with CTU and active chemistry, and prescribed-wall species
transport. In all three full-mechanism regimes, fixed and selected output must
remain byte-identical.

The milestone does not dispatch selected mechanisms through AMR, EB, or MPI;
does not load mechanisms at runtime; and does not place CVODE inside a flow
solver. The 32-species limit, sequential native integration, lack of production
performance evidence, and lack of detailed-fuel or external physical
validation remain explicit boundaries.
