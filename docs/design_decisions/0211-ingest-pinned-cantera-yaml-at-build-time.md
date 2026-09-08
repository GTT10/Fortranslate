# 0211: Ingest one pinned Cantera YAML phase at build time

## Status

Accepted for `0.218.0`.

## Context

The full ten-species H2/O2 runtime path already had qualified NASA7,
third-body/Troe kinetics, Jacobian, implicit reactor, and gas transport, but
its normalized JSON, thermo table, and transport table were maintained as
separate hand-transcribed artifacts. That left source ordering, duplicate
reactions, named colliders, units, and provenance vulnerable to drift.

The pinned Cantera YAML contains both ideal-gas `ohmech` and Redlich-Kwong
`ohmech-RK` phases. A phase name and a deliberately bounded representable
subset are therefore part of the numerical contract.

## Decision

Commit the exact YAML source and select `ohmech` explicitly. Use Cantera 3.2
as a build/test-time parser to create a deterministic normalized JSON bundle.
Accept only ideal-gas NASA7 species at 101325 Pa with complete Lennard-Jones
gas transport and elementary, three-body, or Troe-falloff reactions. Preserve
element composition, source reaction index, duplicate flag, default and named
third-body efficiencies, selected phase, source hash, source metadata version,
runtime API version, and source/target units.

Generate the existing full-H2/O2 reaction/Jacobian kernels together with its
NASA7 and primitive transport loaders. Retain the public thermo and transport
interfaces as thin adapters so no application-level dispatch changes.

## Consequences

The committed YAML, normalized JSON, and generated Fortran form one
byte-checked chain, and the existing live Cantera trajectory/CFD tests guard
the conversion numerically. The elementary subset remains independent and
byte-identical.

This decision does not add runtime YAML parsing, arbitrary mechanism dispatch,
CHEMKIN, NASA9, PLOG, Chebyshev, SRI, Lindemann-only falloff, custom reaction
orders, mechanisms above 32 species, detailed fuels, or a CVODE-class
production integrator. Generator real literals remain canonicalized with
twelve digits after the decimal in exponential notation.
