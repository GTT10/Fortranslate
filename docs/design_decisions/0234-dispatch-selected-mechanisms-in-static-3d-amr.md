# 0234: Dispatch selected mechanisms in static 3D AMR

## Status

Accepted for `0.242.0`.

## Context

The `0.241.0` static, periodic, two-level 3D AMR applications qualified
transactional fixed elementary and full-H2/O2 chemistry in serial and
rank-local sparse MPI execution. Configure-time selected mechanisms already
served regular 3D flow, but the AMR frontends owned fixed mechanism loading,
composition resolution, lifecycle, diagnostics, and MPI startup as one
application-specific path. Reusing that monolith would duplicate AMR state
publication and would allow ranks to initialize from different generated
mechanism provenance.

## Decision

Move serial and MPI static 3D AMR state ownership into mechanism-independent
application drivers and add generated selected-mechanism frontends.

- Frontends read configuration, load and validate generated NASA7, ordered
  reactions, and transport provenance, resolve the exact species-name
  composition, validate the requested temperature range, and pass the
  generated chemistry integrator explicitly.
- The fixed frontends use the same drivers with their established mechanism
  databases and composition mapping. Their public chemistry and
  characteristic-PLM output hashes remain unchanged.
- `thermo_model='selected'` is an explicit opt-in to the AMR configuration
  reader. The fixed frontend still rejects it.
- Before root allocation or state initialization, the sparse MPI driver
  requires exact rank consensus for bundle SHA-256, integrator, composition,
  ordered NASA7 records, ordered reactions including pressure-dependent data,
  and transport records. A testable validator returns collective failure
  without publishing state.
- Ranks with no fine planes remain participants in every consensus, source,
  and candidate-acceptance collective.
- H/O/N diagnostics run only when every species belongs to the qualified
  H2/O2/N2/AR name set. Other valid selected mechanisms retain Euler,
  species-closure, physical-state, and synchronization gates and report
  elemental diagnostics as unavailable.
- Selected static 3D AMR molecular transport and checkpoint/restart remain
  startup errors. The fixed schema-2 restart contract is unchanged and is not
  relabeled as selected-mechanism provenance.

## Consequences

The generated two-species fixture is deterministic in serial and at one, two,
and four MPI ranks, with the four-rank case exercising zero-fine ranks. The
generated full-H2/O2 executable is byte-identical to the fixed executable and
to one-, two-, four-, and eight-rank MPI output at both levels. Direct MPI
tests reject valid-but-different SHA, integrator, composition, NASA7,
reaction equation/rate, and transport contexts.

This decision does not add AMR molecular transport, selected checkpoint
schema, dynamic topology, physical coarse boundaries, 3D EB AMR, distributed
I/O, performance qualification, or external PeleC field parity.
