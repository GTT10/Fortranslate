# 0237: Persist static 3D AMR transport checkpoints

## Status

Accepted for `0.245.0`.

## Context

The [0236 decision](0236-add-static-3d-amr-molecular-transport.md) added
mixture molecular transport to the static, strictly-interior, periodic,
two-level Cartesian 3D AMR hierarchy, but deliberately rejected transport
checkpoint and restart. The fixed schema 2 and selected schema 3 bodies do not
identify the complete gas-transport database, parabolic controls, transport
operator, or cumulative transport diagnostics. Reading either schema for an
active transport run could therefore publish state under a different
transport operator or report only a restart suffix.

The next increment needs to preserve the established transport-disabled
formats while making transport restart rank-neutral, context-bound, and
transactional. It remains a static two-level single-patch decision; it does
not turn formatted root I/O into a crash-safe or scalable persistence layer.

## Decision

Keep the existing static-3D-AMR magic and assign exclusive transport schemas:

| Call mode | Transport disabled | Transport enabled |
|---|---:|---:|
| Fixed | schema 2, unchanged | schema 4 |
| Selected | schema 3, unchanged | schema 5 |

Schema 2 and schema 3 remain readable only through their existing
transport-disabled contracts. Schema 4 requires the fixed public transport
case (`full_h2o2`, `chemistry_enabled=.false.`, and
`transport_enabled=.true.`). Schema 5 requires selected chemistry context
(`thermo_model='selected'`, `chemistry_enabled=.true.`, and
`transport_enabled=.true.`) together with the complete selected context.
Fixed/selected and transport-enabled/transport-disabled cross-schema reads
are rejected without fallback or schema reinterpretation. A 0.244 reader is
not expected to consume schemas 4 or 5.

### Transport context

Schemas 4 and 5 record and validate the complete ordered
`gas_transport_species` database for the checkpoint species. Each record
contains its name, geometry, Lennard--Jones well depth and diameter, dipole,
polarizability, and rotational relaxation value. The serialized context also
binds the parameter convention
`EPSILON_OVER_K_K;SIGMA_ANGSTROM;DIPOLE_DEBYE;POLARIZABILITY_ANGSTROM3;ROT_RELAX_DIMENSIONLESS`
and the operator identity `STATIC_AMR_3D_R_T_H_T_R_V1`.

The transport controls are stored alongside that database: transport
enablement, viscosity, thermal conduction, species diffusion, barodiffusion,
and `transport_cfl`. The existing body continues to bind the numerical,
solver, mesh, patch, physics-flag, and clock records. A transport checkpoint
also stores the cumulative maximum transport diffusivity, minimum accepted
transport limiter theta, maximum reflux correction, and original composite
integral baseline. Nonfinite, out-of-range, missing, reordered, or mismatched
records and controls are rejected before transport state is accepted.

Selected schema 5 additionally retains the schema-3 selected context: bundle
SHA-256, generated chemistry integrator, bundle-order composition, complete
ordered reactions, and chemistry tolerances. The bundle digest remains
provenance; this decision does not claim cryptographic payload authentication.

### Checkpoint boundary and transactional read

Write only at `POST_ACCEPTED_COARSE_STEP`, after the complete committed
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` step, including transport reflux,
average-down, temperature recovery, and diagnostic accumulation. The phase
marker is part of the serialized transport context and must match on read.

Read all context, metadata, level fields, recovered temperatures, the terminal
marker, and strict end-of-stream into private candidates. Publish caller-owned
levels, temperatures, time, step, baseline, reflux history, and transport
diagnostics only after every check and the final close succeeds. Any failed
read leaves the caller targets unchanged.

### MPI ownership and continuation

Before root allocation, transport payload communication, or root I/O, every
rank participates in all-or-none consensus for selected-context presence,
species and transport dimensions, ordered transport records, transport
controls, operator policy, and selected chemistry records where applicable.
Ranks with no fine planes remain participants and must not dereference empty
fine storage.

The checkpoint remains a root-formatted, rank-neutral representation. It stores
no communicator size or ownership map. Rank zero alone reads or writes the
file; a successful root operation broadcasts its status and restored metadata
before the communicator scatters the hierarchy to ownership derived from the
current communicator. A checkpoint written with one rank count can therefore
be restarted with a changed rank count, subject to the same static topology
and context contract. A root read/write failure is collective and no rank may
continue to output or advance after the failure.

### I/O boundary

Validation and publication are transactional at the application/API level,
but formatted replacement with `status='replace'` is retained. This decision
does not claim crash-atomic replacement, durable commit, payload
authentication, or scalable/distributed I/O. Those are separate decisions.

## Consequences

The old fixed schema-2 and selected schema-3 transport-disabled restart
contracts remain isolated and available for compatibility. Active transport
restarts use schemas 4 and 5 and cannot silently omit the database, controls,
operator identity, or cumulative diagnostics. Rank-neutral root I/O permits a
changed-rank continuation while preserving the owner-local stepping model.

The scope remains one static, strictly-interior, periodic, two-level Cartesian
hierarchy with one rectangular fine patch and the established mixture-averaged
transport model. Dynamic topology or regridding, boundary-touching patches,
physical coarse/fine boundaries, 3D EB AMR, scalable I/O, crash-atomic
replacement, payload authentication, performance qualification, and external
PeleC field parity remain outside this decision.

This ADR defines the persistence contract only. Validation counts, hashes,
installation results, and any release-complete statement belong in a later
validation record rather than here. The current 0.244 qualification boundary
is recorded in [validation/0.244.0.md](../validation/0.244.0.md).
The focused 0.245 qualification is recorded in
[validation/0.245.0.md](../validation/0.245.0.md).
