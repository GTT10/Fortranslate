# 0225: Bind selected dynamic EB AMR 2D restart to its topology class

## Status

Accepted for `0.233.0`.

## Context

The selected serial reactive EB AMR 2D application gained static two-level
checkpoint/restart in `0.232.0`. Fixed calls retained schema 3 and selected
static calls received schema 4 with generated bundle, integrator, and
composition context. The shared driver already supports a serial dynamic
two-level hierarchy with one fine patch, but schema 4 was explicitly frozen to
the static selected contract.

Reusing schema 4 for a dynamic selected checkpoint would broaden a released
format after qualification and would not make the static/dynamic topology
class explicit at the envelope boundary. Changing schema 3 would break the
frozen fixed byte stream.

## Decision

Dispatch three exclusive checkpoint envelopes in
`reactive_eb_amr_2d_driver_mod`:

- fixed calls omit selected context and retain schema 3 byte-for-byte;
- selected static calls retain schema 4 byte-for-byte;
- selected dynamic calls use schema 5.

Schema 5 carries the same all-or-none `SELECTED_CONTEXT` block as schema 4,
followed by a `DYNAMIC_BASELINE` record containing the initial composite
integrals. The established checkpoint body then records the
dynamic-regridding policy, the actual fine-patch bounds at the checkpoint,
regrid cadence and count, clock, diagnostics, geometry, and coarse/fine
fields. The writer chooses the
schema only after validating selected context. The reader derives the single
expected schema from the caller's selected/fixed and static/dynamic contract;
static/schema-5 and dynamic/schema-4 reads fail without fallback.

The selected front end admits dynamic regridding only for the existing serial,
exactly-two-level, single-fine-patch hierarchy. Three-level, multipatch, and
dynamic-parent modes remain rejected. Existing lexical alias checks cover the
input, both outputs, checkpoint, and restart paths before persistence begins.

## Consequences

The qualified full-H2/O2 case starts with coarse patch bounds `(2:5,2:5)`,
regrids once to `(5:12,4:12)` before writing its first-step schema-5
checkpoint, and resumes in a separate process to the exact uninterrupted
three-step coarse and fine CSV files. Those uninterrupted selected files also
match the fixed dynamic reference byte-for-byte, and the resumed cumulative
conservation diagnostic equals the uninterrupted value. A malformed baseline,
changed selected composition, and either static/dynamic cross-schema direction
fail transactionally.

Schemas 3 and 4 remain frozen. Schema 5 establishes provenance and topology
compatibility, not a cryptographic digest of its body. Three-level,
multipatch, dynamic-parent, MPI EB AMR, crash-atomic or scalable I/O, schema
migration, and payload authentication remain separate work.
