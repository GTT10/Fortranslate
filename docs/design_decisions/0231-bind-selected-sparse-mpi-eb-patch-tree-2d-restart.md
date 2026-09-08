# 0231: Bind selected sparse MPI EB patch-tree restart to its context

## Status

Accepted for `0.239.0`.

## Context

The `0.238.0` selected sparse MPI reactive EB patch-tree 2D application can
run the same dynamically regridded lifecycle as the fixed application, but it
rejects persistence. The fixed schema-8 checkpoint fingerprints the numerical
and physical model but does not identify the configure-time mechanism bundle,
generated chemistry integrator, or bundle-order composition. Reading such a
file under a different selected context would therefore be ambiguous.

The existing sparse I/O already gathers owner-local fields to one I/O root,
writes the serial patch-tree format, reconstructs distribution for the active
communicator on read, and scatters fields directly to their current owners.
That rank-neutral topology path should remain the sole distribution contract.

## Decision

Keep fixed magic `PELEF_REACTIVE_AMR_EB_PATCH_TREE_2D` and fixed schema 8
byte-for-byte. A call carrying the complete selected context exclusively uses
schema 9. After species order, schema 9 stores:

- `SELECTED_CONTEXT`, the 64-hex bundle SHA-256, generated `explicit` or
  `implicit` integrator, species count, and normalized bundle-order mole
  fractions;
- the existing complete schema-8 fingerprint, geometry, and topology;
- committed clock metadata followed by `COMPOSITE_BASELINE`, `nvar`, and the
  original composite conserved integrals;
- complete operator counters, regrid history, every state and recovered
  temperature field, and `END_CHECKPOINT` followed by strict end-of-stream.

Selected writes require all context, fingerprint, baseline, clock, counter,
and regrid-history arguments. They validate finite normalized composition,
positive density and energy, species closure, nonnegative raw species state,
positive temperature, and EOS recovery before opening the target. Selected
reads require schema 9 and validate the full private candidate before it can
be returned. Fixed readers continue to require schema 5 or 8 according to
their existing call contract and cannot interpret schema 9.

The sparse MPI adapter independently requires exact communicator-wide
agreement on context presence, bundle digest, integrator, and composition.
After that agreement, it collectively requires all eight selected writer
metadata arguments before gather and the selected reader fingerprint before
root I/O. No rank may return locally while peers enter gather or broadcast.
Only the selected root parses or writes the formatted file. A checkpoint
contains neither communicator size nor owner assignments; restart broadcasts
rank-neutral topology, computes a new distribution for the active
communicator, and scatters the validated fields to current owners. The shared
application explicitly selects the fixed or selected I/O call at both write
and restart boundaries.

Before mechanism loading or file creation, the selected front end rejects
lexical aliases among input, output, checkpoint, and restart paths.

## Consequences

A one-rank selected full-H2/O2 process can stop after its first committed
coarse step and resume on two or four ranks to the exact uninterrupted fixed
and selected final CSV bytes. The original composite baseline and cumulative
operator/regrid history continue across the process boundary.

The bundle SHA-256 is mechanism provenance, not authentication of the
checkpoint payload. Root-formatted `status="replace"` I/O is neither
crash-atomic nor scalable. Schema migration, fixed-depth MPI modes, runtime
mechanism loading, CFD CVODE, performance/thread qualification, detailed-fuel
validation, and external PeleC field parity remain outside this decision.
