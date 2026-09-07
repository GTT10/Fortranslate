# 0228: Bind selected multipatch EB AMR 2D restart to its context

## Status

Accepted for `0.236.0`.

## Context

The fixed serial dynamic two-level patch-set lifecycle already stores the
complete committed child-patch list, regrid history, clock, diagnostics, root
field, and every child field under
`PELEF_REACTIVE_EB_AMR_PATCH_SET_2D_CHECKPOINT` schema 3. That byte stream is
frozen. Selected calls were still rejected because schema 3 has no generated-
mechanism identity or original conservation baseline.

Reusing schema 3 would permit a selected restart under the wrong bundle,
integrator, or composition. Replacing the fixed layout would break existing
restart artifacts. A selected restart also needs the original composite-
integral baseline so resumed cumulative conservation diagnostics retain the
uninterrupted reference.

## Decision

Keep the fixed patch-set magic/schema 3 byte-for-byte. A selected patch-set
call uses schema 4 under the same magic. Immediately after the header it
writes:

- an all-or-none `SELECTED_CONTEXT` block containing the configure-time bundle
  SHA-256, generated `explicit` or `implicit` integrator, species count, and
  finite normalized bundle-order composition;
- a `COMPOSITE_BASELINE` block containing the finite `nvar`-component initial
  two-level composite integrals.

The existing schema-3 body then records species, geometry, wall and numerical
policies, dynamic controls, every committed child patch, clock and regrid
diagnostics, root and child state/temperature fields, and the terminal marker.
The reader derives the expected schema from the fixed/selected call contract.
It parses every record, recovered temperature, terminal marker, and end-of-
stream into private candidates before publishing any caller target.

The selected application admits this envelope only for the established serial
dynamic exactly-two-level patch-set path. Bundle-order composition and the
generated integrator are forwarded through every root and child stage. Before
mechanism loading or file creation, the front end derives all possible
`_patchNNNN` child-output paths and checks them for lexical aliasing with every
input, base/root output, checkpoint, and restart path.

## Consequences

A pinned full-H2/O2 process advances two disjoint child patches, stops after
its first coarse step, and resumes in another process to byte-exact
uninterrupted root and child CSV files. The selected uninterrupted files also
equal the fixed dynamic reference. Restoring the persisted baseline makes the
final cumulative conservation diagnostic textually identical. The existing
fixed regression independently pins the schema-3 checkpoint SHA-256 so this
increment cannot silently rewrite the frozen byte stream.

Changed composition, corrupt baseline, fixed/selected cross-schema input,
corrupt terminal marker, truncation, or trailing content publishes no
candidate. An invalid write leaves an existing checkpoint untouched. Derived
child-output aliasing is rejected before the mechanism bundle is loaded.

The selected schema records provenance but does not authenticate its payload.
Selected MPI EB AMR persistence, crash-atomic or scalable I/O, schema
migration, and payload authentication remain separate work.
