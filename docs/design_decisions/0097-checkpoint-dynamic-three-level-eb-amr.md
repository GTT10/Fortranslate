# Decision 0097: checkpoint dynamic three-level EB AMR

## Context

The public dynamic three-level lifecycle can replace its configured finest
seed with a tag-driven rectangle, but the static three-level checkpoint
requires that rectangle to match the namelist. Reusing that stream would
either discard the committed topology or weaken a previously qualified
compatibility contract.

## Decision

Keep the static three-level stream unchanged and select a distinct versioned
magic and schema for dynamic topology. Store the actual middle-to-finest
bounds, accepted regrid count, and every control that affects subsequent tag
planning in addition to the existing three-level mechanism, geometry,
physics, time, and field records.

On restart, require dynamic mode and matching physics and regrid controls.
Rebuild the configured fixed root-to-middle hierarchy, validate the stored
finest rectangle against its two-cell middle margin, and reconstruct that
actual child geometry. Read every level into private candidates, recover
temperatures through the EOS, validate the terminal marker, and publish the
complete hierarchy, step count, and regrid count as one transaction.

## Consequences

A split dynamic run resumes from the topology it committed rather than
repeating initialization from the seed. Restored accepted steps preserve the
`regrid_interval` schedule and restored regrid accounting remains visible to
the application. Checkpoint interval, stop-after-write, final time, maximum
steps, and output paths remain restart-adjustable.

The stream remains serial and formatted. Dynamic root or middle topology,
finest removal or sibling finest patches, parallel I/O, arbitrary depth,
molecular transport, and MPI-distributed ownership remain separate work.
