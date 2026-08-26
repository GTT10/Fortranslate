# Decision 0123: connect root-only EB materialization to existing writers

## Context

The sparse MPI EB hierarchy can be gathered without allocating complete fields
on non-root ranks, but the checkpoint and CSV lifecycle still has to preserve
the established serial formats and report filesystem failure consistently to
all participants.

## Decision

Add a dedicated MPI I/O layer above sparse ownership and below applications.
For checkpoint and CSV publication, gather each numerical entity only to a
caller-selected writer, invoke the existing serial writer only there, then
broadcast its logical result. Publish sender-local gather counts only after the
write result is accepted.

Keep the formatted multipatch checkpoint schema unchanged. For CSV output,
write the configured root path and derive each child path with the established
`_patchNNNN` convention. Writer-only path, configuration, time, and checkpoint
metadata are authoritative; numerical ownership and root selection remain
collective preconditions.

## Consequences

Existing analysis and serial checkpoint readers remain compatible, while
non-root temporary memory stays proportional to owned entities. Open or write
failure becomes visible on every rank and publishes zero accounting.

The formatted checkpoint is still read as a complete hierarchy and is not yet
redistributed directly to a new sparse owner map. CSV publication is a set of
ordinary files rather than an atomic parallel format, so a later-child failure
can leave earlier files on disk.
