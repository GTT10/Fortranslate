# Decision 0122: gather sparse EB output only to root

## Context

Sparse MPI EB physics and regridding keep numerical fields globally
single-copy, but legacy checkpoint and output consumers need a complete ordered
hierarchy. The original materialization boundary broadcasts every owner field
to every rank, multiplying temporary memory by the MPI rank count even when
only one writer will consume it.

## Decision

Add a second collective materialization API with a caller-selected root. Pack
state and temperature together for each root tile and fine child. If its owner
is not the selected root, send that entity exactly once directly to the root.
Only the root allocates the complete coarse arrays and copies the replicated
patch template before filling every child field.

Require every rank to agree on a valid root and valid sparse input before any
payload transfer. Validate the complete reconstruction on the root and publish
outputs and sender-local traffic counts only after communicator-wide
acceptance. On rejection, every complete array remains unallocated and the
returned patch set is empty on every rank.

## Consequences

Serial checkpoint and CSV writers can be wrapped without an all-rank numerical
replica. Traffic is proportional to remote owned entities, and non-root memory
does not grow with total hierarchy size. The existing all-rank materialization
API remains for legacy replicated operators.

This decision supplies the transfer boundary only. Wiring the formatted
checkpoint and CSV lifecycle, rebuilding replicated geometry on restart, and
parallel file formats remain subsequent work.
