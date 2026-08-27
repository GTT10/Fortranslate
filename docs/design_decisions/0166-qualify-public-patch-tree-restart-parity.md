# Decision 0166: Qualify public patch-tree restart parity

## Context

The arbitrary-depth reactive EB tree already had self-describing serial and
sparse-MPI checkpoint APIs, and the public serial application could invoke
them. Its process-level continuation path had not yet been compared with an
uninterrupted application run.

## Decision

Qualify the installed serial executable with three independent invocations:
an uninterrupted reference, a run that stops after its first scheduled
checkpoint, and a continuation that reads that file. Use the same four-level
temperature-tagged case and retain periodic regridding so the restored global
step controls the next topology decision.

Compare final composite output by `(level, patch, i, j)` rather than row
position. Require identical topology and columns before comparing every
numeric field. Inspect the checkpoint envelope and require the stopped output
to precede the configured final time.

## Consequences

The public serial lifecycle now demonstrates process-boundary numerical parity
for arbitrary-depth geometry, state, and counters. Continuation still takes
physics and tagging controls from the restart input, so compatibility hashes
and policy for intentionally changed controls remain future work. A public
sparse-MPI executable remains separate from the already qualified sparse-MPI
checkpoint API.
