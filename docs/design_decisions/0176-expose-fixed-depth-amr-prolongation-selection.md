# Decision 0176: Expose fixed-depth AMR prolongation selection

## Context

The EB-safe limited-linear kernel is conservative and qualified in isolation,
but every public reactive 2D EB AMR lifecycle still invokes PCM. Selecting the
new kernel must cover static initialization and later regrids without silently
changing checkpoint compatibility.

## Decision

Add `prolongation_method` to `&eb_amr`, accepting `pcm` and `linear` and
retaining PCM as the default. Route both values through one dispatcher and
propagate the selection through two-level, sibling-patch, and three-level
fixed-depth initialization and regridding. Select linear in the installed
hot-wall transport case.

Do not serialize an incomplete compatibility contract. Reject linear mode when
a fixed-depth checkpoint or restart path is requested, and reject it at the
arbitrary-depth serial and sparse-MPI boundary until the shared fingerprint
records the method.

## Consequences

Fresh fixed-depth applications can use the qualified higher-order initializer
without changing the robust default or the low-level PCM API. Every current
checkpoint remains unambiguous and transactional. Linear restart and
arbitrary-depth support require a later schema/fingerprint milestone.
