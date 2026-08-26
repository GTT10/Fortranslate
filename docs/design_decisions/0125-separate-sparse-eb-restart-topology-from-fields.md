# Decision 0125: separate sparse EB restart topology from fields

## Context

Root-only checkpoint read and direct scatter avoid broadcasting checkpoint
fields, but their first interface still accepts a full reactive patch set as a
replicated topology template. That type necessarily allocates conserved state
and temperature for every child, preserving a rank-multiplied numerical-field
copy outside the checkpoint reader.

## Decision

Introduce a geometry-only reactive EB patch topology. Each ordered child stores
its EB geometry and coarse/fine patch box, with no conserved state or
temperature components. Validate owner distributions and rank-local sparse
payload shapes directly against that descriptor. Construct it directly from
fine geometries plus an ordered regrid-plan collection, or extract it from an
established full patch set for compatibility.

Add geometry-only direct-scatter and root-only checkpoint-read entrypoints. The
selected root must provide the complete read fields and matching patch set;
non-root ranks provide empty complete inputs. Keep the former full patch-set
entrypoints as compatibility wrappers that extract the new descriptor before
entering the collective operation.

## Consequences

Sparse restart no longer requires replicated child numerical fields. The
formatted checkpoint schema, one-send-per-remote-entity schedule, metadata
broadcast, and transactional failure behavior remain unchanged.

EB geometry arrays and patch boxes remain replicated. The writer/gather and
physics/regrid compatibility boundaries still accept full patch-set templates;
checkpoint-driven topology reconstruction, their geometry-only migration,
arbitrary-depth dynamic EB topology, and decomposition of the level-wide root
physics kernel remain subsequent work.
