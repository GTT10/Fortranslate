# Decision 0160: Plan EB tree tags serially

## Context

The arbitrary-depth reactive 2D EB tree can already replace its complete
topology transactionally, but callers must supply every branching rectangle
and EB child geometry. A solution-driven lifecycle needs to derive those
rectangles from accepted fields without coupling the numerical tree to one
level-set representation or mutating the accepted state during planning.

## Decision

Copy and deepest-first synchronize the accepted tree, then plan normalized
temperature-gradient tags independently for every prospective parent. Cluster
tags with the existing deterministic 2D collection planner and retain
parent-major child ordering. Ask a caller-supplied procedure to construct EB
geometry for each rectangle.

After every nonempty relation, initialize a temporary topology and PCM field
tree from the synchronized root before planning the next relation. Stop at the
level ceiling, when no tags remain, or when all parents are below the tagging
stencil extent. Compose the resulting complete plan with the existing atomic
overlap-preserving rebuild.

## Consequences

The serial arbitrary-depth EB tree can create, reshape, deepen, preserve, and
collapse branching topology from temperature tags. Geometry generation stays
application-owned, while topology, EOS, overlap, and conservation failures
retain the complete accepted tree. MPI ranks do not yet plan or broadcast this
metadata, recompute ownership, or migrate changed topology directly.
