# Decision 0062: weight MPI AMR ownership by level subcycle work

## Context

The sparse MPI AMR owner map greedily balanced stored cell counts. A fine cell
is nevertheless advanced more often than a root cell: cumulative refinement
ratio `r` times for hyperbolic physics and `r^2` times for parabolic transport.
Cell-only balance can therefore place substantially different compute work on
ranks even when their stored-cell counts are similar.

## Decision

Store a 64-bit work count for every patch and rank. Define patch work as its
cell count times the cumulative refinement ratio raised to a selected exponent:

- exponent 0 preserves cell/storage weighting;
- exponent 1 estimates hyperbolic subcycle work;
- exponent 2 estimates parabolic subcycle work.

Assign patches in deterministic hierarchy order to the rank with the lowest
current work total, retaining the lowest-rank tie break. Require every rank to
agree on the exponent and reject values outside 0 through 2. Validate cell,
patch, and work totals together. Explicit-plan and owner-local tag-driven
sparse regrids inherit the old distribution's exponent.

## Consequences

Callers can select a cost model matching the dominant operator without changing
field ownership or communication APIs. The four-level qualification hierarchy
reduces maximum estimated parabolic work for two and four ranks while never
increasing it for one or eight ranks. This remains a deterministic static
greedy schedule; measured runtime feedback, device heterogeneity, and dynamic
work stealing remain future work.
