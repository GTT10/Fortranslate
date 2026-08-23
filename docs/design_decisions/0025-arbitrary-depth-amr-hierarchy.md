# Decision 0025: represent arbitrary AMR depth as adjacent-level relations

## Status

Accepted for PeleF `0.30.0`.

## Context

The first AMR application hard-coded one root and one fine patch. Prolongation,
restriction, flux registers, and reflux were already correct for an adjacent
pair, so duplicating those algorithms for a fixed third level would add another
limit instead of satisfying the arbitrary-level requirement.

## Decision

Keep the qualified two-level relation as the unit of coarse/fine geometry and
store an allocatable sequence of relations for a nested single-patch hierarchy.
Each relation has its own refinement ratio and parent-local bounds. Store each
level field in an allocatable component so different patch sizes are natural.

Prolong from root to deepest level. Synchronize in reverse order: reflux an
interface, average its child onto the covered parent cells, then continue to
the next coarser interface. Define the composite integral as every noncovered
parent region plus the complete deepest level.

## Consequences

- hierarchy depth and refinement ratios are runtime data with no fixed cap;
- the verified two-level algorithms remain the adjacent-pair implementation;
- mixed refinement ratios produce cumulative subcycle schedules explicitly;
- deepest-to-root ordering carries finer information into every coarser level;
- reactive state ownership, recursive operator advancement, recursive
  regridding, composite output, and multiple patches remain separate follow-on
  integrations.
