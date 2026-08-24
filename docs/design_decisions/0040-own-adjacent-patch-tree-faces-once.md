# Decision 0040: own each adjacent patch-tree face exactly once

## Context

Separated child patches can treat every side as a coarse/fine interface.
Independently owned adjacent children cannot: their touching side is a
fine/fine face. Filling both sides from the parent would give inconsistent
ghosts, and refluxing both registers would duplicate a correction on covered
parent cells.

## Decision

Allow adjacency only for patch sets constructed inside the patch-tree engine;
the older separated multipatch API retains its default rejection. Before each
child substep, fill parent-interpolated ghosts and then replace every ghost
fine index covered by a sibling with that sibling's interior state and
temperature. Apply the same lookup to all four PPM/WENO layers.

Let each child advance and return its time-integrated boundary flux. At every
touching pair, replace the left patch's right flux and the right patch's left
flux by their arithmetic mean. Correct the two interface-adjacent fine cells
conservatively so both patches reflect that one owned flux. Mark both register
sides as internal and zero them before coarse/fine reflux. Use the same rule
for hyperbolic and molecular-transport recursion.

## Consequences

Adjacent patches remain separate storage owners while composite conservation
does not depend on two independently calculated interface fluxes. Chains of
small adjacent boxes can source wide ghosts from more than one sibling by
global fine index. Existing separated patch behavior is unchanged.

The qualified implementation is same-process and reconciles the completed
local Runge--Kutta flux integrals. Globally stage-synchronous exchange between
distributed owners belongs to the MPI patch-ownership integration.
