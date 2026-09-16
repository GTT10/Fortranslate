# Regular 3D AMR patch tree

The existing `amr_patch_3d` describes a single child box in parent-cell indices.
`amr_patch_tree_3d_mod` composes those boxes into a tree instead of replacing the
qualified two-level kernels. Geometry/state operations are separate from time
integration so that multi-level conservation can be checked before adding flux
and subcycling logic.

## Representation and invariants

Node 1 is the root (parent 0, level 0). Every other node names an earlier parent;
`patches(i)` and `parents(i)` describe node `i+1`. Parent-before-child order is
required, but breadth-first order is not. Parent dimensions must match the
box's coarse dimensions. Siblings cannot share parent cells; face, edge and
corner adjacency are allowed. Refinement ratios are common at each depth and
may differ between depths. Physical domain boundaries are valid geometry.
There is no compile-time level or child-count limit.

The constructor checks bounds, parent relationships, overlap, level ratios and
integer extent overflow before replacing the destination. Each field stores a
four-dimensional `(component,x,y,z)` array. Components are generic conserved
quantities; this module does not infer thermodynamic temperatures or chemistry.

## Operations

- `initialize_amr_patch_tree_3d`: construct a topology transactionally.
- `prolong_amr_patch_tree_3d`: initialize all descendants with existing PCM
  prolongation from their parents; commit only when every node succeeds.
- `average_down_amr_patch_tree_3d`: traverse nodes in reverse order, reusing
  existing restriction so that descendants are synchronized before parents.
- `composite_integrals_amr_patch_tree_3d`: exclude cells covered by immediate
  children at every node, use the physical volume along its refinement chain,
  and accumulate all remaining cells with compensated summation. Covered
  parent values need not already be synchronized for this inventory.

Invalid topology/state leaves destination objects or inventory vectors unchanged.
The hierarchy owns its arrays; restriction uses a private candidate copy.

## Remaining integration

This foundation does not yet evolve a multi-patch flow. Same-level and
coarse/fine ghost filling, conservative flux correction and subcycling,
thermodynamic recovery, chemistry/transport coupling, dynamic remapping,
MPI decomposition, checkpoint/restart and application input/output must be
connected before claiming a general multi-level solver. EB geometry is separate.
The existing serial/MPI two-level application remains unchanged.
