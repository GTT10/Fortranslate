# Multilevel EB tree transfer, inventories and chemistry — 2026-09-18

Implementation: `8cc66cb7a50f232fa7cdec2b9ed57568a1b65344`, based on
`e2a8f77`. New core module: `eb_reactive_patch_tree_3d_mod`.

## Implemented library operations

The implementation reuses `amr_patch_tree_3d` topology and its parent-before-child
ordering, sibling nonoverlap and per-depth refinement-ratio checks. Each node
stores EB geometry, conserved state and temperature using the existing patch
type. The geometry metric comparison already used by two-level updates is now
exported as `matching_eb_metrics_3d`; its comparison rules are unchanged.

Initialization takes one full-domain finest geometry lattice and the root
state. Successive geometric coarsening constructs each depth, then integer
node origins select local geometry. Root EOS validation and recursive PCM
prolongation initialize fluid state; covered storage inherits parent payload.
Every child must coarsen to matching parent-section metrics. Initialization
publishes the tree only after all nodes succeed.

Bottom-up restriction processes descendants before parents and publishes all
nodes transactionally. Composite integration sums only exposed fluid cells:
solid cells and parent regions owned by immediate children are excluded.
Geometry/state shape failures and nonfinite exposed conserved values preserve
the caller's inventory vector. Integration is an inventory operation, not a
physical-state/EOS acceptance check.

The reaction API applies chemistry to exposed fluid cells only at every node,
then averages down. It returns the composite chemical inventory change only
after complete success. Failures on a later node retain the original entire
tree and chemical ledger. No hydro, transport, halo fill, reflux, tree
checkpoint or multilevel application stepping is added in this increment.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_patch_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_patch_tree_3d$' --output-on-failure --output-log build/eb-patch-tree-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 1.28 s (total 1.31 s), zero failures/skips. Six configurations
combine x/y/z planar walls at 0.31 with first-depth ratios 2/3 and second-depth
ratio 2. The root is 4-cubed, with two disjoint children and one grandchild.

Checks prove PCM composite conservation within
`3e-14 * max(1, abs(initial amount))`; a 25% grandchild state change contributes
exactly once; bottom-up restriction leaves exposed composite inventory exactly
unchanged and propagates deep inventories to the root. Deliberately invalid
hidden parent values neither enter the composite sum nor own reaction solves.
Explicit H2/O2 chemistry at 1e-9 s produces nonzero species exchange, retains
mass/momentum/energy to the same roundoff bound, and reports the exact composite
inventory difference. Invalid grandchild reaction/restriction, mismatched
geometry and invalid initialization exercise output and whole-tree rollback.

Provenance, project contract and whitespace checks passed; no Python files
changed. Logs: `build/eb-patch-tree-checks/build.log` and `ctest.log`.
This core-library check does not establish arbitrary-depth hydro/transport,
selected-runtime packaging, MPI integration, external PeleC parity or physical
validation. Those remain completion work.
