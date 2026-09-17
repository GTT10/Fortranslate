# Multilevel EB temperature tagging — 2026-09-18

Implementation: `078896158d8b520cef0cb4f19e092835dcafeb25`, based on
`adbcdba`. API: `tag_reactive_eb_tree_temperature_3d`.

Exposed fluid cells contribute temperature and exact integer cell extents to
the existing regular AMR deterministic tile builder. Solid cells and
child-covered parent cells do not vote or require valid temperatures; exposed
fluid temperatures must be finite and positive. Fine peaks are not averaged
away. Nondecreasing per-level thresholds select requested refinement levels.
The builder uses breadth-first parent order and k/j/i tiles, then the EB
wrapper removes wholly solid nodes and renumbers descendants.

The caller supplies immutable full geometry on the deeper of the old and
requested lattices. Shared refinement ratios must match. Coarsened geometry
is checked against every old node. Tile dimensions, clipped buffer width,
thresholds and capacity are validated; output topology only publishes on
success. Capacity currently bounds the intermediate regular tile tree before
solid pruning, so a buffered selection may reject even if its pruned tree
would fit. It never silently discards a fluid refinement to meet capacity.
No periodic wrapping, moving geometry or ratio-changing remap is added.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_tagging_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_tagging_3d$' --output-on-failure --output-log build/eb-tree-tagging-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.22 s (total 0.24 s). A five-node
three-level ratio-2 tree with z=0.31 wall contains one exposed grandchild
cell at 1500 K; other fluid cells are 300 K. Solid and hidden parent
temperatures are NaN. Thresholds 1000/1400 select exactly two child levels
and the expected spatial tiles; 1000/1600 select only one level; 2000/3000
coarsen to the root. A maximum-integer buffer is clipped safely and the
result initializes with no wholly solid nodes. Capacity one and invalid hot
cell temperature reject while retaining the previously selected topology.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. Logs:
`build/eb-tree-tagging-checks/{build,ctest}.log`. This is a tagging library
check; application adaptation/restart integration, scaling, convergence and
physical validation remain separate work.
