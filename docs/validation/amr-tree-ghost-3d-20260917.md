# Regular 3D patch-tree snapshot ghosts — 2026-09-17

Parent: `3117fc0`, local `main`. API and remaining scope:
[patch-tree design](../design_decisions/amr-patch-tree-3d.md).

Added a transactional two-cell halo builder. Checked global integer coordinates
identify peers at the same depth even when they have different parents. Missing
fine coverage falls back through ancestors with PCM. Periodic wrap, outflow,
reservoir inflow, slip/no-slip/isothermal momentum reflection, optional viscous
reflection and moving-wall energy correction follow regular-grid conventions.

This is snapshot state exchange, not time evolution. It does not yet add PLM or
temporal interpolation, EOS temperature recovery, isothermal heat flux,
hydro/reflux/subcycling, MPI exchange or EB coupling. Peer lookup is a node scan;
no performance/scaling or physical-validation claim is made.

## Reproduction

Reused GNU Fortran 13.3 Debug build with bounds/FPE checks:
`build/spray-mpi-migration-gnu-debug`. No compiler, MPI or backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_ghost_3d test_amr_patch_tree_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree(_ghost)?_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests: **2/2 passed, 0 failures, 0 skips**, 0.05 s.
The new test constructs three levels/five boxes, with adjacent grandchildren
belonging to different parents and holes requiring ancestor fallback. An
independent dense-grid oracle prolongs the entire grid at each depth, overlays
fine boxes, then calls the existing regular-grid boundary helper. Every interior
and halo component is exactly equal for every node under periodic, mixed physical
Euler and mixed physical viscous/moving-wall modes. Coarse and fine values are
intentionally different to expose wrong donor selection.

Missing inflow data, unpaired periodic faces, invalid node, nonfinite ancestor
state and 64-level global-coordinate overflow preserve the previous output.
The prior tree conservation test also passed. The final build/test includes a
reordered coordinate expression that avoids intermediate integer overflow.
Build, mechanism provenance, project contract, Python AST syntax for two
verification scripts and whitespace checks passed. No full matrix, installed
binary audit, external solver parity or experimental comparison was run.

Logs and hashes: `build/amr-tree-ghost-3d-checks/`.

```text
8d85c4c7707c4595b0b1f345a6d2fe09128e8ff3d97cf7167dbda9c760a7c08a  src/amr/amr_patch_tree_ghost_3d_mod.F90
161098b6b40472a5abed34a4f2337bc0a966f21071ae64b6eae7a55be9a09aef  tests/unit/test_amr_patch_tree_ghost_3d.F90
```
