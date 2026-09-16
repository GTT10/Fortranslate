# Coarse-time interpolation for 3D patch-tree ghosts — 2026-09-17

Parent: `e7b509e`, local `main`.

The tree ghost API now accepts paired coarse future fields and per-level time
fractions. Only donors below the target resolution are interpolated; the actual
donor level chooses the fraction even when fallback skips the immediate parent.
Current peers and interiors remain exact. Endpoint fractions use direct copies.
Physical boundary transforms follow interpolation, with reservoir precedence
unchanged. Future peer/finer fields are not required. Invalid paired arguments,
shapes, finite values or fractions leave the destination unchanged.

This prepares the ghost API for subcycling; it does not schedule steps or evolve
flow. Callers supply consistent time brackets and fractions. PLM interpolation,
hydro/reflux, thermodynamic recovery and application integration remain open.
See the [API design](../design_decisions/amr-patch-tree-3d.md).

## Reproduction and evidence

Reused GNU Fortran 13.3 Debug build with bounds/FPE checks,
`build/spray-mpi-migration-gnu-debug`. No compiler/MPI/backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_ghost_3d test_amr_patch_tree_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree(_ghost)?_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests: **2/2 passed, 0 failures, 0 skips**, 0.07 s.
The expanded halo test compares all values with an independent whole-domain
PCM/overlay oracle at depths 1 and 2, for fractions 0, 1 and `[0.25,0.75]`.
Periodic and mixed physical viscous/moving-wall cases cover both parent-level
and root fallback, different-parent fine peers, and exact fine interiors.
Future fine fields remain unallocated. Missing paired arguments, wrong fraction
count, extrapolation, NaN fraction and an absent coarse peer endpoint reject
without altering output. Existing snapshot/overflow and tree conservation
checks also pass within these two tests.

Build, provenance preflight, project contract, Python AST syntax for two
verification scripts and whitespace checks passed. No full matrix, installed
binary audit, external solver parity, scaling or physical validation was run.
Logs and hashes: `build/amr-tree-temporal-ghost-checks/`.

```text
bfdbb95219d81918a1a90e3f224b9d3d1dad5d70e4f0412da7920ed6c1e9a033  src/amr/amr_patch_tree_ghost_3d_mod.F90
1a3f77ee168c76f6cff637b06595bf1dd2883e95c7e24103b485cc34c7779684  tests/unit/test_amr_patch_tree_ghost_3d.F90
```
