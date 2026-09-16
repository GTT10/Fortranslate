# Conservative regular 3D patch-tree remapping — 2026-09-17

Source parent: `dcd54f172a1b374e54d0cefed7b4e72c22ede3d2` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Remap module SHA256:
`10ab5ab18044f240ed33206328c5b20c125d3345df7e83dc78c28263b38624d2`.
Test SHA256:
`d82104fe2fb1d08f196ac1e30b95a6bf85155b23f62c96ee4cc4942747c5dee5`.

## Scope

A generic component remap computes volume intersections of exposed old/new
cells, followed by new-tree restriction. It supports changed locations, counts,
depths and refinement ratios with unchanged root dimensions/domain. Covered old
parents do not contribute. Integer least-common-multiple coordinates prevent
floating-point face-position ambiguity, with overflow guards. Component and
coverage sums are compensated. The reactive wrapper uses remapped temperature
only as an EOS starting estimate and publishes state/temperature together.

This is PCM conservative projection, not slope reconstruction. The caller owns
topology publication and regrid scheduling. Exposed-cell descriptor storage
avoids a dense global finest grid, but the overlap search is quadratic; no
large-scale performance claim, MPI remapping, application integration or
physical validation follows from this check.

## Completed focused checks

Reused GNU Fortran 13.3 Debug/OpenMPI 4.1.6 build, bounds and invalid/zero/overflow
traps, MPI/spray enabled. No toolchain change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_remap_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree_remap_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final build and **1/1 focused test passed; 0 failures/skips, 0.04 s**.
The initial test also passed; the final run added an explicit LCM-overflow case.
No tolerance or golden hash changed.

The test independently overlays old tree cells on a 24-cubed dense raster,
then averages its cells into every new cell. The old/new three-level trees
change first refinement ratio from 2 to 3, move/add branches, and refine and
coarsen different portions. Checks include all-cell oracle agreement,
composite component conservation, complete coarsening, exact identity remap
against restricted input, invalid-topology rollback, and checked common-grid
overflow using deep binary/ternary trees with small local arrays.

A nonuniform elementary-H2/O2 mixture additionally checks conserved component
inventory, species closure, positive recovered temperature, EOS correction of
the volume-averaged temperature guess, and preservation of both outputs after a
nonphysical remapped candidate. Provenance/project contract and Python AST of
both checker scripts passed, as did whitespace checks.
Logs: `build/amr-tree-remap-checks/build.log` and `ctest.log`.
