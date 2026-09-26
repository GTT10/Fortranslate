# Recursive particle copy refresh — 2026-09-27

Implementation `b6c95a2b`, based on `7f464038`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`refresh_mpi_amr_level_copies` consumes physical population snapshots indexed
by logical level. Before the selected level moves, it builds nonaggregated
virtual copies from all logical descendants and next-level ghost copies from
that level's own physical population. Selecting by logical ownership preserves
particles that have crossed a geometric coarse/fine boundary but have not yet
been synchronized into a different physical level.

Descendant collection validates live IDs and removes dead parcels; retained
routing places its virtual copies on the requested parent level. The next-level
ghost route filters against grown child boxes without reclassifying the source
population by geometric level. The finest level returns allocated empty virtual
and next-ghost populations. The source snapshots remain read-only; both outputs
publish only after all operations succeed. MPI agrees level, grow, capacity
and EB presence before entering the level-dependent operations.

The caller must supply descendants at the selected interval's start time,
retain incoming ghosts across child substeps and invoke refresh at the correct
parent boundary. This routine does not inspect clocks or perform physical
synchronization, virtual aggregation, recursive SDC dispatch or source prediction.
Those connections remain open.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_amr_tree_2` and `_4` passed:
2/2, 0.49/0.50 s; CTest total 1.04 s. No build/test failures or skips occurred.
Provenance, Python syntax, project contract and whitespace checks passed.

The three-level periodic-x fixture moves logical fine particles outside the
fine valid box into parent geometry and logical parent particles into the fine
valid box. Refresh still selects the former as parent virtuals and the latter
as fine ghosts, with expected IDs, counts, owners, nodes, positions and component
inventories. A capacity rejection preserves both previously published copies;
source physical populations remain unchanged. Refresh at the finest level
returns empty outputs. Existing initialization, retained routing, collection,
virtual/ghost and migration checks also run.

The new optional EB route compiles but was not independently run here. This is
not full recursive evolution, external solver numerical parity or physical
validation. No full matrix or installed audit was performed.

## Reproduction

Artifacts: `build/spray-copy-refresh-checks/`, including build, provenance,
contract, CTest summary and copied detailed `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
