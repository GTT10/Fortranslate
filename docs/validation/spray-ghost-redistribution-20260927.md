# Redistribution of evolved ghost particles — 2026-09-27

Implementation `3b792794`, based on `b934db23`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`redistribute_mpi_amr_ghost_particles` shares a private routing implementation
with ghost creation. Creation retains its immediate-coarse physical selection.
Redistribution validates component inventories, removes zero-mass entries,
then uses the existing private physical migration checks and target grown-box
lookup without that selection. Thus a ghost that entered a fine valid cell is
retained. Copies outside all grown boxes are removed. Destination exchange
sorts IDs; a post-exchange lookup reconstructs node indices and checks ownership.
All three outputs publish together.

This operates only on temporary ghost copies. It neither updates physical
population inventories nor creates missing ghosts from physical particles.
The caller schedules regeneration, synchronization and recursive lifecycle.
Optional EB validation retains the existing nonperiodic restriction. Private
validation may move copies through finest owners before the target-level
exchange; it does not publish a physical population.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Four tests passed in 2.02 s total:

- `spray_mpi_eb_tree_2`: 0.49 s.
- `spray_mpi_amr_tree_2`: 0.48 s.
- `spray_mpi_eb_tree_4`: 0.50 s.
- `spray_mpi_amr_tree_4`: 0.49 s.

The added regular three-level fixture changes the fine node's owner, moves
some ghosts inside its valid cells, leaves periodic-image ghosts in its halo,
moves some outside the grown region and marks others dead at the domain upper
boundary. Capacity zero rejects while preserving checked counts, coordinates,
components and node outputs. Successful redistribution retains the expected
IDs/components in sorted order on the new owner and removes dead/escaped copies.
Repeating with zero grow retains only the valid-cell ghosts. Existing regular
and EB creation, virtual-population and migration tests also run.

The new redistribution API's EB path was compiled but not independently called
by the added fixture. Provenance, build, Python syntax, project contract and
whitespace checks passed. No build/test failures or skips occurred. No full
matrix, installed audit, recursive application, external solver parity or
physical validation was performed.

## Reproduction

Artifacts: `build/spray-ghost-redistribution-checks/`, including CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree test_spray_mpi_eb_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr|eb)_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
