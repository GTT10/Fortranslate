# Fine-level ghost particle creation — 2026-09-27

Implementation `462705b5`, based on `21e0004c`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`locate_amr_tree_ghost_cells` looks up a requested fine level's boxes enlarged
by `grow` fine cells. It returns node zero when no box matches and local cell
indices that can be outside the valid patch. It tries valid boxes first,
then enlarged boxes in tree order. Periodic images participate in lookup
without changing the particle's coordinates.

`create_mpi_amr_ghost_particles` first validates/migrates a private physical
population. It selects parcels with finest covering level `target_level-1`
that intersect the grown target boxes, partitions copies by the selected
fine owner and exchanges them. The returned target-node vector is rebuilt
and checked after receive sorting. Parcel/component/node outputs publish
together after success; source arrays are read-only. Optional nodal geometry
uses EB migration and solid-point checks. Target level, halo width and EB
presence agree collectively. Capacity exhaustion or other rejection retains
all previous outputs.

The pinned `Source/Particle.cpp::setupGhostParticles` calls AMReX
`CreateGhostParticles(level,ngrow)` and adds copies at `level+1`.
This implementation provides the immediate-coarse selection and fine-owner
copy construction. It retains source IDs in a separate population; ghosts
must not enter physical inventories or checkpoints. Tree-order tie selection
has not been compared with AMReX locator ordering. Ghost gas fields,
level-local evolution/deposition, source transfer, lifecycle refresh/removal
and recursive application routing remain open. Calling the existing composite
coupler on these copies would double-count particles and is not supported.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and
`/usr/bin/mpiexec`.

Four CTests passed, 2.05 s total:
`spray_mpi_amr_tree_2`, `_4`, `spray_mpi_eb_tree_2`, `_4`.
The regular three-level fixture selects particles immediately outside a fine
patch and across a periodic boundary into a one-cell halo. It checks fine
owner counts, node IDs, unchanged coordinates/components, capacity rollback,
rank-disagreeing width rollback, empty zero-width output and exclusion of
particles already at the requested level. EB checks fine-owner halo creation
and retained output after a solid-point rejection. Existing virtual population
and migration assertions also run. Overlapping-grown-box tie cases were not
independently exercised.

Provenance, target builds, Python syntax, project contract and whitespace checks
passed. No build/test failures occurred. Existing compiler warnings remain.
No full matrix, selected-mechanism build, installed audit, recursive application
run, temporal-order check or physical validation was performed.

## Reproduction

Artifacts: `build/spray-ghost-particle-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree test_spray_mpi_eb_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr|eb)_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
