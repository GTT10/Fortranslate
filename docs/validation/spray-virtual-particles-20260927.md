# Coarse virtual particle creation — 2026-09-27

Implementation `eeb860ff`, based on `51ba17f2`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Population lifecycle primitive

`create_mpi_amr_virtual_particles` selects active physical parcels with a
finest covering level deeper than `target_level`. It copies them without
aggregation and uses level-bounded MPI migration to assign coarse owners.
Both parcel and component outputs publish only after success; input arrays
are read-only. Existing virtual outputs survive invalid inputs, collective
configuration disagreement, migration capacity failures and EB rejection.
Optional `nodal_phi` enables nonperiodic EB surface checks. The finest target
level produces empty outputs.

The caller must supply each active physical particle once, with sorted local
IDs. Source IDs are retained in the separate output for correspondence; virtuals
must not be merged into physical populations or counted in physical inventories.
The caller must regenerate/clear these temporary copies at level boundaries.

The pinned `Source/Particle.cpp::setupVirtualParticles` recursively brings
finer real/virtual populations to coarse levels. AMReX's pinned
`AMReX_ParticleContainerI.H::CreateVirtualParticles` supports a `None`
aggregation branch that copies particles. This implementation provides that
nonaggregated coarse representation directly from all finer physical parcels.
It does not implement cell aggregation or AMReX reserved virtual IDs, nor does
this test establish reference-solver trajectory parity. Ghost populations,
level-specific interpolation/deposition and recursive time integration remain
unconnected; this primitive alone is not a recursive spray solver.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and
`/usr/bin/mpiexec`.

Four CTests passed, 2.03 s total:
`spray_mpi_amr_tree_2`, `_4`, `spray_mpi_eb_tree_2`, `_4`.
The regular three-level fixture checks descendant-only selection at level one,
all descendants at the root, an empty finest-level population, coarse owner
counts, source IDs/coordinates/components, source preservation, capacity
rollback and rank-disagreeing target rollback. EB checks coarse ownership,
unchanged active parcels and retained outputs after a solid-point rejection.
Existing ownership/migration assertions also run.

Provenance, target builds, Python syntax, project contract and whitespace checks
passed. No build/test failures occurred. Existing compiler warnings remain,
including real comparisons and an unrelated application loop-subscript warning.
No full matrix, selected-mechanism build, installed audit, recursive application
run or physical validation was performed.

## Reproduction

Artifacts: `build/spray-virtual-particle-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree test_spray_mpi_eb_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr|eb)_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
