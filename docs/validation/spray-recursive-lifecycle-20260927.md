# Recursive spray post-level lifecycle — 2026-09-27

Implementation `f6243f74c6c846b1298c1fb52518cc44670fcc4a`, based on `96302623`,
local main, version 0.246.0.

## Reference and implementation

Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/Particle.cpp::postTimeStepParticles` removes virtuals every level step,
removes incoming ghosts on the final parent subcycle, and redistributes physical
levels at/above the current level when `(iteration < ncycle && level < finest)`
or at root (also after injection, which is not connected here).
Pinned AMReX `1546c1ccbf41c704a33790f4bf445a2fcee1372e`,
`AMReX_ParticleContainerI.H::locateParticle` first searches valid boxes and then
falls back to grown boxes at the minimum requested level. Copies of these pinned
sources are retained with the artifacts below; external working HEADs were not
used as the reference.

The gas recursive driver now offers `post_level(depth,time,interval,iteration,
ncycle,comm,ok)` after descendants and gas synchronization. Callback presence is
MPI-agreed. One-rank rejection rolls back the outer gas transaction. Callers must
keep callback-owned state private until the driver succeeds.

The spray wrapper checks subtree clocks, redistributes only synchronized physical
populations at the reference schedule with grow=iteration, empties its virtuals,
and empties incoming ghosts only on the final subcycle. Copies do not participate
in physical redistribution. All populations remain private until the coupled
outer transaction commits.

`redistribute_mpi_amr_subtree_particles` validates unique live subtree IDs, selects
the finest geometric level at/above the minimum and uses grown minimum-level
boxes for particles that left its valid region. It preserves coarser logical
populations and publishes all affected levels together. Escapes beyond the grown
support or capacity failures reject without deleting physical inventory.

## Independently completed checks

Artifacts: `build/spray-recursive-lifecycle-checks/`.
Existing Ninja Debug build `build/spray-amr-reactive-debug`, fixed H2O2 mechanism,
GNU Fortran 13.3.0 (`/usr/bin/gfortran`), OpenMPI/OpenRTE 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`, verified from CMakeCache).
The shell-default `mpirun` belongs to Intel MPI; it was not used by these tests.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_amr_tree_[24]|spray_mpi_eb_sdc_recursive_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, build, syntax, contract and whitespace checks passed. Build warnings
remain, including existing real comparisons, an application loop subscript warning
and the test callback executable-stack linker warning. No build or test failures,
no skips. CTest **3/3 passed**, 18.50 s total:

- Tree redistribution: 2 ranks 0.48 s, 4 ranks 0.50 s. Three-level crossing
  promotes/demotes physical parcels, retains a fine-only grown-box fallback,
  preserves untouched coarse populations and component inventory, checks every
  physical ID exactly once, and rejects capacity overflow atomically.
- Coupled recursive core: 2 ranks 17.47 s. Existing interior fine parcel with
  chemistry and transport retains residual 1.5189353956769058e-15 against the
  unchanged 3e-12 tolerance. Additional gas callback checks verify child 1,
  child 2, root order, exact interval metadata, collective rejection from one
  rank at root, unchanged gas/temperature/ledgers, and accepted retry.

`build.log`, `build-final.log`, `ctest.log` and `LastTest.log` record these checks.
This verifies the new scheduling and subtree routing. It does not yet establish
coupled conservation while a physical parcel crosses a coarse/fine interface,
full three-level gas/particle evolution, virtual aggregation, application
subcycling, injection/regrid/restart integration, or external solver parity.
