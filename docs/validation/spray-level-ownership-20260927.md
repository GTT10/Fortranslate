# Level-bounded particle ownership — 2026-09-27

Implementation `a317a0d4`, based on `26ddff6d`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Provenance preflight H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation and scope

The serial particle locator/partitioner and regular/EB MPI migration APIs
accept optional `maximum_level`. It is zero based and must name an existing
level. The finest covering node at or below that ceiling wins, including
covered parents. Default omission retains the existing finest-level selection.
Periodic coordinate wrapping and half-open cell conventions are unchanged.

MPI validates the effective ceiling collectively before partition/exchange;
post-exchange ownership verification uses the same ceiling. Failure retains
caller parcel/component arrays and serial partition outputs. The EB wrapper
continues to validate the original finest-resolution nodal surface, including
solid-point rejection, even when ownership is assigned to a coarse node.
No parcel attributes or coordinates are changed by this operation.

The pinned reference's `Source/Particle.cpp::particleMKD` creates virtual and
ghost particles for level-specific source construction. Assigning a private
particle copy to coarse owners without altering the gas tree is a prerequisite
for that lifecycle. This change does not create/aggregate those copies or
implement their gas stencils, clocks, source transfer or recursive scheduling.
It is not evidence of completed recursive spray SDC or reference-solver parity.

## Focused results

GNU Fortran 13.3.0, Ninja Debug, OpenMPI/OpenRTE 4.1.6, fixed-mechanism build
`build/spray-amr-reactive-debug`, `/usr/bin/mpifort`, `/usr/bin/mpiexec`.

Five CTests passed, 2.02 s total:

- `spray_amr_tree_partition`: three-level lookup at ceilings zero and one,
  exact local cell indices, periodic points and invalid-level rollback.
- `spray_mpi_amr_tree_2` and `_4`: root/intermediate owner distributions,
  retained IDs/coordinates/components, rank-disagreeing and invalid ceiling
  rollback, plus the existing default migration, capacity, domain/owner
  disagreement, payload, crossing, regrid, duplicate-ID and empty tests.
- `spray_mpi_eb_tree_2` and `_4`: coarse covered-owner migration followed by
  default finest migration, retained IDs/components, plus existing geometry,
  solid-point, capacity, crossing and empty-population checks.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures occurred. Existing real-comparison warnings
remain. An exploratory CTest listing reported unrelated unbuilt targets;
none were selected or counted as executed. Full matrix, selected-mechanism
builds, installed binaries and physical validation were not run.

## Reproduction

Artifacts: `build/spray-level-view-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_amr_tree_partition test_spray_mpi_amr_tree test_spray_mpi_eb_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(amr_tree_partition|mpi_(amr|eb)_tree_[24])$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
