# Synchronized SDC population initialization — 2026-09-27

Implementation `041a594e`, based on `65f20f7d`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`initialize_mpi_amr_sdc_populations` builds three separate arrays of
`spray_level_population`, indexed by level+1, from one synchronized physical
snapshot. It selects exact-level physical particles, nonaggregated virtual
representations of deeper particles and immediate-coarse ghost copies in grown
fine boxes. All copies precede particle motion. Virtual node assignments use
capped-level lookup; physical and ghost assignments use their existing routing.
Root ghosts and finest virtuals are allocated empty populations.

Level count, grow width, capacity and EB presence are collectively agreed.
All three arrays publish only after every creation and node lookup succeeds.
The original parcel/component inputs remain read-only. The capacity applies to
each rank-local population separately. Existing source validity and routing
contracts remain in force, including the nonperiodic EB restriction.

This initializes an outer synchronized step only. During recursive evolution,
virtuals must be refreshed and ghosts retained/cleared at the appropriate level
boundaries. The helper does not implement that evolution, advance particles,
assemble source halos or connect recursive application dispatch. Only physical
populations may be supplied to synchronized final collection.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_amr_tree_2` and
`spray_mpi_amr_tree_4` passed: 2/2, 0.49/0.51 s; total 1.05 s.
No build/test failures or skips occurred.

The three-level periodic-x fixture verifies exactly one physical copy per ID,
level/owner assignments for all populations, positions/velocities/temperatures
before motion, liquid component payloads, expected root/parent virtual counts
and finest ghost counts. A reduced capacity fits the physical populations but
fails later during root virtual creation; all three previously published
population arrays retain their IDs, nodes and component inventories.
Existing selection, collection, virtual, ghost and migration checks also run.

The new initializer's optional EB path is compiled but was not independently
executed here. No recursive particle evolution, full matrix, installed audit,
external solver numerical parity or physical validation was performed.
Provenance, Python syntax, project contract and whitespace checks passed.

## Reproduction

Artifacts: `build/spray-sdc-population-checks/`, including build, provenance,
contract, CTest summary and copied detailed `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
