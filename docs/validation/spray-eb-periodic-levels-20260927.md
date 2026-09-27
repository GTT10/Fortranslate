# Periodic EB particle migration and level copies — 2026-09-27

Implementation `8038ba925c62d05124cf140fe7b12ead98bfd63d`, based on `05c451f5`,
local main, version 0.246.0.

## Change

`migrate_eb_tree_spray_parcels` accepts optional periodic axes. The general AMR
migration performs periodic ownership lookup; the EB wrapper then folds private
coordinates into the physical domain before checking fluid occupancy. Every
requested periodic axis must have exactly matching opposite nodal level-set
faces. Existing replicated topology, geometry, periodic-control, payload and
receive-capacity checks remain active. Publication remains transactional.

Physical level selection, virtual copies, grown ghost routing and physical
collection now forward those axes instead of rejecting all periodic EB input.
Ghosts retain canonical coordinates, while the existing periodic grown-box
lookup selects their destination across a seam. This uses the existing
single-destination copy convention; it does not introduce multiple images of
one particle within a destination population.

Recursive particle stage/source integration and film-copy periodic support still
remain to be connected. The application continues to reject recursive periodic
input in this increment.

## Independently completed checks

Artifacts: `build/spray-eb-periodic-level-checks/`. Existing fixed-mechanism Ninja
Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final source: **2/2 passed**, 0.48/0.50 s on two/four ranks, 1.03 s total.
Provenance, build, Python syntax, project contract and whitespace checks passed.
Final logs are `build-final.log`, `ctest-final.log` and `LastTest.log`.

The initial build rejected an overlong call line; it was split. The first test
run failed on both rank layouts because the general migration wraps ownership
lookup while preserving original coordinates. Adding private coordinate wrapping
before EB intersection resolved it. Those failures remain in `build.log` and
`ctest.log`. No skips; existing compiler warnings remain.

The fixture verifies multiple positive/negative wraps, unchanged IDs, velocity
and component inventories, finest physical ownership, coarse virtual ownership,
physical collection and grown ghost routing. A coarse particle at the right
periodic face supplies the left fine patch's ghost box. Nonmatching nodal seams
and rank-disagreed periodic axes reject without publishing wrapped coordinates.
Existing nonperiodic ownership, solid-point and capacity rollback checks also
pass in the same two tests.

No coupled periodic gas/particle step, film image evolution, restart, application
runtime, selected mechanism or physical validation was run.
