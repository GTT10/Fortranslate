# MPI EB spray SDC kick/drift — 2026-09-27

Implementation `2e1e1abc`, based on `b1a42c3f` on local `main`.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics pin `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change and reference boundary

Particle coupling now distinguishes the physical kick duration from motion:
`kick_drift` updates drag, temperature, evaporation and breakup for half the
supplied interval and drifts for the full interval; `kick` performs the same
half update without drift or wall impacts. Physical substep limits scale with
the half update; motion CFL applies only when moving. Default `coupled` retains
its existing behavior. Mode disagreement across ranks fails transactionally.

The common-step MPI EB wrapper exports old/new half-kick increments as rates
by dividing by half the full interval. Each new kick is recomputed from one saved
old kick/drift endpoint. Failed gas updates cannot publish parcels, masses,
IDs, gas, temperature or inventories. Accepted wall/outflow increments include
the old stage and only the final new stage. Particle feedback is an internal
exchange and is not added to the external forcing inventory.

Pinned PelePhysics `Source/Spray/SprayParticles.cpp`, `SprayFuelData.H` and
`Drag.H` use half-duration physical updates with full-duration drift and
half-duration source normalization. Pinned PeleC `Source/Particle.cpp` calls
old `particleMKD` at sub-iteration zero and new `particleMK` at the final
sub-iteration. PeleF reevaluates a private new kick on every SDC iteration;
this is a bounded coupled algorithm, not exact upstream scheduler parity.

The new wrapper is library-only, common-step and specular-contact only. It does
not yet wire films, computational multiplicity splitting, particle acceleration,
inlet/isothermal boundary data, recursive particle scheduling, application
selection or spray SDC restart. Gas forcing is supported by the underlying gas
API. The application SDC spray rejection remains in place.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled
`build/spray-amr-reactive-debug`; `/usr/bin/mpifort`, `/usr/bin/mpiexec`,
OpenMPI/OpenRTE 4.1.6. Existing real-comparison warnings and the linker
executable-stack warning for internal Fortran callbacks remain.

- Provenance, builds, Python syntax, project contract and diff whitespace passed.
- Particle tests **2/2 passed**, 1.49 s: `spray_mpi_amr_step_2` and
  `spray_mpi_eb_step_2`. Interior half-kick gas/source/particle updates exactly
  match a legacy half-duration interval. Kick preserves position; kick/drift
  moves by full interval times post-kick velocity. A subsequent kick preserves
  that endpoint. Rank mode disagreement retains state and source.
- Multi-substep EB kick/drift combined inventory relative residual:
  **2.5849394142282115e-16** (bound 3e-12).
- Coupled wrapper test **1/1 passed**, 3.54 s total:
  `spray_mpi_eb_sdc_2`, two ranks, two SDC iterations, two AMR levels,
  chemistry, all transport components, mixed-liquid heating/evaporation/drag
  and EB reflection. Combined gas/parcel/wall/chemical inventory relative
  residual **1.1968665763966046e-15** (bound 3e-12).
- Invalid hydro coefficient fails after the old particle callback, retaining
  gas, temperature, parcels, component masses, ID and all output ledgers.
  Nonzero particle feedback remains separate from zero external forcing.

No build or test failures occurred. No full matrix, installed-binary audit,
selected-mechanism rerun, upstream trajectory comparison, order study or
physical spray/LES validation was performed.

## Commands and artifacts

Artifacts: `build/spray-kick-drift-checks/`, including `build.log`,
`ctest.log`, `detail.log`, `wrapper-build.log`, `wrapper-test-build.log`,
`wrapper-ctest.log`, `wrapper-detail.log`, `provenance-final.log` and
`contract-final.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_step test_spray_mpi_amr_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr_step|eb_step)_2$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
