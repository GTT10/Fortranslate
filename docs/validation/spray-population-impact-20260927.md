# Physical impact films across particle SDC stages — 2026-09-27

Implementation `73c63b6b15e31c8152063238724e9ef7fc592f58`, based on `ac728d0c`,
local main, version 0.246.0.

## Change and reference

`advance_mpi_eb_population_stage` accepts the complete optional film/settings
group and routes it only to the physical population during a positive drift.
Virtual and ghost particles retain reflection without creating impact films,
splash children or IDs. The private population-halo wrapper keeps film state
private through particle redistribution, so a failure after impact flight also
rolls back films and the watermark. The three-population wrapper publishes only
after all populations and source accumulation succeed.

This physical-only splash condition follows pinned PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/SprayParticles.cpp`, lines 164–165: active particles exclude virtual
and ghost populations, and splash requires active particles and movement.
The source was inspected with `git show` in
`/mnt/sn850/Projects/gtt2/external/PeleC/Submodules/PelePhysics`.

`advance_mpi_eb_population_sdc_stage` additionally caches the full post-impact
film snapshot in `spray_sdc_population_stage%midpoint_films` and
`%midpoint_film_mass`. The old stage reserves IDs and stores the private snapshot;
intermediate iterations return cached sources and do no impact work. The final
half kick has zero drift and publishes the cached film inventory together with
particles, without duplicate deposits or wall inventory. Film-option presence
must remain consistent across ranks and stages.

This does not yet connect newborn film placement or thermal evolution to the
recursive gas/film coordinator. That integration must reconcile impact snapshots
with film old/new clocks rather than overwrite an independently evolved film
state. The recursive application remains disabled for impacts.

## Independently completed checks

Artifacts: `build/spray-population-impact-checks/`. Existing fixed-mechanism
Ninja Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and
OpenMPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final source: **2/2 passed**, 0.94 s on two ranks and 0.98 s on four ranks,
1.97 s total. Provenance, build, Python syntax, project contract and whitespace
checks passed. The first build failed on an overlong extended signature; splitting
the declaration fixed it (`build.log`). Final build output is `build-sdc.log`,
and final test output is `ctest-final.log` with preserved `LastTest.log`.
No test failures or skips; existing compiler warnings remain.

For the preceding 1 m/s deposition and 100 m/s splash fixtures, population-stage
physical film component masses, IDs and normal wall impulse exactly match the
direct halo kernel. Copies remain single reflected parcels. Zero receive capacity
forces a population redistribution failure after private impacts; films, component
inventories, ledgers and IDs remain unchanged and a normal-capacity retry passes.
The SDC sequence then verifies deferred publication in old/intermediate stages,
retained film snapshots, unchanged IDs after the old stage and exact final film
inventory without duplicate wall impulse.

Existing kernel inventory and source-transfer rollback checks also pass. No
coupled recursive impact application, film thermal ordering, selected mechanism,
restart, physical validation or full external-solver parity was run.
