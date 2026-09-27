# Recursive gas/particle/stationary-film coupling — 2026-09-27

Implementation `e4e6d4234211fedd7fefd496041d2839f18edacf`, based on `a94dfdd4`,
local main, version 0.246.0.

## Change

`advance_mpi_eb_spray_sdc_subcycle` accepts the optional `films`, `film_mass`
and `film_settings` group plus `film_steps`. MPI ranks must supply the same
complete group. The wrapper builds a private physical film hierarchy with empty
virtual/ghost populations and connects the recursive film stage and post-level
cleanup helpers to the gas driver.

The gas source callback adds particle and film rates. Predictor construction
adds both old rates and retained halo remainders, using combined parent endpoint
rates. Final synchronization subtracts the combined provisional leaf source and
adds only the actual physical particle and film exchanges. Physical film wall
heat joins the wall inventory once per accepted level interval; virtual/ghost
film exchanges do not enter that inventory.

After recursive synchronization, only physical films are flattened and validated
through the existing migration operation. Gas, temperature, particles, films,
component masses, IDs and ledgers publish only after all stages and final checks
succeed. `film_steps` counts the two film stages per accepted level interval,
including empty physical populations; it is zero on failure.

This connects existing stationary films to coupled recursive evolution. Particle
collisions still do not create films or splash particles in this path. The
application's recursive wall-film rejection remains until that lifecycle is
connected. Source synchronization remains PeleF's conservation extension, not
proof of exact pinned PeleC numerical parity.

## Independently completed checks

Artifacts: `build/spray-recursive-film-checks/`. Existing fixed-mechanism Ninja
Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive(_film)?_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, 24.55 s total. The existing particle-only test took 17.11 s;
the new film test took 7.40 s. Both use two MPI ranks. Provenance, build, Python
syntax, project contract and whitespace checks passed. No failed builds/tests or
skips; existing compiler warnings remain. Test output is preserved in
`ctest.log` and `LastTest.log`; final build output is `build-test.log`.

The new fixture places one mixed film on each physical level and runs a reacting,
transport-enabled gas interval with an evaporating particle. Both films heat and
evaporate, retain their physical owner counts, and complete six film stages
across one root and two child intervals. Combined gas/particle/film mass, species,
momentum and energy satisfy the unchanged `3e-12` normalized inventory criterion;
the maximum residual is **1.2621683035207835e-15**. The particle-only residual is
**1.1112755134928231e-15**.

An invalid gas control and then a film material conductivity range that excludes
the film temperature each reject the entire interval. The second failure occurs
after the particle stage has been attempted. The test verifies unchanged gas,
temperature, particle/film component inventories, film temperature and mass,
IDs, clocks reported to the caller and ledgers, followed by a successful retry.

No four-rank, selected-mechanism, application, restart, impact/splash, temporal
accuracy or physical-validation run was performed for this increment.
