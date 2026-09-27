# Newborn film placement and thermal ordering — 2026-09-27

Implementation `f2cddb59774f663b34f8980471f90329d4306414`, based on `3c95dd39`,
local main, version 0.246.0.

## Change

`migrate_mpi_amr_wall_films` accepts paired `logical_level` and `grow` controls
instead of `maximum_level`. It retains films on the generating level, using the
same valid-box-first, grown-box-second lookup as particles. Deposits beyond every
grown box reject rather than disappear. Root placement uses the root owner;
physical domain bounds and EB wall/probe validation remain enforced. Payloads,
nodes and probes publish together, with collective capacity/control validation.
Omitting the new controls retains finest physical ownership.

`append_mpi_eb_film_births` adds newly created deposits to the physical film
midpoint after its old half kick. It preserves thermally advanced existing films,
old source rates, old wall inventory and virtual/ghost copies. The combined
physical midpoint is redistributed at its generating logical level. Failure
preserves the midpoint; an attachment marker rejects duplicate calls. The next
old film stage resets that marker.

The final film half kick then updates the newborns exactly once. This ordering
follows pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/SprayParticles.cpp`: `CreateSBDroplets` runs after the existing
particle loop. It does not prove full numerical parity with that implementation.
The source was inspected with `git show` in
`/mnt/sn850/Projects/gtt2/external/PeleC/Submodules/PelePhysics`.

These APIs still need to be connected to the recursive particle/film callback:
extract only appended impact deposits, update the private film midpoint, and
avoid publishing a stale particle-side film snapshot over evolved film state.
Application routing remains disabled for recursive impacts.

## Independently completed checks

Artifacts: `build/spray-film-logical-placement-checks/`. Existing fixed-mechanism
Ninja Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and
OpenMPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(film|sample)_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**4/4 passed**, 3.01 s total. Film placement took 0.49/0.50 s on two/four
ranks; sample/stage checks took 0.99/0.98 s. Provenance, build, Python syntax,
project contract and whitespace checks passed. No failed builds/tests or skips;
existing compiler warnings remain. Final outputs are `build-births.log`,
`ctest-births.log` and preserved `LastTest.log`.

The placement fixture verifies exact film payloads and probes on a selected
fine level's grown boxes, root routing, return to finest physical owners,
capacity rollback, missing-grown-box rejection and rank-disagreed controls.
The thermal fixture checks failed attachment rollback, unchanged old rates,
unheated newborn component masses, duplicate attachment rejection and final
newborn masses exactly matching a single half-kick reference. Existing copies
retain their original population counts.

No coupled impact-to-thermal application, recursive-driver handoff, restart,
selected mechanism, physical validation or full external-solver parity was run.
