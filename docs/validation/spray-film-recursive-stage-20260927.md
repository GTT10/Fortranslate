# Recursive EB film clocks and cleanup — 2026-09-27

Implementation `a2f32902c49bc92f29f7d6233029eaf141ab6890`, based on `120c2b41`,
local main, version 0.246.0.

## Change

`advance_mpi_eb_recursive_film_stage` connects film copy construction and the
existing old/new film stages to level clocks. The parent must finish its
provisional interval before a child starts; child intervals must remain inside
that interval. Current and descendant clocks must agree at an old-stage start.
MPI ranks agree on time and stage controls before proceeding.

Old stages build private virtual and outgoing ghost films. Old/final stages
fill gas and temperature halos, including ancestor temporal endpoints and inflow
metadata. Intermediate iterations return cached source rates without rebuilding
the halos. Populations, clocks, sources and wall inventory publish together only
after a successful stage. The full hierarchy must still remain private to the
caller's outer gas transaction.

`finish_mpi_eb_film_level` verifies that the selected interval and descendants
have finished. It clears virtual films after each level interval and incoming
ghost films only after the final child cycle. Stationary films require no drift
redistribution. This helper is intended for a post-gas-synchronization callback.

## Independently completed checks

Artifacts: `build/spray-film-recursive-stage-checks/`. Existing fixed-mechanism
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

Final source: **2/2 passed**, 0.88 s with two ranks and 0.91 s with four ranks,
1.83 s total (`ctest-final.log`, preserved `LastTest.log`). Provenance, build,
Python syntax, project contract and whitespace checks passed. No failures or
skips in this increment; existing compiler warnings remain. An earlier test run
also passed before additional owner-array and cleanup-clock validation; final
build output is in `build-final.log`.

The two-level fixture checks root copy creation, rejection of a premature child
stage and premature parent cleanup, two child intervals with old/intermediate/
final stages, ghost retention after the first child cycle and release after the
second, unchanged parent physical film payloads during child evolution, and
rejection of a child interval beyond its parent endpoint. Root cleanup succeeds
after both child intervals. Existing particle and film kernel/stage checks run
in these same tests.

These are library sequencing checks with supplied gas states. No coupled
recursive gas/film application interval, collision-created film or splash,
selected-mechanism run, restart or physical validation was performed. Gas source
and predictor callback integration, final physical collection, outer publication
and application routing remain open; recursive film input is not enabled here.
