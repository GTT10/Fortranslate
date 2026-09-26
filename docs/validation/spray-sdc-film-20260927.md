# Stationary-film thermal evolution in MPI EB SDC — 2026-09-27

Implementation `cab1d10a`, based on `92a3569d` on local `main`.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Scope and algorithm

`advance_mpi_eb_spray_sdc` accepts optional `films`, `film_mass` and
`film_settings` as an all-or-none group agreed across ranks. `film_steps`
returns the sum of old and final-new accepted film substeps, or zero on failure.

Each old/new particle source evaluation is followed by a half-interval film
thermal/evaporation update using its private gas state. The accepted particle
and film integrated sources are added with overflow checks and normalized by
the half interval. Wall heat stays in the wall inventory. The old evaluation
saves a film midpoint; every new evaluation starts from that same midpoint.
Only the final new film state, steps and heat are retained. Final publication
is atomic with gas, reaction history, particles, IDs and all inventories.

This advances existing stationary films only. Particle wall contacts remain
specular and do not deposit new film or generate splash. Film impact coupling,
application selection and recursive particle/film SDC remain open. The
application's SDC wall-impact rejection is unchanged. This is not a claim of
pinned upstream film scheduling parity or temporal order.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed H2O2 tests-enabled build
`build/spray-amr-reactive-debug`; `/usr/bin/mpifort` and `/usr/bin/mpiexec`,
OpenMPI 4.1.6. Existing compiler comparison and internal-callback executable-stack
warnings remain. No build/test failures occurred.

CTest **2/2 passed**, 9.57 s total:

- `spray_mpi_eb_sdc_2`: 5.97 s, no-film regression including final splitting.
- `spray_mpi_eb_sdc_film_2`: 3.56 s, two ranks/two AMR levels, chemistry,
  molecular transport, mixed-liquid particle heating/evaporation/drag,
  EB reflection and two existing mixed-liquid films on different leaf owners.
  Films heat and evaporate with at least two accepted half-interval substeps.
  Combined gas/particle/film/wall/chemical relative inventory residual is
  **8.8274766512609944e-16**, within the unchanged 3e-12 bound.

Rank-disagreeing film presence is rejected collectively. Invalid hydro after
old source work preserves gas, temperature, particle state, film masses,
film temperatures, film population and budget outputs. Successful migration
places each film on its finest owner. External forcing is exactly zero.
Provenance, Python syntax, project contract and diff whitespace passed.
No full matrix, selected-mechanism rerun, installed audit, order study or physical
spray/film validation was performed.

## Reproduction

Artifacts: `build/spray-sdc-film-checks/` (`provenance.log`, `build.log`,
`ctest.log`, `detail.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film_reactive test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc(_film)?_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
