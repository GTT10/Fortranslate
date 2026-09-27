# Periodic common-step film ownership — 2026-09-27

Implementation `31c20d24`, based on `9679053d`, local main, version 0.246.0.

## Change

The common particle coupling and stationary-film thermal update now forward
their periodic axes to film migration before and after exchange. EB survivor
migration also receives the axes. Previously these migration calls used the
default nonperiodic ownership even when the sampling/deposition path was
periodic. Exact wall coordinates and transactional publication remain unchanged.

## Independently completed checks

Artifacts: `build/spray-periodic-common-film-checks/` (`provenance.log`,
`build.log`, `ctest.log`, `LastTest.log`). Existing fixed-mechanism Ninja Debug
build, GNU Fortran 13.3.0, OpenMPI 4.1.6, `/usr/bin/mpifort`, `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film_step pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(periodic_film_step_(2|4)|periodic_impact_application)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, no failures or skips: periodic thermal update on two/four ranks
0.65/0.61 s, existing periodic impact application 53.38 s, total 54.69 s.
Provenance, build, Python syntax, project contract and whitespace passed.
Existing compiler warnings remain.

The new thermal fixture places one mixed film exactly at the upper x endpoint
of a three-level tree. It checks that the endpoint is retained while the film
belongs to the wrapped finest-level owner. Both films heat and evaporate, gas
source reconstruction passes at the existing `3e-13` criterion, and combined
gas/film/wall inventories satisfy the existing `2e-12` criterion. Covered gas,
empty populations and rejection rollback checks also run.

The existing application regression covers periodic deposition and splash,
film state and two-to-four-rank continuation. Cases:
`build/spray-amr-reactive-debug/tests/spray_mpi_eb_periodic_impact_application/impact_nfkrpecs/`.

This uses planar EB geometry. It is not an oblique-seam thermal qualification,
selected-mechanism matrix, physical-validation or external-solver parity result.
