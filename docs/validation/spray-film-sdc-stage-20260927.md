# Old/new SDC film stages — 2026-09-27

Implementation `5359fc74d2b1d2c43e074f2fd16ad2b530eff483`, based on `7780749c`,
local main, version 0.246.0.

## Change

`spray_level_film_population` stores films, component masses, logical nodes and
fluid-side probes together. `advance_mpi_eb_film_sdc_stage` advances an ordered
physical/virtual/ghost group using caller-supplied stage-time gas halos:

- stage 0 applies one old half step privately and caches midpoint populations,
  old source rates, retained halo rates, optional physical rates and wall heat;
- intermediate iterations return the cached old rates without another exchange;
- the final iteration applies one new half step and publishes all populations
  together. The physical wall ledger contains both accepted half steps.

All three populations contribute to the provisional same-level and retained
halo rates; only physical films contribute optional composite export and wall
heat. Depleted films, masses, nodes and probes compact together before the next
half step. Stage ordering, time, level and optional-export mode are validated
collectively. Arithmetic guards precede source addition, conversion to rates and
wall-ledger accumulation. A late copy failure retains the caller's populations,
source fields, wall ledger and stage cursor.

This is a stage integration primitive. It does not yet create the film hierarchy
or copies, connect collision-created films, invoke the recursive gas callback,
synchronize physical-film feedback at the outer step, or enable recursive films
in the application. Those are still required; the project is not complete.

## Independently completed checks

Artifacts: `build/spray-film-sdc-stage-checks/`. Existing fixed-mechanism Ninja
Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0, `/usr/bin/mpifort`
and `/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, 0.80 s on 2 ranks and 0.82 s on 4 ranks, 1.67 s total.
No build/test failures or skips. Provenance, build, Python syntax, project
contract and whitespace checks passed. Existing compiler warnings remain.

The added planar-EB fixture checks the following with identical films in all
three populations:

- combined rates include all copies while the physical export and wall heat
  count the physical population only;
- stage 0 leaves caller film populations unchanged;
- midpoint films match one direct half step;
- an invalid ghost-film material is ignored by the intermediate iteration,
  which returns cached rates, but rejects the final iteration after earlier
  populations have advanced privately;
- that late failure retains all population/source payloads and permits retry;
- accepted final films match two direct half steps, with physical wall heat
  counted once per half step;
- complete old-stage depletion removes all film metadata, and intermediate/new
  stages succeed with allocated empty populations.

The existing film halo and particle stage regressions run in the same tests.
No full recursive gas/film interval, application run, impact/splash case,
three-level film evolution, selected mechanism or physical validation was run.
