# EB impact deposition/splash in common-step SDC — 2026-09-27

Implementation `cd517a3a`, based on `8ed16a23` on local `main`.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

`advance_mpi_eb_spray_sdc` accepts optional `film_impacts` (default false).
True requires the film population, component masses and settings; presence and
activation agree collectively before source callbacks. It forwards the private
film state and shared ID watermark to the existing particle interval impact
handler. The old kick/drift stage creates deposited film and splash children.
New stages are kick-only, so they cannot repeat the contact. Newly deposited
film participates immediately in the existing half-interval thermal update.
`model_rejected` forwards empirical-model rejection without publishing state.
Final computational splitting remains after the accepted SDC stages.

Default specular contacts and existing-film-only thermal evolution are retained.
This is library integration; application film selection and recursive particle/
film SDC remain open. It does not establish upstream scheduling/trajectory parity
or physical validity of the empirical impact model.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, fixed H2O2 tests-enabled build
`build/spray-amr-reactive-debug`, OpenMPI 4.1.6 via `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Existing comparison and executable-stack warnings remain.
Final CTest **2/2 passed**, 7.18 s total:

- `spray_mpi_eb_sdc_film_2`: 3.58 s, default specular/thermal regression.
- `spray_mpi_eb_sdc_impact_2`: 3.55 s, two ranks/two levels, two SDC iterations,
  chemistry and transport, mixed-liquid heating/evaporation/drag and one EB
  splash event. Four children receive IDs 2–5 with parent ID 1 and leave the
  wall into fluid. One newly deposited film migrates to its finest owner and
  heats above its initial 350 K; both film source halves run. Combined
  gas/particle/film/wall/chemical relative inventory residual is
  **1.3456240479579768e-15**, below the unchanged 3e-12 bound.

Invalid hydro after the old source evaluation retains the original particle,
ID watermark, empty film population, gas, temperature and budgets. Film-group
presence disagreement rejects collectively. The test does not exercise splash
outflow or KH interactions under SDC; their preexisting split tests remain separate.

The initial impact fixture omitted allocation of its source-output array and
failed the accepted-step assertion. Allocating the required array fixed the
fixture; no solver tolerance or model guard changed. No build failures occurred.
Provenance, Python syntax, project contract and diff whitespace passed. No full
matrix, selected-mechanism rerun, installed audit or physical validation was run.

## Reproduction

Artifacts: `build/spray-sdc-impact-checks/`, including `ctest-initial.log`,
`ctest.log`, `detail.log`, `build.log`, `build-final.log`, `provenance.log`
and `contract.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_reactive test_spray_mpi_eb_film_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_(film|impact)_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
