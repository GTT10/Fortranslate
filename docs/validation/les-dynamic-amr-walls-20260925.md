# Dynamic LES AMR wall integration — 2026-09-25

Implementation `ad0d9c5`, based on `512bd59`. Serial/MPI AMR timestep selectors and SSPRK2 operators now
pass patch-local physical faces and wall temperatures to the Dynamic LES kernel.
The coupled reacting-spray API and application accept slip, adiabatic no-slip
and isothermal walls, including prescribed tangential wall motion. Coarse/fine
and interior patch faces keep their existing donor fluxes. Nonzero pressure
outlets remain rejected in Dynamic mode.

Moving walls retain the existing gas-input requirement for molecular transport
and viscosity. Thermal/work fluxes enter the existing conservative composite
ledger, restriction/reflux and transactional state publication. No checkpoint
schema or numerical tolerance changed.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled test and
application targets built; selected H2O2 tests-disabled application built.
Selected wall execution was not performed.

- `spray_les_mpi_amr_tree_2`: **passed**, 109.33 seconds. Adds mixed slip,
  adiabatic moving and isothermal moving walls on root-only and nested trees.
  Checks nonzero boundary transfer, composite conservation, timestep/diffusivity
  and state agreement with serial, concentrated/distributed ownership and
  rollback. Existing periodic/inflow/outflow and static cases remain covered.
  Timeout is 180 seconds (previously 120) for the added boundary mode.
- Wall application checker: **passed, exit 0**, using Gaussian Dynamic LES,
  molecular viscosity, moving walls at 2 m/s and isothermal walls at 1100 K.
  Continuous and 2-to-4-rank continuation agree; scheduled regridding, injection
  continuation and changed-context/truncation/overwrite rejection checks pass.

The initial two-test CTest invocation had **one pass and one failure** (109.99
seconds total). The application failed input validation in 0.61 seconds because
its new moving-wall fixture omitted molecular viscosity. The fixture was fixed
to enable transport and viscosity, then its same checker was run directly and
passed. The full CTest pair was not repeated. There were no skips or numerical
test failures. Provenance, Python syntax, project contract and whitespace pass.

This is focused wall integration evidence, not a near-wall turbulence model,
external PeleC parity, physical validation, selected-fuel wall qualification or
a full configuration matrix. EB, pressure-outlet Dynamic LES and LES subcycling
remain outside this change.

## Commands and artifacts

Logs: `build/les-dynamic-amr-wall-checks/` contains `ctest-initial.log`,
`pelef-dynamic-walls-build.log`, `pelef-dynamic-walls-selected.log` and
`pelef-dynamic-walls-application.log`.
Successful application artifacts:
`build/dynamic-walls-application/restart_9py560pb/`.
The standalone serial test target was built; serial operators were exercised
inside the MPI comparison test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_tree test_les_amr_tree pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_amr_dynamic_walls|spray_les_mpi_amr_tree_2)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-debug/pelef_mpi_amr_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/dynamic-walls-application --les --dynamic --dynamic-walls
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
