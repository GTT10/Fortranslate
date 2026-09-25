# Coupled MPI EB spray domain walls — 2026-09-26

Implementation fbc666d, parent 4b7c3df, local main, version 0.246.0.
The MPI EB spray application now accepts gas slip/no-slip walls and routes
particle faces to existing specular reflection. Explicit faces still require
particle_outflow=true. Film impacts and KH children with these wall faces are
explicitly rejected. Existing restart context binds wall labels.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug; existing fixed
tests-enabled build/spray-amr-reactive-debug. No backend changes.

- Four-rank coupled CTU common/recursive tests: **2/2 passed**, 2.27 and 4.51 s
  (CTest total 6.83 s). A particle reflects at x-low and the embedded wall,
  heats and evaporates; chemical species exchange is nonzero. Total inventory
  including gas boundary, particle/EB wall and chemistry ledgers satisfies the
  existing 3e-12 relative bound. Invalid gas solver after private particle work
  preserves public gas, particles and ledgers.
- Application checker: **4 successes and 6 expected rejections**, no failures
  or skips. Two/four-rank continuous gas and particle outputs agree within the
  existing 3e-12 checks. Two-rank checkpoint resumes on four ranks with exact
  continuation rows. A parcel starting near z-high returns inside with negative
  z velocity; embedded-wall reflection and thermal transfer remain active.
  Rejections cover missing restart, overwrite, changed context, hydro method,
  wall labels and truncated particle payload.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_reactive_ctu-wall(-subcycle)?$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-wall-application-checks/application --wall --restart
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Existing compiler
warnings remain. No tolerances or golden data changed.
Artifacts: build/mpi-eb-wall-reactive-checks/{build.log,LastTest.log};
build/mpi-eb-wall-application-checks/{build.log,application.log} and
application/eb_spray_c_fct_5_/ beneath that directory.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Selected spray application, regridding/rebalancing with wall reflection,
film impacts, KH children and isothermal tree walls remain separate work.
No full matrix, installed audit, external solver parity or physical nozzle/
wall-interaction validation ran. Complete PeleC parity remains unproven.
