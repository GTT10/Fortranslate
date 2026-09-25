# MPI EB KH children at domain walls — 2026-09-26

Implementation `2458429`, parent `c1450ea`, local main, version 0.246.0.
Owner-local particle stepping now permits KH child generation with explicit
nonperiodic outflow/slip/no-slip faces. Parents reaching a reflecting domain
face skip KH release on that substep, then use the existing EB/domain drift.
Eligible parents and children use the normal conservative shedding, feedback,
reflection and migration transaction. Periodic KH children and explicit-face
film combinations remain rejected. Application wall guards now match this scope.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, fixed H2/O2,
existing tests-enabled `build/spray-amr-reactive-debug`.
The direct application checker passed **3 successful executions and 1 expected
context rejection**: uninterrupted MPI4, MPI2 checkpoint, MPI4 continuation.
The fixture starts parents next to the upper x wall at 1000 m/s, with no-slip
x faces, other faces outflow and KH children enabled. It executes three regrids,
enables measured rebalancing and compares gas, particles, lineage/pending KH
state and budgets across restart at the unchanged 3e-12 bound.
Final population is 58 with contiguous IDs and x velocities -1000 m/s, inside
the domain. Energy absolute ledger residual is 1.23e-17 and momentum residuals
are below 1.33e-23. This is bounded conservation/restart evidence, not physical
breakup qualification. The application CTest entry was registered but not run
again separately.

Initial PCM execution produced particles but failed when writing the gas snapshot:
the shared single-level checkpoint configuration only accepts explicit boundaries
with PLM/CTU. The final fixture selects PLM/CTU and PeleC fluxes. No tolerance
or solver validation was weakened. PCM tree boundary/checkpoint compatibility
remains an open implementation issue; the failed run is retained.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-wall-khrt-checks/application --khrt --domain-wall
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace passed. Build warnings remain.
Artifacts under `build/mpi-eb-wall-khrt-checks`: `build.log`, `application.log`,
`application-initial.log`, successful `application/impact_8x3f3p84`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

No selected-mechanism, full-matrix, installed, external-solver or physical
qualification ran. Complete PeleC parity remains unproven.
