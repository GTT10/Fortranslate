# Periodic MPI EB impact and film integration — 2026-09-26

Implementation `6b1ffae0`, parent `ad9976c9`, local main, version 0.246.0.
MPI impact flight accepts paired periodic faces with matching opposing nodal
geometry. It resumes remaining flight on the opposite face, with zero periodic
exterior transfer. EB impacts use periodic lengths for wet-film thickness.
The coupled split forwards periodic axes to film heating/evaporation; existing
wrapped EB field sampling/deposition supplies conservative gas exchange.
Application restrictions are removed for this connected path. Physical outflow
and specular walls retain their own budgets. Outer-domain film creation is not
added; films remain on stationary embedded surfaces.

## Focused execution

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, fixed H2/O2,
existing tests-enabled `build/spray-amr-reactive-debug`.
The direct application checker passed **6 successful executions and 2 expected
context rejections**. Deposition and splash each run uninterrupted MPI4,
MPI2 checkpoint and MPI4 continuation, with CTU, three regrid stages and
measured-rebalance controls. Parents cross the upper x periodic face before
hitting the y-normal EB. Two films are formed at x approximately 9e-10 m;
the deposited mass heats above its initial 350 K, and splash creates eight
children. Gas, particle, film and budget continuation agree within the unchanged
3e-12 comparison, including components, IDs and film state. Maximum absolute
budget residual is 1.496064e-17.

No failed attempts or skips occurred. The new CTest application entry is
registered, but the checker ran directly and was not counted a second time.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-periodic-film-checks/application --periodic
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Build warnings remain.
No tolerance or golden-data changes. Artifacts under
`build/mpi-eb-periodic-film-checks`: `build.log`, `application.log`,
`application/impact_7ipxkjou`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

No selected-mechanism, full matrix, external comparison, scaling or physical
impact/film qualification ran. This MPI tree path does not qualify the separate
single-level serial periodic-film path. Complete PeleC parity remains unproven.
