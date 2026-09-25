# Distributed periodic film thickness — 2026-09-26

Implementation `c45432ad`, parent `3976e450`, local main, version 0.246.0.
`sample_mpi_wall_film_thickness` accepts optional `periodic_length(3)`, zero
for nonperiodic axes. Positive lengths must be finite, at least 1e-12 and at
most 1e9; ranks agree on them before exchanging queries. For each film, the
query computes a bounded range of image offsets, applies the existing normal,
plane and circular-footprint tests to each image, and sums heights globally.
This includes footprints larger than a periodic length, not just the nearest
image. The candidate thickness publishes only after all queries succeed.
One film/query combination may enumerate at most one million candidate images;
excessive image spans reject transactionally before integer overflow/iteration.

The single-event MPI EB impact API forwards the optional periodic lengths to
its wet-film thickness query. Periodic flight, film field exchange and application
inputs still require integration; this API change does not enable them by itself.

## Completed focused check

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, existing fixed
H2/O2 tests-enabled `build/spray-amr-reactive-debug`.
The existing `test_spray_mpi` executable completed successfully at four ranks
(exit 0), including new analytic film-image checks. Owners and querying ranks
include empty populations. A query one period away receives the same summed
remote film height; a 1.2-wide footprint with period 0.5 receives exactly three
image contributions at its center. Rank-dependent lengths reject and preserve
output; existing invalid-normal/tolerance checks and regular MPI spray checks
also execute. The executable does not print a PASS banner. MPI EB impact target
also built with the new forwarding argument; no separate impact run was counted.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 4 build/spray-amr-reactive-debug/tests/test_spray_mpi
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. No failed attempts,
skips, tolerance changes or golden-data changes. Build warnings and X authorization
diagnostics remain. Artifacts: `build/mpi-periodic-film-query-checks/{build.log,build-eb.log,test.log}`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

No complete periodic-film application, full matrix, external comparison, scaling
or physical validation ran. Complete PeleC parity remains unproven.
