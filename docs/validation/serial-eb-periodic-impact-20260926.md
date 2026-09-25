# Serial EB impact boundary flight — 2026-09-26

Implementation `4d22e6dd`, parent `7818a514`, local main, version 0.246.0.
The serial `advance_mixture_level_set_impacts` interval now accepts optional
`outer_faces` alongside the exterior budget. It validates periodic pairing and
opposing nodal values, processes the first outer event before continuing flight,
and supports periodic wrap and specular walls as well as existing outflow.
Both outer and EB events consume the contact-depth limit. Failed intervals
preserve particles, films, IDs and ledgers. Periodic crossing creates no exterior
transfer; reflection records impulse separately from EB deposition.

Periodic footprint-image enumeration is now shared by serial impact thickness
and the existing MPI film query through `count_periodic_film_images`. It retains
the one-million candidate-image cap and includes overlapping periodic copies.
MPI combines image multiplicity with film height before its distributed sum.
Single-level particle/film coupling and application input remain to be connected.

## Focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, existing fixed
H2/O2 tests-enabled `build/spray-amr-reactive-debug`.
Serial impact CTest passed **1/1**, total 0.05 s. New cases cross the upper z
face before depositing on x=0.5: periodic flight deposits at z=0.09 and reflected
flight at z=0.91. Both count two events, conserve inventories at unchanged 1e-12
checks, and preserve candidate state/diagnostics when limited to one event.
Existing impact/rebound/splash and outflow checks also execute.
The existing MPI spray executable completed at four ranks (exit 0), including
analytic remote periodic-film and multiple-image sums plus mismatch rollback.
It emits no PASS banner. No failed attempts or skips occurred.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact test_spray_mpi --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 4 build/spray-amr-reactive-debug/tests/test_spray_mpi
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace passed. Build warnings and runtime X
authorization diagnostics remain. No tolerances or golden values changed.
Artifacts: `build/serial-eb-periodic-impact-checks/{build.log,ctest.log,LastTest.log,mpi.log}`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
No application, selected-mechanism, full-matrix or external/physical validation
was run for this serial kernel change. Complete PeleC parity remains unproven.
