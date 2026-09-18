# Collective EB contact dispatch — 2026-09-18

Implementation `580106d`, base `c753e79`.

The new MPI contact dispatcher joins nodal first-contact geometry, remote film
thickness queries and the existing collective impact resolver. It calculates
the same impact parameters as the planar/serial paths, reserves child IDs
collectively and reports a global wall-received budget. Contact flags and
fractions remain associated with the original population. Model-range failures
propagate collectively and preserve caller state.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled
configurations were reused. The test enables runtime checks and invalid/zero/
overflow floating-point traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **0.97 s** (two ranks 0.46 s; four ranks 0.47 s).
Selected runtime library compilation, provenance, Python syntax, contract and
whitespace passed. No failures or skips. Selected execution was not run.

Each rank sends one synthetic water parcel toward x+y=1 at 200 m/s. A stationary
film held only on the final rank supplies thickness/drop-diameter 0.03 to all
queries. Every impact produces four wet-splash children; IDs follow rank-prefix
reservation and the shared watermark becomes five times the rank count. The
contact fraction is 0.5 within 1e-14. Global liquid/film mass, momentum and energy
plus wall receipt close within relative 3e-12, component-wise floor 1e-20.
Parcel and film component payloads equal their single-component masses exactly.

Insufficient parcel capacity and a rank-disagreed nodal value preserve payloads,
watermark and diagnostic sentinels. Increasing only the remote film to ratio
0.2 triggers collective model rejection with no publication. A further dry
impact with a parcel only on the final rank confirms empty-rank participation
and shared child-ID reservation.

Artifacts: `build/spray-mpi-eb-impact-checks/{build.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
This is a single-contact batch. Films are sampled before the batch; newly
created patches do not change other batch contacts. Remaining flight, dead-parent
compaction, final migration, source/reactive composition and application
persistence are not connected here. No broad matrix, scaling, installed audit,
external parity or physical validation was run.
