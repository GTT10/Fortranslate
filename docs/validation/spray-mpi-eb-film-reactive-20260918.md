# MPI EB reactive particle/film composition — 2026-09-18

Implementation `691c601`, base `4734078`.

The existing MPI EB reactive wrapper optionally carries stationary mixed films
through particle half / film half / reactive gas / film half / particle half
composition. Optional-state presence is collectively checked. Film source
updates and finest-owner migration occur on private state; full publication
includes gas, particles, films, component arrays and the three global ledgers.
Film substep counts are reported separately and reset to zero on failure.
Films are supplied by the caller; particles remain specular in this wrapper.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja, existing fixed tests-enabled and selected H2O2 tests-disabled
configurations. The new MPI test enables runtime checks and floating-point traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film_reactive test_spray_mpi_eb_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(film_reactive|reactive)_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **4/4 passed**, total **4.27 s**. New film tests: two ranks 1.06 s,
four ranks 1.11 s. Existing no-film tests: 1.05/1.02 s. Selected library build,
provenance, Python syntax, contract and whitespace passed. No test failures or
skips; selected execution was not run.

Two development build errors were corrected: an overly broad CMake insertion
duplicated neighboring targets, and a new test call exceeded the Fortran line
length. Their logs are retained. Neither changed production numerical semantics
or tolerances.

The two-level cut-cell fixture has one reflecting mixed parcel and two mixed
films initially held on one rank. Films migrate to root/fine owners and heat
and evaporate in both source halves. HLLC, implicit H2O2 chemistry and all
molecular transport switches are active. Composite gas/particle/film inventories
close against exterior, wall and chemistry ledgers within relative 3e-12, floor
1e-30. The chemistry species ledger is nonzero.

An invalid gas solver after the first particle and film halves confirms rollback
of gas, temperatures, particle state, film count/components/mass/temperature and
diagnostics. Rank-disagreed film presence rejects collectively before updates.
Successful steps report at least two film substeps; failures report zero.

Artifacts: `build/spray-mpi-eb-film-reactive-checks/{build-initial.log,build-line-error.log,build.log,build-final.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
Distributed impact-generated films, MPI application input/output/restart, scaling,
external parity and physical validation remain unqualified. No broad matrix or
installed-binary audit was run.
