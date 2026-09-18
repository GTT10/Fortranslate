# MPI AMR spray uniform acceleration — 2026-09-18

Implementation `34b6cd54688fcf74209a0c1aa6b504fca574c88d`, based on `ac3e6d1`.

The MPI AMR reacting-spray API accepts optional uniform acceleration and a
separate source ledger. Coefficients are finite, bounded by `1e6` per component
and agree collectively. Nonzero acceleration requires a correctly sized ledger.
Gas and parcels receive symmetric half kicks inside the particle halves and
outside the LES/gas operators. Each kick changes momentum by `mass*dv` and total
energy by `dot(momentum_old + 0.5*mass*dv, dv)`, preserving internal energy.
Gas EOS recovery checks temperature consistency. Covered coarse cells update
but only exposed cells contribute to the volume integral; parcel multiplicity
weights particle impulse/work. MPI reduces the separate source inventory.
All caller state and diagnostics remain transactional. Attached films are not
accelerated. Application inputs and source-ledger restart are not included.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `test_spray_mpi_amr_reactive` and selected H2O2 tests-disabled
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.
Final CTest **1/1 passed**, 60.88 seconds (60.92 s total), zero failures/skips.
The initial version also passed the same test in 60.46 s before adding the
later-failure rollback assertion and final finite-source guard; this is one
distinct test rerun after changes. Provenance, Python syntax, project contract
and whitespace checks passed. No failed builds/tests occurred.

The nested three-node mixed-liquid fixture runs on two ranks. In addition to
the existing default-disabled conservation/physical-boundary checks, uniform
acceleration `[100,-200,300]` checks gas-plus-particle change against external
impulse/work at the existing `3e-12` relative tolerance. Total impulse also
matches `total_mass*dt*acceleration` at that tolerance. Rank mismatch rejects
with unchanged source and state; an invalid gas solver after the first kick
also preserves caller state and source sentinels. No tolerances changed.
Force with active LES, films, chemistry or physical boundaries was not separately
qualified. No external solver parity, full matrix, scaling or physical validation.

## Commands and artifacts

Logs in `build/spray-amr-force-checks/`: `provenance.log`, `build.log`,
`build-final.log`, `build-selected.log`, `build-selected-final.log`, `ctest.log`
and `ctest-final.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_force_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
