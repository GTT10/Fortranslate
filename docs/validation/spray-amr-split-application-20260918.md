# MPI AMR split-parcel application restart — 2026-09-18

Test implementation `d023038d0c5009a0e3c7d81c53fd80f00ea0d4da`, based on
`ae1faf0`. Solver splitting implementation remains `3a14bfb`.

The AMR application restart checker now exercises computational parcel
splitting through existing CLI options, scheduled injection, coarsening and
checkpoint continuation. This closes the application-level evidence gap from
the preceding coupled-library increment; no numerical implementation changed.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, existing fixed configuration.
`pelef_mpi_amr_spray_3d` rebuilt. CTest **1/1 passed**,
`spray_mpi_amr_split_application`, 4.45 seconds (4.49 s total), zero failures/skips.
Selected executable was not rebuilt or run in this increment.

The mixed-liquid case uses ETAB, high surface tension, multiplicity limit 5,
and disabled drag/heating/evaporation. Three injections create 12 original
parcels. Computational splitting produces 40 children and 52 final parcels,
each with multiplicity in `(0,5]`. All 52 IDs are unique, with watermark 52.
The intermediate checkpoint contains 20 parcels and watermark 20, including
the next injection at the checkpoint time. Total represented mass remains
`1e-10 kg` and each parcel retains water mass fraction 0.25 within `2e-12` relative.

Continuous 2-rank gas/parcel outputs match 2-to-4 rank continuation within the
existing `2e-12` relative / `1e-25` absolute tolerances, including IDs and parent
IDs. Changed multiplicity limit rejects restart; existing changed-input,
truncated-payload and immutable-checkpoint rejection checks pass. These results
cover this scheduled topology/ownership pattern, not arbitrary rank-dependent
ID allocation order, all breakup models, or physical breakup validation.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, tolerance changes, full matrix, external solver parity
or scaling qualification.

## Commands and artifacts

Logs: `build/spray-amr-split-application-checks/{provenance,build,ctest}.log`.
Run artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_split_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_split_application$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
