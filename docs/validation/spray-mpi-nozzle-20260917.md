# MPI nozzle and size-distribution inputs — 2026-09-17

Source: parent `c860c15` plus this change. Source hashes and logs are in
`build/spray-mpi-nozzle-checks/`. Existing GNU Fortran 13.3.0 Debug / Open MPI
4.1.6 configuration; no compiler/backend changes.

The MPI CLI now passes scalar/per-port finite outlet radius, legacy/solid/hollow
profile, swirl/spread and size-distribution options to the existing injector.
Size options expose Weibull, normal, lognormal and ChiSquared models, mass/number
weighting and truncation bounds. Numerical model formulas are unchanged.
Nondefault active configurations enter checkpoint identity; defaults preserve
previous context strings. These are parcel sources, not resolved nozzle geometry.

Completed checks:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

MPI application checks: 2/2 passed, 40.79 s; no failures/skips. New cases exercise
solid/lognormal and hollow/swirl/number-weighted Weibull ports, verify nonzero
outlet radii in their planes, bounded nonconstant diameters reconstructed from
pure-component droplet masses, and changed-rank checkpoint continuation (2↔4).
Final gas and particles match continuous runs at relative 2e-11 / absolute 1e-25.
Changed swirl rejects restart before history output. A scalar finite-radius normal
distribution case also produces four distinct drop masses. Earlier schedule,
profile, conservation and restart checks remain passing. Provenance, project
contract, Python syntax and whitespace checks pass.

Retained outputs:

- `build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_2/run-kwl0e1jd/`
- `build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_4/run-rt9ezftn/`

No new ChiSquared-specific, selected-mechanism, full-matrix, physical validation
or production-scaling result is claimed. MPI breakup/wall-film/SGS and particle
AMR/EB integration remain open.
