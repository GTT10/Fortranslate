# Spray CVODE build repair — 2026-09-16

Parent `c46c325` fails to build `test_spray_integration` with spray and SUNDIALS
enabled. The CVODE branch of `advance_spray_chemistry_3d` contained mixed wall-film
setup referencing variables that belong to `advance_spray_reactive_3d`, not to
this chemistry-only procedure. Compiler diagnostics report undeclared
`component_mass`, `wall_films`, `film_component_mass` and `mixture_film_mass` and
an invalid allocation. Removing that misplaced block restores compilation; the
mixed-film setup in the coupled driver remains in place.

GNU Fortran 13.3.0, Release, MPI OFF, SUNDIALS ON using the existing
`build/port-spray-cvode` configuration and SUNDIALS 7.2 Fortran install at
`build/sundials-7.2-fortran/install`.

```sh
cmake --build build/port-spray-cvode --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-cvode -R '^spray_integration_(cvode|native)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Before: build fails with the diagnostics above. After: build passes, and 2/2
focused chemistry entries pass, no failures/skips (13.80 s). Each entry advances
125 cells three times and checks actual reaction progress, nonnegative species,
species closure, inert nitrogen and unchanged conserved hydrodynamic quantities.
Provenance, contract and whitespace checks pass. No Python source changed.

Before/after build logs, final tests, preflight and source hash are retained in
`build/spray-cvode-repair-checks/`. No full matrix, new mixed-wall CVODE physical
qualification or MPI reactive spray integration is claimed.
