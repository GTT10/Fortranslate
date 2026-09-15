# Collective cone injection — 2026-09-16

Tested parent: `1a2c00e` plus this change. Source/test hashes and build/test/preflight
logs are in `build/spray-mpi-injection-checks/`.

The MPI cone injector reserves candidate ID ranges, uses the existing size/nozzle
injector and commits particles, optional component columns and the ID mark only when
all ranks succeed. Mixed injection evaluates effective density/cp/energy from the
component table and local fractions. Returned budgets are local conserved/species
amounts. Geometry, count, capacity or material failure on one rank returns zero
budgets without consuming IDs or replacing arrays elsewhere. Injection does not
implicitly migrate; the caller handles spatial ownership and scheduling.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; configuration/library boundary is retained
from spray-mpi-migration-20260916.md:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) passed in 0.81 s total, no failures/skips. New checks
cover unequal local injection counts including zero, rank-local blend fractions,
globally disjoint IDs, injected total/species mass, component inheritance and an
invalid-axis failure after ID reservation with collective state/ID rollback.
Existing migration/splitting/ID cases also pass. Provenance, contract and whitespace
checks passed; no Python source changed. No complete matrix, distributed CFD spray
application or physical/scaling qualification was run.
