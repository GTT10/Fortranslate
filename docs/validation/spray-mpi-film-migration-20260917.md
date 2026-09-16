# MPI mixed wall-film patch migration — 2026-09-17

Source: parent `3c761cb` plus this change. Source hashes and build/test/preflight
logs are in `build/spray-mpi-film-migration-checks/`. Existing GNU Fortran 13.3.0
Debug / Open MPI 4.1.6 configuration.

`migrate_mpi_wall_films` routes stationary film patches and their whole-patch
component masses to Cartesian cell owners. Closed upper boundaries select the
last cell while the stored physical position stays exact. Out-of-domain positions
are rejected; no periodic wrapping is applied. Position, unit normal, mass,
temperature, diameter, height and material identity travel together. The routine
accepts 1--16 components and up to one million local patches, validates common
geometry/partition/capacity, exchanges counts before payloads and commits only
after all ranks accept their received inventories.

Receive ordering is source rank then original local order. Films currently have
no global ID; partition-independent ordering and lateral film transport are not
claimed. This is a migration API, not yet MPI wall-impact, gas-film feedback,
nonperiodic gas boundary or checkpoint/application integration.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Both focused MPI entries passed (2/2; 9.97 s; no failures/skips). New film cases
check collective capacity rejection without changing patch/component arrays,
upper-x and upper-y wall ownership, exact geometry/material/composition payloads,
a second redistribution leaving empty ranks and global component inventory
conservation at relative 2e-12. Existing MPI integration cases also pass.
Provenance, project contract and whitespace checks pass. No Python files changed.
No full matrix, physical wall-film validation or production scaling was run.
