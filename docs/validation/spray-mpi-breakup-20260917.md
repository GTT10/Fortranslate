# Representative MPI mixed-drop breakup — 2026-09-17

Source: parent `8e3e326` plus this change. Hashes and logs are in
`build/spray-mpi-breakup-checks/`. Existing GNU Fortran 13.3.0 Debug / Open MPI
4.1.6 build; no compiler/backend changes.

The distributed mixed-particle coupling API now applies representative ETAB or
KH-RT breakup after each frozen-gas particle update. It rebuilds effective liquid
properties and the mixture film at the updated composition/temperature, calls
the existing breakup model and rescales per-drop component masses. Multiplicity
preserves whole-parcel inventories. Candidate particle, composition and gas state
remain private through collective feedback acceptance and end-of-interval migration.
KH child shedding is still explicitly rejected because it requires new distributed
ID reservations and child inventories. The MPI CLI still uses breakup disabled;
checkpoint persistence for active breakup is a separate remaining connection.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Both MPI entries passed (2/2; 9.87 s; no failures/skips). New ETAB/KH-RT cases use
1000 m/s parcels over 1e-8 s with heat, drag and evaporation disabled to isolate
breakup. ETAB starts from deformation 1.1. Both models reduce physical drop mass,
increase multiplicity and retain parcel counts. Every component inventory agrees
within relative 2e-12, and global gas/particle mass, momentum, energy and species
inventories agree within relative 2e-12 with a 1e-12 denominator floor. Existing
migration, coupling rollback, reaction, transport and checkpoint checks in the
same executable also pass. Provenance, contract and whitespace checks pass.

No Python files changed. No active-breakup checkpoint/application run, child
shedding, full matrix, physical validation or scaling qualification is claimed.
