# MPI hydro with physical boundaries — 2026-09-17

Source: parent `4b81db2` plus this change. Source hashes and logs are in
`build/spray-mpi-boundary-hydro-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`advance_mpi_spray_hydro` accepts optional face labels, inlet states/temperatures
and a boundary-change output. With face labels supplied, each SSPRK2 stage fills
two ghost layers and uses the serial characteristic-PLM face-flux kernel.
Only global x faces and local transverse physical faces contribute to the global
surface integral; MPI interfaces do not enter the budget. The returned ledger is
replicated on every rank and positive into the gas. Gas, temperature and ledger
commit only after both stages and the collective reduction succeed. Omitting
face labels retains the existing periodic slab path. Boundary data use the same
caller contract as the ghost API; labels and mode are checked collectively.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

The first build caught an undeclared spacing variable in the new test; adding its
unit-grid parameter fixed compilation. Final focused tests passed 2/2 in 12.65 s,
with no failures or skips. Uneven 2/4-rank partitions compare periodic,
wall/outflow and inlet/wall updates against the full serial boundary solver:
conserved fields within 2e-12 scaled error, temperature within 1e-8 K and boundary
ledger within 1e-10 relative error with a 1e-12 floor. Independent volume inventory
changes agree with the surface ledger within 2e-12 of initial absolute inventories.
An excessive step fails during numerical advancement, preserves gas/temperature
exactly and returns a zero ledger. A rank-local face mismatch also rolls back.
Existing periodic hydro, coupling, film and restart API checks passed in these
same entries. Provenance, contract and whitespace checks passed; no Python changed.

This does not connect nonperiodic molecular transport, SGS, particle movement or
application inputs. No full configuration matrix or physical validation is claimed.
