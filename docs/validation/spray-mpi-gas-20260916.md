# Thermodynamically checked MPI spray feedback — 2026-09-16

Tested parent: `5ac63ba` plus this change. Build/test/preflight logs and source
hashes are saved in `build/spray-mpi-gas-checks/`.

`deposit_mpi_spray_gas` validates initial gas thermodynamics, routes inventory
increments through the existing owner deposition, then recovers temperatures in
private owner blocks. Collective acceptance commits state and temperature together.
The caller supplies a common ordered species table and retains responsibility for
particle candidates, retry and the combined particle/gas transaction.

GNU Fortran 13.3.0 Debug and Open MPI 4.1.6, using the existing matching MPI build
configuration documented in `spray-mpi-migration-20260916.md`:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) passed in 0.81 s, zero failures/skips. Added checks
cover deposited heating against separate deposition and local recovery, temperature
increase, and exact preservation of all owner states/temperatures when one owner
has unrecoverable energy or negative species. Existing MPI primitive checks also
pass. Provenance, project contract and whitespace checks passed. No Python source
changed. No full matrix, coupled MPI spray application, physical or scaling
qualification is claimed.
