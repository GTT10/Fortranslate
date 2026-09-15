# Distributed mixed particles and gas reaction — 2026-09-16

Tested parent: `63a8645` plus this change. Build/test/preflight logs and source
hashes: `build/spray-mpi-reaction-checks/`.

`advance_mpi_mixture_reaction` places cell-local chemistry between two periodic
mixed-particle half-steps. All four caller inventories (gas, temperature, particles,
component mass) remain unchanged until every stage succeeds. Interval, tolerances
and chemistry method agree collectively. Ordered model tables remain a shared
caller convention. This does not yet advance gas advection or transport.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, existing matching MPI configuration
as documented in `spray-mpi-migration-20260916.md`:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) pass in 9.07 s, no failures/skips. The added
1400 K reactive mixture case advances hydrogen through native implicit chemistry
while retaining global gas-plus-particle mass, momentum and total energy
(`2e-12 * max(1e-12, abs(initial))`) and H/O/Ar/N atom inventories
(`2e-10 * max(1e-12, abs(initial))`). A negative Arrhenius prefactor on only the
last rank rejects chemistry after the first particle half-step; initial gas,
temperature, IDs, positions, velocities, particle temperature and liquid component
masses remain exact. Earlier MPI primitive and coupling checks also pass.

Provenance, project contract and whitespace checks pass. No Python source changed.
No new MPI CVODE qualification, full reactive CFD application, changed-rank
trajectory parity or physical/scaling validation is claimed.
