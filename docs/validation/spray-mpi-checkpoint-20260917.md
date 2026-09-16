# MPI spray checkpoint collection and redistribution — 2026-09-17

Parent `cb734ab` plus this change. Build/test/preflight logs and source hashes:
`build/spray-mpi-checkpoint-checks/`.

The new checkpoint API gathers x-slab gas and complete mixed-particle payloads on
root, writes the transactional mixed/ID schema, then can read and redistribute
into the current communicator's decomposition. Root temporary storage is bounded
to 100 million gas/temperature real values and one million particles. Candidate
arrays protect caller state on failed reads. A new optional gather destination on
particle migration reuses its existing serialization and global receive sorting.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; unchanged matching MPI configuration
as documented in `spray-mpi-migration-20260916.md`:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

2/2 entries pass (2 and 4 ranks), 9.46 s, no failures/skips. New checks write a
world-communicator checkpoint and read it on each rank's `MPI_COMM_SELF`, proving
2-to-1 and 4-to-1 redistribution with exact gas/temperature and particle
position/velocity/mass/temperature/multiplicity/component payloads. They verify
global ordering, time/step/event metadata and `huge(int64)-19`, then read back on
the original communicator and compare owner blocks and particle inventories.
Wrong context rejects without altering caller gas or particle/component state.
Earlier MPI tests pass. Provenance, contract and whitespace checks pass; no Python
source changed.

This proves API serialization and redistribution, not resumed trajectory parity
or MPI CLI restart. It is root-mediated I/O, not a scalable production parallel
checkpoint format. No full matrix or physical/scaling qualification was run.
