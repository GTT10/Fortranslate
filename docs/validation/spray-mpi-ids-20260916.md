# Collective particle ID reservation — 2026-09-16

Tested parent: `07bd853` plus this change. Source/test hashes and build/test/preflight
logs are in `build/spray-mpi-id-checks/`.

`reserve_spray_ids` validates local request counts and a shared 64-bit high-water
mark, sums requested IDs and uses an exclusive rank prefix for disjoint intervals.
It commits the new mark and first local ID only after successful collectives and
overflow checks. Zero requests return first_id=0. The caller owns checkpointing of
the mark and wiring to injection/breakup. Rank-layout-independent identity is not
claimed; allocation order follows communicator ranks.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; configuration and resolved-library boundary
are as recorded in spray-mpi-migration-20260916.md:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) passed in 0.82 s total, no failures/skips. New cases
cover unequal requests including zero, successive reservations, negative-count
collective rollback, 64-bit overflow and inconsistent existing marks. The existing
migration payload/ownership/capacity checks also pass. Provenance, contract and
whitespace checks passed; no Python source changed. No full matrix, integrated MPI
spray application or physical/scaling qualification was run.
