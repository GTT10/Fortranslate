# Distributed particle feedback deposition — 2026-09-16

Tested parent: `5b4ea7f` plus this change. Source/test hashes and build/test/preflight
logs are in `build/spray-mpi-deposition-checks/`.

The MPI deposition API routes eight-cell stencil indices and weighted inventory
increments to owners, divides by volume and adds to local field candidates. Sending
rank and local request order determine accumulation order. Invalid input, bad routing
or arithmetic overflow rejects collectively; old local fields remain intact. This
is a generic field operation, not gas thermodynamic recovery or a completed coupled
particle step. Sampling/deposition now also validate spacing with empty input ranks.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; configuration/library boundary is retained
from spray-mpi-migration-20260916.md:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) passed in 0.83 s total, no failures/skips. New checks
compare positive/negative component contributions exactly with serial accumulation
across uneven blocks and periodic/physical boundaries. Empty depositing ranks and
one-rank invalid-position rollback also pass. Existing sampling, migration, injection,
splitting and ID cases pass. Provenance, contract and whitespace checks passed; no
Python source changed. No full matrix, distributed reactive spray application or
physical/scaling qualification was run.
