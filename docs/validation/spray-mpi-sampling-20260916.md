# Owner-routed particle field sampling — 2026-09-16

Tested parent: `8fec6bd` plus this change. Source/test hashes and build/test/preflight
logs are in `build/spray-mpi-sampling-checks/`.

The MPI sampler builds the existing eight-cell stencil, sends requested indices to
Cartesian owners and returns field values in the request order. The weighted sum
uses the serial stencil order. No global field gathering or root processing is used.
Each rank supplies exactly its own cell block; global layout and boundary metadata
must agree. The block helper handles uneven remainder-cell allocation. Invalid
input/received indices reject collectively, retaining prior sample arrays. The API
is bounded to 125000 positions, 256 fields and one million received requests per
rank. Distributed feedback and coupled stepping remain open.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; configuration/library boundary is retained
from spray-mpi-migration-20260916.md:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI entries (2 and 4 ranks) passed in 0.82 s total, no failures/skips. New cases
compare two synthetic fields exactly against serial interpolation across uneven
blocks and periodic/physical boundaries. Empty requesting ranks and a one-rank
out-of-domain position with collective result preservation also pass. Existing
migration, injection, splitting and ID cases pass. Provenance, contract and whitespace
checks passed; no Python source changed. No full matrix, distributed reactive spray
application or physical/scaling qualification was run.
