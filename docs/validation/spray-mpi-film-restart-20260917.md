# MPI wall-film restart — 2026-09-17

Source: parent `d6edce1` plus this change. Source hashes and logs:
`build/spray-mpi-film-restart-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

The MPI checkpoint API accepts optional film patches, film component masses and
a replicated wall-budget vector as a complete group. Mode, shape and budget
agreement are checked collectively. Film migration now supports an optional
common gather rank; checkpoint writes gather copies on root and use the existing
film/component/boundary-budget schema. Reads deserialize on root, redistribute
patches to current Cartesian owners, broadcast budgets and commit gas, particles,
films, components and metadata together. Calls without films retain the previous
schema and behavior. There is no new per-rank checkpoint file format.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Both focused entries passed (2/2, 10.50 s; no failures/skips). Each rank writes
one wall film with distinct temperature/material identity, mixed components and
wall energy budget 1.25 J. Every rank then independently reads the distributed
file using MPI_COMM_SELF, checking exact gas recovery, all film payloads and the
budget (2→1 and 4→1). A subsequent world-communicator read restores original film
ownership. Write operations preserve the caller's local film arrays. Deliberately
inconsistent film presence across ranks rejects collectively before file I/O;
state and wall budgets remain unchanged. Existing no-film restart checks pass.
Provenance, contract and whitespace checks pass; no Python files changed.

An initial overlong Fortran signature was wrapped before successful compilation.
No tolerance changes or golden data changes were needed. This is checkpoint API
evidence, not yet MPI wall-impact/application continuation. Root collection is
bounded to one million films and remains unsuitable as a scalable parallel-I/O
claim. Physical wall-film validation and the full matrix were not rerun.
