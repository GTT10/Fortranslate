# Distributed wall-film overlap queries — 2026-09-17

Source: parent `d8d0aef` plus this change. Source hashes and logs:
`build/spray-mpi-film-query-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`sample_mpi_wall_film_thickness` queries distributed stationary circular patches
using the same normal, plane-distance and radial-footprint criteria as the existing
serial wall-impact code. Query batches are broadcast by source rank; each owner
sums its matching film heights and a reduction returns the total to the querying
rank. No film inventory is replicated. Local limits are one million patches and
125000 queries; the common plane tolerance and query normals are validated before
communication. Failure preserves the previous output array.

The algorithm visits all owners for each source batch and scans local patches;
it is a bounded implementation, not a spatially indexed/scaling claim. Floating
summation order may change with rank count. Crossing detection, converting thickness
to film/diameter ratio and calling the impact event API remain driver work.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

2/2 focused MPI entries passed in 9.99 s; no failures/skips. Tests include an
empty query source and an empty film owner, remote overlapping patches with known
summed heights, radial and plane exclusions, opposite-normal filtering, mismatched
plane tolerance and an invalid normal on one rank. The latter two reject
collectively without changing prior results. Existing MPI integration tests also
pass. Provenance, project contract and whitespace checks pass; no Python files
changed. No full matrix or physical wet-wall validation was run.
