# Collective mixed-liquid wall-impact events — 2026-09-17

Source: parent `c1a493a` plus this change. Source hashes and build/test/preflight
logs: `build/spray-mpi-impact-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`resolve_mpi_wall_impacts` resolves a strictly increasing list of local parent
indices using supplied locations, normals and Kv/Tstar/Kcrit/film-ratio/We/Re.
It calls the existing classifier to count splash children, reserves global IDs,
and runs the existing mixed-liquid event model into candidate arrays. All ranks
accept before parents, children, deposited films, component matrices and ID state
are committed. The returned five-component wall ledger is globally reduced and
replicated; its sign is transfer out of particles/film into the wall. Capacity or
late model failure leaves inputs unchanged and returns zero wall ledger.

This event API does not detect crossings or query distributed wet-film overlap.
Those remain caller responsibilities. Consumed parents remain as zero-mass
records for the caller to remove, and new films/children still require spatial
migration. Existing wet-impact correlation limits remain unchanged. No automatic
wall-boundary application support is claimed.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

New tests use one mixed parent per rank except the empty final rank. A dry impact
with Kv=150, Tstar=1, Kcrit=80, film ratio=0, We=100 and Re=1000 produces four
children and one deposited patch. Tests check rank-prefix IDs, parent lineage,
child composition, per-component inventory closure at relative 2e-12, and global
mass/momentum/energy/species including the wall ledger (relative 2e-12 with a
1e-12 floor). Setting Re=0 forces a late event-model failure after ID reservation;
parent/component/film/ID state remains unchanged and the ledger is zero.

No Python files changed. No full matrix, wet-wall physical validation or production
scaling is claimed.

Final result: 2/2 focused MPI entries passed in 9.66 s, no failures/skips.
Provenance, project contract and whitespace checks passed.
