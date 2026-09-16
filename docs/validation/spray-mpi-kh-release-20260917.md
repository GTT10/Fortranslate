# Collective mixed-liquid KH child release — 2026-09-17

Source: parent `711af3b` plus this change. Hashes and logs are retained in
`build/spray-mpi-kh-release-checks/`. Existing GNU Fortran 13.3.0 Debug / Open MPI
4.1.6 build.

`advance_mpi_khrt_shedding` advances the existing KH accumulator and constructs
four symmetric children per releasing parent. It builds effective mixed-liquid
properties, preserves component fractions, reserves globally unique IDs and
returns gas energy increments indexed by original parent. Capacity, ID or model
failure leaves parcels, component masses, ID watermark and the prior energy
array unchanged. Empty ranks participate in every collective. Child IDs depend
on rank-prefix reservation; cross-partition ID parity is not claimed.

This is an event API. A coupled driver must deposit its energy ledger into a
candidate gas state and commit gas/particles together. The existing full MPI
coupler still rejects child shedding pending that connection. No gas deposition,
CLI or shedding-checkpoint integration is claimed by this increment.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Focused tests use one mixed parent on each rank except the empty final rank,
gas density 10 kg/m3, gas speed 300 m/s, Reynolds number 2000 and interval 1e-6 s.
They check collective capacity rollback, four-child release, rank-prefix IDs,
lineage, negative gas energy transfer, inherited composition and component
inventory conservation at relative 2e-12. Coupled inventory accounting includes
the returned energy ledger; total mass, momentum, energy and species use relative
2e-12 with a 1e-12 floor. Existing MPI integration checks also run.

No Python files changed. No full matrix, physical validation or scaling
qualification is claimed.

Final focused result: 2/2 passed in 9.66 s, no failures/skips. Provenance,
project contract and whitespace checks passed.
