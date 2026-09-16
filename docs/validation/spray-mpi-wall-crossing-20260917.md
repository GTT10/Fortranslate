# MPI planar wall crossing and impact conditions — 2026-09-17

Source: parent `cbfbad9` plus this change. Hashes and logs:
`build/spray-mpi-wall-crossing-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`impact_mpi_planar_walls` compares previous in-domain positions with unwrapped
candidate positions, selects the earliest crossed face with the serial axis tie
order and processes enabled stationary walls. It computes mixture density and
boiling temperature, temperature-dependent viscosity, We/Re/Kv and wall-temperature
criterion, then queries distributed film thickness and calls the collective
mixed-impact resolver. Non-crossing parents stay unchanged. Failed validation or
impact resolution preserves particle/film/component/ID inputs and the prior hit
mask; wall output is zero. The returned hit mask indexes the original parents.

Ordered components, molecular weights, physical settings and pressure samples are
caller inputs. Queries use the film state at the start of the event batch; films
deposited by other simultaneous events do not change that batch's impact criteria.
This differs from serial sequential wet-film accumulation. Crossings of disabled
faces, post-impact migration/removal and remaining-time advancement stay with the
caller. The gas/particle coupled driver and application are not yet wired to this
nonperiodic crossing API.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

2/2 focused MPI entries passed in 10.07 s; no failures/skips. New cases compare
one lower-z crossing and one non-crossing mixed parent per active rank, with an
empty final rank. Deposition gives exact serial parent mass, component masses,
film position/height and deposited composition; the globally reduced wall ledger
matches summed serial impact results within relative 1e-13 with a 1e-12 floor.
No splash IDs are consumed by pure deposition. An out-of-domain previous position
on one rank rejects collectively and preserves the existing hit mask and ID mark.
Existing distributed overlap, splash, film, restart and coupling checks also pass.
Provenance, contract and whitespace checks pass; no Python files changed.

No full matrix, arbitrary wall geometry, moving-wall impact, wet-impact physical
validation or production scaling is claimed.
