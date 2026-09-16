# Distributed mixed-film heating and evaporation — 2026-09-17

Source: parent `d1208a1` plus this change. Hashes and build/test/preflight logs:
`build/spray-mpi-film-coupled-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`advance_mpi_film_coupling` samples owner gas and temperature at stationary film
patches, evaluates existing mixed-film heat/evaporation rates and advances all
patches with a common globally bounded substep. Candidate gas feedback must pass
thermodynamic recovery on every rank. Rejected trials halve the common substep;
late failure preserves the complete input gas, film and component state.
Terminal patches are removed, survivors migrate to Cartesian owners and the
wall-energy ledger is reduced globally. `boundary_change` is the same global
vector on every rank, with energy entering gas+film positive; do not sum it again.
Failure returns zero ledger/substeps.

The API uses existing boundary-aware interpolation/deposition and current
component/film models, including trace depletion. It permits empty ranks and
bounds local film sampling to 125000 patches. Ordered mechanism/liquid tables and
physical model settings are caller-consistent inputs. Synchronous patch feedback
differs from the serial sequential-patch update; trajectory parity is not claimed.
This API does not yet connect wall impact, gas advection boundaries, film restart
or the MPI application lifecycle. Lateral film transport remains absent.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

2/2 focused MPI entries passed in 10.30 s, no failures/skips. A 1e-9 kg mixed
patch on the lower-x wall, initially 330 K, uses a 400 K wall and conductivity
0.6 W/(m K). Non-owner ranks start empty. A 1e-3 s interval with a one-substep
limit rejects with exact gas, temperature, component and patch rollback and zero
output ledger. A 1e-8 s accepted update increases patch temperature, reduces film
mass by evaporation, and closes gas+film+wall mass/energy/species budgets within
relative 2e-12 with a 1e-12 denominator floor. Existing MPI integration tests
also pass. Provenance, contract and whitespace checks pass; no Python files changed.

No full matrix, physical wall-film validation, multi-patch trajectory parity or
production scaling is claimed.
