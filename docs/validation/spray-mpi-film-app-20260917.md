# MPI application wall impacts, films and lineage restart — 2026-09-17

Source: parent `3a2a8b1` plus this change. Source hashes and logs:
`build/spray-mpi-film-app-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration. Fixed H2/O2 frontend exercised here; selected
mechanisms share the application implementation but were not separately rerun.

The `mpi_spray` namelist now accepts derived `impact` settings. Enabled faces
activate coupled deposition/splash and optional stationary-film heating/evaporation.
All wall-model fields and the film evolution identity extend checkpoint context.
Total inventories include films; history appends film count and mass. Per-rank
`.films.csv` files retain geometry, material, thermal state and component masses.
The example is `cases/mpi_spray_3d/wall_film.nml`.

Impact checkpoints use the opt-in `_RET_V1` schema extension. It stores the
existing parent ID field while allowing the referenced parent to have been
consumed. It requires an ID watermark and a lineage-carrying film/split schema;
parcel validation retains `0 <= parent_id < child_id <= last_id`. Strict-parent
schemas retain their prior parent-presence requirement, and the new schema cannot
be silently read as an old one. Root gathering and rank-neutral redistribution
remain unchanged. No historical parent record is fabricated.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_(application|migration)_[24]$' --output-on-failure
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 -m py_compile tests/spray/test_mpi_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

The application/API run passed 4/4 in 128.89 s. New deposition and splash examples
create nonzero film mass and heat films above 330 K. Low-speed parents disappear;
faster impacts leave more than four descendants. Film component masses sum to
total mass within relative 2e-12. Same-rank split/restart reproduces every gas,
particle and film CSV value exactly. Rank-count changes 2-to-4 and 4-to-2 agree
within relative 2e-11 and absolute 1e-25, including splash parent IDs; mass/energy
history errors remain below 2e-11. A conductivity change rejects restart before
history creation. Detailed output paths are retained in `app-artifacts.log`.

The first application run passed deposition restart but rejected splash checkpoint
writing: consumed parents were absent from the live particle set. The explicit
retired-parent schema addresses that lifecycle requirement; the failed run is
retained in `ctest-parent-presence-failure.log`. The existing film schema already
stored parent IDs, so no independent split-state extension was retained. A build
also caught overlong new declarations; they were wrapped. No numerical tolerance
or reference hash was changed.

The final focused guards passed 2/2 in 17.12 s, with no failures or skips.
Python syntax, provenance, contract and whitespace checks passed. These cases check strict missing-parent rejection,
self-parent rejection even in retired mode, schema mismatch rejection and successful
retired-parent round-trip. Synchronous MPI impact ordering and pre-batch film
sampling still differ from serial sequential impact. No arbitrary future-child-ID
parity, physical wet-wall validation, production scaling or full matrix is claimed.
