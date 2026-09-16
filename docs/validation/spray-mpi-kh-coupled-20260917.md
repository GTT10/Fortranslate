# MPI KH children in coupled stepping and restart — 2026-09-17

Source: parent `148fa2a` plus this change. Source hashes and build/test/preflight
logs: `build/spray-mpi-kh-coupled-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

The mixed-particle coupler now advances KH shedding after frozen-gas transfer,
adds the children's transverse kinetic-energy debit to owner-routed gas feedback,
and grows subsequent substep sampling arrays for the children. All updates,
including the ID watermark, stay private until the entire interval succeeds.
The reaction/hydro/transport wrapper also carries a private watermark through both
particle half-steps. Consumed parents are removed before final spatial migration.
Child generation requires the optional `last_id` argument; candidate populations
are limited to the existing 125000-per-rank sampling/deposition bound.

MPI application inputs now enable children and their reference-diameter threshold.
Checkpoint flags bind breakup/RT/shedding mode collectively and preserve parent
IDs and KH inventories through the existing serial schemas. CSV output exposes
these fields. Default-disabled contexts are unchanged.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi pelef_mpi_spray_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_(application|migration)_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

API tests start with an empty rank, force whole-interval substep exhaustion and
verify restoration of gas, temperature, component masses, parent state and IDs.
Accepted child-release updates conserve global gas/particle mass, momentum,
energy and species within relative 2e-12 with a 1e-12 floor, and children survive
migration. The event API's tighter component-specific conservation checks remain.

Application tests use 1000 m/s injection, KH children and a 1e-25 m3 reference
threshold, through 2e-8 s. They compare continuous output against a checkpoint at
1e-8 s followed by 2↔4-rank continuation, including child IDs, parent lineage,
KH inventory and gas state (relative 2e-11 / absolute 1e-25). The final population
exceeds four parcels and has positive parent IDs. Corrected mass/energy errors
remain below 2e-11. A changed threshold rejects restart before history creation.

Future child IDs use rank-prefix allocation. This application's releasing
population stays on one owner under either decomposition; arbitrary repartitioned
simultaneous releases are not claimed to produce identical future IDs. Physical
spray validation, wall-film/SGS integration, AMR/EB particles and production
scaling remain open; no full configuration matrix was rerun.

Final result: 4/4 focused entries passed in 70.96 s; no failures/skips.
Provenance, project contract, Python syntax and whitespace checks passed.
