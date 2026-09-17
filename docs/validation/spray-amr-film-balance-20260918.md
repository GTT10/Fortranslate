# MPI AMR film-aware ownership balancing — 2026-09-18

Implementation: `38b9f07c2bbf79f9929cf4ec017c752e5ea2487d`, based on `1c20093`.

The optional film arrays/component masses/capacity group in
`rebalance_mpi_amr_spray` extends atomic publication to stationary films.
Presence agrees collectively and film inventories are validated before balancing.
After private gas and particle redistribution, film migration uses the same new
node owners. Any failure preserves all original arrays and owners and returns
`changed=.false.`. If gas costs select no change, film ordering/storage remain
untouched. The cost model remains gas work; no particle/film timing was added.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_balance --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_balance_(2|4)$' --output-on-failure --output-log build/spray-amr-film-balance-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips: two ranks **0.56 s**, four ranks **0.57 s**,
total **1.16 s**. The existing three-level test concentrates gas and particles
on rank zero, adds a mixed film at the nested lower wall and supplies costs
`[100,80,60]`. The new owner receives the film with exact position, normal, mass,
temperature, diameter, height, material and component masses. Existing exact gas
redistribution and particle ownership/composition assertions remain active.

A zero film capacity rejects after private gas/particle redistribution and
preserves owners, gas/temperature, particles and film state. Particle-capacity
failure is also checked with films present. Zero costs produce a successful
no-change result, preserving all inventories. The fixture retains earlier
film-free remap checks as setup/regression evidence.

Fixed and selected builds, provenance, syntax, project contract and whitespace
checks passed. No build/test attempt failed; existing compiler warnings remain.
Selected evidence is compilation only. No Python source changed.

Artifacts: `build/spray-amr-film-balance-checks/{build,build-selected,ctest}.log`.
Executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_balance`.

This is a library ownership operation. Impact-created AMR films and application/
checkpoint film integration remain open, as do film-aware cost measurement,
EB interaction, physical validation and production scaling.
