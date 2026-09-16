# MPI film evolution in coupled splitting — 2026-09-17

Source: parent `d3854e5` plus this change. Source hashes and logs:
`build/spray-mpi-film-split-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`advance_mpi_mixture_reaction` accepts an optional film array, component inventory
and film settings together. Physical boundary labels are required with film state.
Film mode is agreed collectively before mode-dependent calls. Private film copies
advance by half an interval before particle/chemistry/flow splitting and by another
half after it, using the existing owner-routed heating/evaporation API. Film wall
heat joins the global boundary ledger once; film substeps join the reported count
with overflow checks. Only whole-interval success commits gas, temperature,
particles, film state, component inventories and optional particle ID watermark.

This evolves supplied stationary films. Particle impact still needs to create
films inside the coupled interval; mixed deposition/splash integration and
application film controls remain open. The default no-film path is unchanged.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests passed 2/2 in 16.97 s, with no failures or skips. Provenance,
contract and whitespace checks passed.

The new film/particle/chemistry case starts with one 1e-9 kg, 330 K mixed film
on rank zero and no free particles. The 400 K wall heats the film; it also loses
mass through evaporation. Mass, momentum and energy including wall input agree
within 2e-12 of initial inventories, and elemental inventories within 2e-10.
A deliberately invalid reaction on one rank fails after the first film half-step
and exactly preserves original film temperatures/masses, composition, gas and
temperature, with zero output ledger and substeps. Gas hydro/transport/SGS are
not enabled in this new case; their no-film split checks remain in the same entries.

The first new test used its host reaction array before initialization and crashed.
Loading its own mechanism fixes the fixture; the failure is retained in
`ctest-uninitialized-fixture.log`. No tolerance was relaxed. No Python files changed.
No full configuration matrix or physical film validation is claimed.
