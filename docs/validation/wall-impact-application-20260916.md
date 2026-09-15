# Wall-impact application and checkpoint — 2026-09-16

Source: parent `0becacc` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`, fixed H2/O2 mechanism, MPI disabled.

## Implemented behavior

The CLI enables selected stationary wall faces, writes film CSV/history state,
and persists films plus wall-only child lineage in an opt-in checkpoint schema.
Restart validates film material, temperature, mass/volume, position and enabled
wall attachment before checking total inventories. Existing film-free layouts
remain. Per-port surface tension/viscosity can be used with wall impact alone.
The writer refuses nonempty film arrays unless their schema is enabled.
Static patches do not yet flow, exchange heat or evaporate.

## Focused evidence

Four CTest entries passed in 20.56 s, zero final failures/skips: checkpoint,
materials, wall coupling and application CLI. CLI contains 20 tests and completed
in 19.57 s CTest wall time. The new CLI fixture checks low-speed deposition and
high-speed four-child splash, film mass in history, lineage and exact gas,
parcel, film CSV plus final-history equality after restart. Changed contact
angle, negative film height, a non-wall impact face and zero surface tension
reject. The integration test additionally checks transactional wrong-context
film reads, exact film mass/height roundtrip and wall-only lineage restoration.

Initial compile attempts exposed an overlong declaration and a misplaced optional
argument declaration; both were corrected. The first CLI run rejected the new
film-output path constructor because its character lengths differed; an explicit
constructor length fixed it. No tolerances or golden outputs were changed.

The supplied example completed 10 steps to 1e-6 s: five parcel records (one spent
parent plus four live children), one film patch. Deposited mass was
2.6392009464910826e-11 kg; airborne mass 7.360799053508917e-11 kg. Maximum relative
system mass/energy errors were 1.3573132322374105e-16 and
3.083373497639191e-16. This isolates bookkeeping and is not physical wall-impact
qualification; the declared wet residual-deposition closure remains.

Provenance, Python syntax, project contract and whitespace checks passed.
Artifacts: `build/wall-splash-example/` includes CSVs, checkpoint, run log and
`result.json`; build output is `build/port-spray-debug/film-app-build.log` and
CTest details are in `build/port-spray-debug/Testing/Temporary/LastTest.log`.
CLI tests otherwise use temporary directories. No broad matrix was run.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(application_cli|integration_(checkpoint|materials|wall_coupling))$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
mkdir -p build/wall-splash-example
cd build/wall-splash-example
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/wall_splash.nml
```

Use a fresh output directory for the example. Source hashes and the wet-impact
closure are recorded in `references/spray_wall_impact.json`.
