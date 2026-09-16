# MPI port-specific injection profiles — 2026-09-17

Source: parent `84d34d5` plus this change; hashes and build/check logs in
`build/spray-mpi-port-profile-checks/`. GNU Fortran 13.3.0 Debug and Open MPI
4.1.6, using the existing fixed-mechanism MPI build and serial Debug build.

Each of 16 ports can now supply up to 64 piecewise-linear flow samples. A port
profile overrides the shared table; otherwise shared flow or fixed port mass
applies. Validation requires coverage only of intervals using the selected table.
Zero-dose events consume no IDs. Every supplied table enters restart identity.
The common checkpoint context capacity increases from 16384 to 131072 characters
without changing existing context strings or checkpoint schemas.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_[24]$' --output-on-failure
cmake --build build/port-spray-debug --target test_spray_integration --parallel 4
ctest --test-dir build/port-spray-debug -R '^spray_integration_checkpoint$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

MPI application tests: 2/2 passed, 37.99 s, no failures/skips. Shared flow supplies
port 1 with two 7.5e-11 kg doses; port 2's own profile gives zero then 1.5e-10 kg.
Total injected mass is 3e-10 kg; parcel counts are 2/4/7/7 through 2e-8 s and
corrected mass/energy errors remain below 2e-11. At a simultaneous-event checkpoint,
2↔4-rank continuation matches continuous gas and particles (relative tolerance
2e-11, absolute 1e-25). Changed profile identity and incomplete time coverage are
rejected before history output.

A 16-port × 64-point zero-flow case writes a context longer than 16384 characters
and reloads on changed ranks with all 16 events consumed and ID watermark zero.
This covers table capacity/persistence, not 16-port physical interaction.
Existing scalar, shared-profile and multiport fixed-dose cases also pass.
Serial checkpoint test: 1/1 passed, 0.06 s total; no failures/skips. Provenance,
project contract, Python syntax and whitespace checks pass.

Application artifacts:

- `build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_2/run-v7qe814g/`
- `build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_4/run-xhdd7oad/`

No selected-mechanism profile rerun, full matrix, physical injector validation or
production scaling is claimed. MPI finite-nozzle/size-distribution CLI controls,
breakup, wall films, SGS and particle AMR/EB integration remain open.
