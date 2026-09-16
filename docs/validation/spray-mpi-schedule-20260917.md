# MPI repeated cone schedule and restart — 2026-09-17

Parent `a2e2c0a` plus this change. Build/test/preflight logs and hashes:
`build/spray-mpi-schedule-checks/`. Per-case inputs, logs, histories, snapshots and
checkpoints remain in fresh application test `run-*` directories.

Repeated events use a start, interval and event count, with existing mass and
parcel count per event. The time step lands on each pending event. Injection,
migration and ID allocation use candidate state. The event count and accumulated
injected budgets survive checkpoint/restart. Context v2 binds the schedule and
changes the reference convention to initial gas plus separately accumulated
injection; v1 application checkpoints are therefore rejected.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; existing matching MPI configuration:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

2/2 application entries pass, 12.48 s, no failures/skips. New cases inject four
particles at 0, 7.5e-9 and 1.5e-8 s and stop at 2e-8 s. History counts are
4/8/12/12; the final watermark is 12 and injected mass is 3e-10 kg. Mass and
energy drift after accounting for injection remain below 2e-11. A split run
saves exactly at the second event, then resumes on changed ranks (2-to-4 and
4-to-2); final fields and particle records match continuous runs at relative
2e-11/absolute 1e-25. Event history proves no repeated second injection.
Existing application, overwrite and restart tests pass. Provenance, contract,
Python syntax and whitespace checks pass.

No independent multiport or tabulated-rate schedule, physical validation, full
configuration matrix or production scaling is claimed.
