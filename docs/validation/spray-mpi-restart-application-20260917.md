# MPI spray application save/restart — 2026-09-17

Parent `26f4d49` plus this change. Build/test/preflight logs and source hashes:
`build/spray-mpi-restart-application-checks/`. Per-run namelists, logs, checkpoint,
histories and snapshots remain in fresh `run-*` directories under
`build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_{2,4}/`.

The MPI executable now accepts `checkpoint_path` and `restart_path`. It restores
gas, mixed particles, time/steps, budgets and the ID watermark, skips initial
injection and retains the original reference budget. Context identity binds model
and numerical settings while allowing different rank counts, stop times, step
limits and file paths. Incompatible input is rejected before history output.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6; existing matching MPI configuration:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

2/2 application entries pass (8.12 s), no failures/skips. Each compares a two-step
continuous run with save after one step and resume on the same rank count (exact
final gas and particle CSV numeric strings), then on a changed count (2-to-4 and
4-to-2, relative 2e-11/absolute 1e-25 field tolerance). Tests verify resumed step
numbers, no duplicate injection, four particle IDs, retained energy budget and
rejection of altered liquid heat capacity before creating history. Earlier input,
output, conservation and overwrite checks remain active. Provenance, project
contract, Python syntax and whitespace checks pass. Test resource declarations
now reserve four processes for each entry because both launch four-rank variants.

These are bounded two-step synthetic-mixture checks. They do not establish
changed-rank bitwise parity, long-run physical/scaling qualification, scheduled
injection restart, or a full compiler/backend matrix.
