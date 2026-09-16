# Independent MPI injection ports — 2026-09-17

Source: parent `4f75e75` plus this change. Source hashes and check logs are
in `build/spray-mpi-multiport-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 build; no compiler/backend configuration changes.

Up to 16 independent ports now supply timing, count, mass, position, direction,
speed, cone angle, diameter, temperature and two-component composition. A merged
schedule orders by time then port index. The saved global event ordinal restores
this ordering; port definitions enter checkpoint context. A shared rate profile
is integrated over each port's interval separately. Scalar-mode checkpoint
identity remains unchanged. The common liquid-property table remains shared.

Completed checks:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

Both application entries passed (2/2; 28.69 s; no failures/skips), including
existing scalar/profile/schedule/restart cases. The new two-port case produces
event counts 1/3/4/4 and parcel counts 2/7/10/10 through 2e-8 s. The simultaneous
event at 7.5e-9 s allocates IDs 3/4 to port 1 before 5/6/7 to port 2. Pure-component
inventories retain this assignment. Total injection is 6e-10 kg; corrected mass
and energy errors remain below 2e-11.

Checkpoint continuation after the simultaneous events matches continuous gas and
particles on changed ranks (2↔4), relative tolerance 2e-11 / absolute 1e-25.
Changed port mass rejects restart before history output. Shared-profile operation
produces parcel counts 0/5/8/8 and total injection 2.25e-10 kg, including an initial
zero-dose event. Provenance, project contract, Python syntax and diff checks pass.

Retained application outputs:

- `build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_2/run-v7jrhrq4/`
- `build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_4/run-qfacy2jf/`

No selected-mechanism multiport rerun, full configuration matrix, physical nozzle
validation or production scaling is claimed. Independent per-port rate tables,
MPI breakup/wall-film/SGS coupling and particle AMR/EB remain open.
