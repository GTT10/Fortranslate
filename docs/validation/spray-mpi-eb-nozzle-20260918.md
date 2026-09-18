# MPI EB finite outlets and size distributions — 2026-09-18

Implementation `2507c3f`, base `b546059`.

MPI EB spray ports now forward finite radius, nozzle model and size options to
the existing deterministic injection kernel. The normalized port settings enter
the existing restart context. Defaults retain the previous point and single-size
input behavior. No new distribution or nozzle formula was introduced.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja; reused fixed tests-enabled and selected H2O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_nozzle.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-nozzle-checks/runs
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The application driver passed all five MPI launches: complete two/four-rank,
two-rank stopped, four-rank resumed and deliberate negative-radius rejection.
No unexpected failures or skips. Both executable builds, provenance, Python
syntax, project contract and whitespace checks passed. Selected execution was
not run. The driver is registered as `spray_mpi_eb_nozzle_application`; this
check invoked it directly rather than through CTest.

Two ports use radius 0.1 mm, eight parcels per event and swirl 0.4 radians.
The first uses a solid outlet and mass-weighted bounded Weibull distribution;
the second a hollow outlet with spread 0.02 radians and number-weighted bounded
lognormal distribution. Each injects at time zero and 20 ns. Final-time records
allow direct checks before any gas exchange: outlet radius, hollow-ring location,
relative transverse swirl, temperature, 5–20 micrometre size bounds, nonuniform
sizes, equal represented masses or multiplicities, and total event mass 1e-10 kg
(with absolute error at most 1e-23 kg). Earlier records have thermal exchange.

Gas and parcel output agree across two/four ranks and after two-to-four-rank
restart to relative 3e-12, scale floor 1e-30. The application retains its coupled
inventory check. Invalid negative radii fail template creation before output.
The check does not newly qualify normal or chi-squared sampling, wet impacts,
physical nozzle flow, stochastic dispersion, scaling or experimental physics.

Artifacts: `build/spray-mpi-eb-nozzle-checks/{build.log,build-selected.log,application.log,provenance.log,contract.log}`;
inputs, CSVs, checkpoints and per-run logs under `runs/nozzle_fsxiuffo` there.
No broad matrix, installed audit or external solver parity was run.
