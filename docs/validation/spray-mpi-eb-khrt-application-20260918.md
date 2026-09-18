# MPI EB reactive/application KH children — 2026-09-18

Implementation `a70cf50`, base `27c2411`.

The reactive wrapper accepts a shared watermark without films for specular KH
flight. KH mode requires a watermark and retains it privately through both
particle halves and the gas update. The MPI application removes its KH-mode
rejection, passes shared IDs on both wall paths, and accepts four-child release
groups when validating restart IDs. Existing particle payloads already retain
lineage, deformation, RT time, initial count and pending shed volume.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja; reused fixed tests-enabled and selected H2O2 tests-disabled builds.
The reactive fixture retains bounds/runtime checks and invalid/zero/overflow traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d test_spray_mpi_eb_impact_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_reactive_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-khrt-app-checks/runs --khrt
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Reactive CTest **2/2 passed**, total **3.22 s** (two ranks 1.57 s; four ranks
1.61 s). The application driver passed four launches: full four-rank,
stopped two-rank, resumed four-rank and deliberate context mismatch. No unexpected
failures or skips. Fixed/selected executable builds, provenance, syntax, contract
and whitespace passed. Selected execution was not run. The application driver is
registered as `spray_mpi_eb_khrt_application`; it was invoked directly here.

The reactive fixture adds a 1000 m/s mixed-liquid KH case without films. A gas
solver failure after private first-half release preserves the incoming particle
population/composition, gas, watermark and diagnostic sentinels. Success produces
children and closes composite mass/species/momentum/energy to relative 3e-12,
scale floor 1e-30. Existing reactive impact/film checks remain in the same test.

The application starts two water parcels at 1000 m/s with 1e-16 kg total dose,
reference diameter cubed 1e-25 m3, and drag/heating/evaporation disabled to isolate
breakup energy exchange. Implicit gas chemistry, molecular transport and HLLC
hydro remain enabled. Prescribed regrids move refinement, coarsen and refine.
The step-two V1 checkpoint contains 34 parcels and watermark 34, including child
lineage and nonzero KH initial counts. The final population is 66 with IDs 1–66.
Gas, particles and cumulative budget columns agree with changed-rank continuation
within relative 3e-12, scale floor 1e-30. A changed normalized wall setting rejects
context even with impact films disabled. Rebalance controls are enabled, but an
accepted ownership move is not an asserted gate of this fixture.

Artifacts: `build/spray-mpi-eb-khrt-app-checks/{build.log,build-test.log,ctest.log,application.log,build-selected.log,provenance.log,contract.log}`;
application inputs, outputs and checkpoints under `runs/impact_vfyqp831` there.
This remains bounded source-first KH/EB coupling. No physical breakup validation,
broad matrix, installed audit, external solver parity or scaling check was run.
