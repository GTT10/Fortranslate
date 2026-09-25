# MPI EB spray thermal gas walls — 2026-09-26

Implementation `9b1eb78`, parent `b6672bc`, local main, version 0.246.0.
The coupled MPI EB spray API now forwards optional wall temperatures to the
full gas split. Isothermal gas labels map to stationary specular particle walls.
Both particle half-steps and the gas split remain within one transaction.
The shared application forwards tree wall temperatures and permits isothermal
input; explicit-boundary wall films and KH child generation remain rejected.
The existing gas checkpoint context binds all six temperatures.

## Completed focused evidence

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, existing fixed
H2/O2 tests-enabled `build/spray-amr-reactive-debug`. No backend changes.

Two four-rank coupled tests passed: common 2.31 s, recursive 4.62 s, CTest total
6.98 s. A two-level planar EB fixture with 2000 K x walls, implicit chemistry,
CTU and particles checks combined inventories within the unchanged 3e-12 bound,
nonzero chemistry, heating/evaporation, outer-wall and EB reflection, and full
state/particle/ledger rollback after private particle work and invalid hydro.

The application checker passed **4 successful executions and 11 expected
rejections**. It compares fixed-mechanism MPI2/MPI4 CTU gas/particle outputs
within 3e-12 and MPI2 checkpoint to MPI4 continuation against uninterrupted MPI4.
Continued gas and sorted particle rows match exactly; budget rows meet 3e-12.
The fixture uses 2000 K x/z walls, y outflow, automatic timestep selection,
artificial viscosity .037 and particles crossing the upper z wall.
A changed wall temperature rejects with `Restart wall temperature mismatch` and
no output. Existing boundary/method/coefficient/context/corruption/overwrite
checks also pass. The new application CTest entry was registered; the checker
was executed directly and is not counted as an additional CTest run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_reactive_ctu-isothermal(-subcycle)?$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-isothermal-spray-checks/application --isothermal --restart --difmag .037
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace passed; no unexpected failures or
skips. No tolerances or golden values changed. Build warnings remain.
Artifacts: `build/mpi-eb-isothermal-spray-checks/{build.log,ctest.log,LastTest.log,application.log}`
and `application/eb_spray_fshwjonc` below that directory.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

## Limits

No direct particle-wall heat-transfer model is added: particles reflect, while
gas conduction supplies the thermal boundary transfer. Selected spray, thermal
regridding/rebalancing, external solver parity, full matrix and physical wall-flow
validation were not run. Complete PeleC parity remains unproven.
