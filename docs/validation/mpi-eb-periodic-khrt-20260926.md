# Periodic MPI EB KH children — 2026-09-26

Implementation `9d22aa3b`, parent `659a9d70`, local main, version 0.246.0.
Periodic EB particle stepping now permits KH child generation. Before release,
a private copy of each parent follows the existing segmented drift across
periodic seams. A new optional per-parcel wall-contact output distinguishes
EB/domain-wall reflection from periodic crossing or outflow. Contacts suppress
release for that substep; accepted parents/children use the normal drift,
feedback and owner migration. Failed drift leaves the contact output unchanged.
Periodic film combinations remain rejected. No particle heat model is added.

## Focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, existing fixed
H2/O2 tests-enabled `build/spray-amr-reactive-debug`.

The drift CTest passed **1/1**, total 0.05 s. It checks contact output for EB
reflection, three periodic wraps without contact, contact after a periodic wrap,
and rollback when the event budget is exhausted, alongside existing trajectory
and impulse/energy checks.

The direct application checker passed **3 successful runs and 1 expected
context rejection**: uninterrupted MPI4, MPI2 checkpoint and MPI4 continuation.
The CTU fixture starts at x=0.003999999 m moving at 1000 m/s toward a periodic
x face, includes three regrid stages and measured-rebalance controls, and retains
KH lineage/pending state across rank changes. Final population is 66, all with
x=1.9999e-5 m inside the wrapped domain and positive x velocity. Gas, particles
and budgets agree across continuation within the unchanged 3e-12 bound.
No unexpected failures or skips occurred. The application CTest entry was added;
the direct checker is not counted as an additional CTest execution.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_drift pelef_mpi_eb_spray_3d --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_drift$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-periodic-khrt-checks/application --khrt --periodic
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Build warnings remain.
No tolerances or golden values changed. Artifacts under
`build/mpi-eb-periodic-khrt-checks`: `build.log`, `ctest.log`, `LastTest.log`,
`application.log`, `application/impact_2gvn94ij`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

The eligibility preview adds an extra drift evaluation for periodic KH parents;
scaling/performance was not measured. No selected-mechanism, full matrix,
installed, external-solver or physical breakup qualification ran. Complete
PeleC parity remains unproven.
