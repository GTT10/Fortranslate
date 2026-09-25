# MPI EB spray prescribed inflow — 2026-09-26

## Implementation

Source 2d6a2df, parent eefe747, local main. The MPI EB reactive-spray wrapper
forwards explicit boundary_faces and inflow_state to the gas full split.
Explicit faces accept only inflow/outflow and require particle_outflow=true;
validation is collective before particle stepping. Existing particles can
leave either open face. A gas inlet does not inject particles. Gas and both
particle halves retain a single transaction, including boundary inventories.

The shared MPI tree application now permits prescribed inflow with spray,
checks particle_outflow after reading spray input, and passes generated inlet
states through both ordinary and impact-model wrapper calls. The existing
prescribed-inlet gas envelope and spray context bind restart.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 via /usr/bin/mpiexec, Ninja Debug,
Python 3.12.3; existing tests-enabled build/spray-amr-reactive-debug.
No compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_reactive_ctu-inflow(-subcycle)?$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-spray-inflow-checks/application --inflow --restart
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Library: **2/2 passed**, 7.21 s total (common 2.43 s, recursive 4.74 s),
0 failed/skipped. Four ranks include two idle ranks. The fixture checks coupled
chemistry, transport, particle heating/evaporation, embedded-wall reflection,
composite inventory closure at 3e-12 relative tolerance, and full transaction
rollback after gas rejection or a one-rank invalid inlet density.

Application: **4 successful executions and 6 expected rejections**. Two- and
four-rank final gas/particle data agree; two-rank checkpoint continuation on
four ranks reproduces the uninterrupted result. Budget fields agree within
the existing 3e-12 relative tolerance. Expected failures cover overwrite,
unsupported raw restart, changed spray context, changed hydro, changed
prescribed velocity and truncated particle payload. No unexpected failures.

Provenance, syntax, project contract and diff checks passed. No tolerances or
golden data changed. Existing compiler warnings remain. Mechanism SHA:
0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Artifacts: build/mpi-eb-spray-inflow-checks/{build.log,build-app.log,LastTest.log}
and application/eb_spray_ijunid7y under that directory.

## Limits

This increment ran the fixed-mechanism application. Selected execution,
active impact/film coupling and actual outer particle escape were not rerun
with prescribed gas inlets. Periodic/wall MPI tree particle dispatch remains
open. No full matrix, installed audit, external solver parity or physical
spray validation was run. Complete PeleC parity remains unproven.
