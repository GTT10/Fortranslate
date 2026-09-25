# Molecular transport at SDC outflow boundaries — 2026-09-26

Implementation `145a2e7d`, parent `58c52b0d`, local main, version 0.246.0.

The SDC built-in source evaluator now calls the existing physical-boundary
transport flux API for outflow configurations. It forms the cell divergence
from all faces, including index-zero faces, and retains old/new transport
boundary rates. Their time average is added to the boundary inventory and
removed from the volume-source inventory, preventing duplicate accounting.
The periodic transport route is unchanged. Application input now permits
outflow plus molecular transport. Other physical boundary types remain open.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `pelef_reactive_3d`, passed.
- **2 application cases passed**: nonuniform reacting hotspot with all-outflow
  faces and external force, transport enabled/disabled, 4³ cells and one step
  to 1e-9 s with two SDC iterations. Enabling transport changes maximum
  pointwise temperature by **6.588276278307603e-5 K**. The enabled run's
  normalized Euler inventory residual is **1.75720e-17** and transport theta
  is 1. Both runs pass the application's boundary/source-accounted gate.
- Provenance, Python syntax, project contract and whitespace passed. No
  failed or skipped runs; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2
python3 build/sdc-outflow-transport-checks/check.py
python3 -m py_compile tools/check_project_contract.py build/sdc-outflow-transport-checks/check.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-outflow-transport-checks/` contains checker, inputs, CSVs,
build/run logs, provenance and report. The selected build and CTest were not
rerun for this change. Prior outflow-transport input-rejection evidence is
superseded by this implementation. This checks active transport and aggregate
inventories, not convergence, full reference transport parity or physical
ignition validation. Prescribed pressure, inflow/walls and AMR/MPI SDC remain open.
