# EB tree Mach inlet application and restart — 2026-09-25

Implementation `e268525`, base `5e19201`.
The shared gas application driver now accepts six-face
`inflow_mach`, `inflow_temperature` and `inflow_pressure` in `&eb_tree` and
forwards active controls to serial/MPI common or recursive full stepping and
automatic timestep selection. Fixed and selected mechanism front ends share
the driver. Omitted controls retain the prior optional-argument path.

Finite Mach values 0..100 enable an inlet; negative values retain outflow.
Active temperature must be positive; pressure is zero for extrapolation or
positive up to 1e15 Pa. Inactive temperatures/pressures must be zero. Runtime
EOS reconstruction still determines thermodynamic validity. Active controls
prefix application checkpoints with `PELEF_EB_TREE_MACH_INLET_V1` and 18
round-trip real values; all values must match on continuation. Old no-inlet
headers remain readable when no inlet is enabled. The MPI EB spray front end
explicitly rejects active inlets until its coupling path is connected.

## Focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug,
Python 3.12.3. Fixed serial/MPI executables and selected H2/O2 serial/MPI
executables built in their existing separate build directories.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_mpi_eb_tree_mach_3d$' --output-on-failure --output-log build/eb-mach-application-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The fixed five-patch plane-cut application fixture enables automatic timestep
selection, recursive stepping, chemistry and transport. Low-x Mach 0.2 uses
400 K and extrapolated pressure; high-y Mach 0.1 uses 400 K and 120000 Pa.
The checker compares exact serial/2/3-rank final CSVs, same-rank continuation,
2-to-3-rank continuation and serial/MPI checkpoint interchange. It retains
existing conservation, checkpoint-corruption, mode-mismatch and output-protection
checks and adds changed Mach/temperature/pressure and invalid inlet input cases.
No tolerance was relaxed.

**1/1 passed**, test **93.66 s**, total **93.70 s**. Maximum observed budget
residual was **3.01067808737867141e-15**. Changed inlet logs specifically report
`Restart inlet mismatch`. Both builds, provenance, Python syntax, project
contract and whitespace passed; no failures or skips. Compiler warnings remain.
Logs: `build/eb-mach-application-checks/{build-fixed.log,build-selected.log,
ctest.log,LastTest.log}`. Persistent fixtures and per-run logs:
`build/spray-amr-reactive-debug/tests/mpi_eb_tree_mach_3d/`.

Selected execution, coupled spray inlet execution, full matrix, production
scaling, reference solver comparison and physical nozzle validation were not
performed for this increment.
