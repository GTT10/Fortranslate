# Molecular transport in periodic SDC — 2026-09-26

Implementation `88c15957`, parent `d35f8258`, local main, version 0.246.0.

The built-in periodic SDC source evaluator now accepts molecular transport as
a second source family alongside external force. It evaluates existing periodic
transport face fluxes and their divergence at the old state and each candidate
new state. Viscosity, conduction, species diffusion and barodiffusion flags are
forwarded. Source prediction and final assembly average old/new rates through
the existing SDC loop. Minimum transport theta is accumulated across evaluations.
The application retains its existing transport timestep restriction and now
permits nonreacting periodic SDC with transport. Custom callbacks remain an
alternative to built-in sources. Example: `cases/reactive_boundaries_3d/sdc_transport.nml`.

## Completed checks

- Fixed tests-enabled and selected tests-disabled Debug builds passed in
  `build/spray-amr-reactive-debug` and `build/spray-amr-reactive-selected-debug`.
  GNU Fortran 13.3.0, Python 3.12.3, Ninja, pinned H2/O2 bundle; warnings remain.
- **4 application runs passed**: fixed/selected × transport enabled/disabled.
  A 4³ periodic hotspot in a 0.004-m cube uses 900-K base temperature, 250-K
  hotspot rise, 300000 Pa, dt=1e-7, three SDC iterations and external forcing.
  Transport changes maximum pointwise temperature by **0.004942116 K** (fixed)
  and **0.002495693 K** (selected pure N2), exceeding the inactive-path gate 1e-6 K.
  Normalized all-component inventory residuals are **2.62535e-17** and
  **1.28977e-17**; minimum transport theta is 1 in both enabled cases.
- `unit_reactive_ppm_step_3d`: **1/1 passed**, 5.29 s test / 5.34 s total.
- Initial fixture used unsupported problem name `hotspot`; the parser rejected
  it. Correcting it to `reactive_hotspot` produced the runs above. Initial failure
  log is retained. No final failures or skips.
- Provenance, Python syntax, contract and whitespace passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d test_reactive_ppm_step_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2
python3 build/sdc-transport-checks/check.py
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py build/sdc-transport-checks/check.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-transport-checks/` holds checker, inputs, CSVs, build/run/CTest
logs, provenance and report. The existing PeleF transport discretization and
limiting are reused; these checks establish routing, active transport and
inventory, not pinned full-solver parity or convergence. Chemistry, physical
boundaries, AMR/MPI/EB SDC and physical validation remain open.
