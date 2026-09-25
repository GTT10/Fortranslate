# Periodic CTU step with conserved sources — 2026-09-26

Implementation `e3bbee3e`, parent `dba63a29`, local main, version 0.246.0.

`advance_reactive_euler_ppm_ctu_3d` accepts a conserved source rate held fixed
over the step. The existing CTU flux API converts it for half-step face
prediction; the final conserved candidate receives `dt*source` exactly once.
The source option requires the PeleC pressure policy. Shape/finiteness checks
precede the zero-time identity path. Final EOS failure rolls back the complete
state, temperature and limiter diagnostic.

## Verification

`unit_reactive_ppm_step_3d` passed **1/1**, 3.42 s test / 3.47 s total.
New checks cover nonuniform periodic source inventory and species closure,
uniform density/momentum/energy/species injection with exactly one increment,
and complete rollback under an excessive negative energy source. Existing
multistep hydro and splitting checks also ran. No failures or skips.
Provenance, Python syntax, project contract and whitespace passed.
GNU Fortran 13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug build;
existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-source-step-checks/{provenance.log,build.log,ctest.log,unit.log}`.
The underlying flux source comparison is recorded in
[the 38-case compiled-reference check](ppm-conserved-source-20260926.md).

This establishes a library step for supplied frozen rates. It does not implement
time-dependent source quadrature or the reference SDC iteration. Physical-boundary
step routing, application source generation, MPI exchange and physical
validation remain open.
