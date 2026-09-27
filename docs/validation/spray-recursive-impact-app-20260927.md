# Recursive EB impact/film application input — 2026-09-27

Implementation `e425250cc012b07a63476b7689fd6e3307becc4a`, based on `19df4039`,
local main, version 0.246.0.

## Change

The MPI EB spray application now routes `time_integrator='sdc'`,
`subcycling=.true.` and `wall_impacts=.true.` to the coupled recursive
particle/film driver. Existing `walls` controls select film thermal evolution.
Film state, splash IDs and impact-model rejection diagnostics pass through the
same outer transaction. The obsolete blanket rejection of recursive films is
removed; periodic recursive migration remains unsupported. Required impact
viscosity, surface tension and wall-option validation remain unchanged.

The application regression removes the obsolete film-rejection expectation and
adds a dedicated two/four-rank impact case. The case uses one root interval and
two child intervals, gas chemistry and molecular transport, particle transfer,
EB splash and newborn-film evolution. Documentation now describes the enabled
input and keeps film-bearing recursive restart outside the checked scope.

## Independently completed checks

Artifacts: `build/spray-recursive-impact-app-checks/`. Case outputs:
`build/spray-amr-reactive-debug/spray-sdc-recursive-impact-app/impact_38_iz0va/`.
Existing fixed-mechanism Ninja Debug configuration; GNU Fortran 13.3.0 and
OpenMPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_impact_app$' --output-on-failure
python3 -m py_compile tests/application/check_recursive_eb_spray_impacts.py tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 CTest passed**, executing both two- and four-rank application runs, 37.71 s
test time and 37.75 s total. Provenance, build, Python syntax, project contract
and whitespace checks passed. No failures or skips; existing compiler warnings
remain. Final outputs are `build-final.log`, `ctest-final.log` and `LastTest.log`.
An earlier run also passed before the explicit empirical-model rejection message
was added to the recursive branch; it is retained separately.

Each rank layout produces **four live splash particles and one film**, with
unique IDs, positive mass and closed component inventories. Gas, particle and
film CSV values agree within `3e-12` relative tolerance between layouts.
Combined budget residual is **1.3299325971753836e-15** on both layouts against
the unchanged `3e-12` criterion. Near-zero momentum entries use the existing
mass/energy reference scale. The application's absolute composite error is
`1.0408340855860843e-17`.

This is a one-interval fixed-mechanism application check. The existing broader
recursive application regression was syntax-checked but not rerun here. No
film-bearing recursive restart, regrid/rebalance, long-time evolution, selected
mechanism, physical validation or full external-solver parity was run.
