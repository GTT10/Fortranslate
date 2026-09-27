# AMR application CVODE SDC — 2026-09-27

## Implementation

Source: `e97198002e71e9820931860833616f71882f5e69`, parent
`185c0fe480c4d44a7e78ec0a083738148229238d`. A fresh configuration built the
implementation working tree before commit; the commit contains those same
implementation/test files. This is not a separate post-commit checkout audit.

Fixed and selected serial/MPI AMR tree frontends now accept
`sdc_reactor='cvode'` with SUNDIALS. They pass the forced CVODE callback through
the optional LES wrapper, application driver and common/recursive SDC stages.
Native remains the default. The driver rejects disagreement between the named
selector and callback presence. Fixed and selected binaries link their own
mechanism-compatible CVODE runtime; their core libraries are not mixed.

Existing full-configuration checkpoint identity already includes `sdc_reactor`.
No checkpoint schema change was necessary. Changed reactor selection rejects
restart before final CSV publication. The input reader continues to reject
CVODE when the build lacks SUNDIALS or chemistry/SDC is disabled.

## Focused checks

New build: `build/amr-cvode-app-debug`, GNU Fortran 13.3.0, Open MPI 4.1.6,
SUNDIALS 7.2.0 official static Fortran interfaces, Debug, MPI/spray/SUNDIALS
ON, selected full H2/O2 bundle. Build-time MPI dependency audit passed ABI 40.
Compiler warnings include existing real comparisons/conversions and
internal-procedure executable-stack linker warnings; no build errors.

| Test | Scope | Result |
| --- | --- | --- |
| `application_amr_sdc_cvode` | Fixed mechanism, common step | Passed, 47.06 s |
| `application_amr_sdc_cvode_recursive_selected` | Selected mechanism, recursive step | Passed, 37.13 s |

**2/2 passed, zero failures/skips; total 84.25 s.** Each test uses a two-level
periodic tree, 4³ root cells and one refined parent cell. The active 1100 K
reactor heats the gas; uniform external forcing and molecular transport are
included. Each checks exact final gas/budget CSV bytes for serial versus
2-rank execution, 2-to-4-rank continuation, and MPI-to-serial continuation.
The checkpoint stops after step 1 of a multi-step trajectory. Restart with
`native` instead of `cvode` fails with the expected context diagnostic and
publishes no final CSV.

Maximum normalized Euler inventory residuals were
`4.9222869465612917e-18` (common) and `4.9222869207118976e-18` (recursive),
against the unchanged `3e-12` criterion. Species conversion is not counted
as an Euler conservation error.

The LES wrapper is built and exercised with LES disabled. Enabled LES plus
CVODE, EB/spray forced chemistry, scheduled/adaptive regridding, other
mechanisms, installed binaries, a full build matrix and external physical
validation were not tested here. Common selected and recursive fixed are
connected but not separately run in this focused pair.

## Reproduction and artifacts

From the repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/amr-cvode-app-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_SUNDIALS=ON -DSUNDIALS_DIR=/home/gtt13/projects/pelec/build/sundials-7.2-fortran/install/lib/cmake/sundials -DPELEF_MECHANISM_BUNDLE=/home/gtt13/projects/pelec/mechanisms/h2o2_full.json -DPELEF_MECHANISM_SOURCE=/home/gtt13/projects/pelec/mechanisms/h2o2_cantera.yaml
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_amr_tree_reactive_3d pelef_mpi_amr_tree_reactive_3d pelef_amr_tree_reactive_3d_selected pelef_mpi_amr_tree_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/amr-cvode-app-debug -R '^application_amr_sdc_cvode' --output-on-failure
python3 -m py_compile tests/application/check_amr_sdc_cvode_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All checks completed successfully. CMake was regenerated after registering the
tests; no implementation recompilation was needed. Pinned mechanism source SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/amr-cvode-app-checks/` (`configure.log`, `build.log`,
`reconfigure-build.log`, `preflight.log`, `contract.log`, `ctest.log`, `LastTest.log`).
Inputs, CSVs, checkpoints and per-run logs:

- `build/amr-cvode-app-debug/amr-sdc-cvode/amr_cvode_1193ucei/`
- `build/amr-cvode-app-debug/amr-sdc-cvode-selected/amr_cvode_3rjx24mp/`
