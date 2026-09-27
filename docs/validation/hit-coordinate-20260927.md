# HIT coordinate conversion — 2026-09-27

Implementation: `90085c9b2908540b24a4ca25e7981b9fe6e2e2ea`, based on
`5caa1182`. The optional `reactive_hit.coordinate_scale` converts file
coordinates before interpolation. Its default is 1; nonpositive/nonfinite
values and multiplication overflow are rejected. Velocity scaling is unchanged.

The gamma-law example now uses the reference's centimetre lengths converted
to metres: 2π cm domain, 0.5 cm Taylor length and file coordinate factor 0.01.

## Completed checks

GNU Fortran 13.3.0 Debug, OpenMPI 4.1.6, SUNDIALS 7.2.0; existing unchanged
`build/amr-cvode-app-debug` configuration with MPI, spray and SUNDIALS enabled.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/amr-cvode-app-debug -R '^application_gamma_law$' --output-on-failure
python3 -m py_compile tests/application/check_gamma_law_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All passed. CTest: 1/1, zero failures/skips, 1.10 s test / 1.15 s total.
Seven successful application runs and five expected input rejections.
The new centimetre-coordinate input with factor 0.01 produced byte-identical
CSV to the metre-coordinate input with default factor 1. Zero factor was
rejected without publishing output. Existing gamma-law state/closure checks
remain active; maximum EOS residual was `5.7411274570087557e-16`.

Build and preserved CTest logs: `build/hit-coordinate-checks/`.
Inputs, outputs and individual logs:
`build/amr-cvode-app-debug/gamma-law-app/gamma_1570bmrj/`.

This is a manufactured unit-conversion/application check. The 16³ case was
not run, and no full PeleC trajectory or physical turbulence validation was
performed. Binary input uses the same post-read conversion but was not
separately rerun here.
