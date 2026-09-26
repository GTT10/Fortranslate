# Recursive spray SDC boundary forwarding — 2026-09-27

## Source and scope

Implementation: `e51d8a1a`, based on `912bc1fc`, local main, version 0.246.0.
The recursive wrapper forwards inflow state/Mach/temperature/pressure to gas and
particle sampling, and wall temperatures to gas SDC. Population drift uses
outflow for inflow faces and no-slip reflection for isothermal walls, matching
the existing common-step wrapper. Gas halo sampling retains the original faces.

This is a library boundary integration check, not application dispatch, particle
wall-impact qualification, external solver parity, or physical spray validation.
Particles remain inside the fine patch. Subtree lifecycle integration, breakup,
films and particle forcing remain open.

## Independently completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, existing Ninja Debug build
`build/spray-amr-reactive-debug`; fixed H2O2 mechanism. Artifacts:
`build/spray-recursive-boundary-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_(inflow|isothermal)_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, final build, syntax, contract and whitespace checks passed.
CTest: **2/2 passed**, no skips or test failures, 15.09 s total:

| Case | Time | Relative composite inventory residual |
| --- | ---: | ---: |
| inflow | 7.62 s | 1.3919441503977229e-15 |
| isothermal | 7.42 s | 1.9455220335484211e-15 |

Both retain the existing 3e-12 inventory tolerance, check accepted level counts
[1,2], feedback/heating/evaporation, and rejected-hydro whole-step rollback.
Two preceding compile attempts failed: an accidental callback signature edit,
then an overlong test call line. Both were corrected before running tests;
`build.log`, `build-final.log`, and successful `build-retry.log` preserve evidence.
`ctest.log` and `LastTest.log` record the completed tests.
