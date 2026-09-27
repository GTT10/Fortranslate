# Recursive periodic EB particle coupling — 2026-09-27

Implementation `cd916766aaa1303198c9fb3d376beeb2e9f99e46`, based on `d0242fa1`,
local main, version 0.246.0.

## Change

The recursive gas/particle wrapper accepts paired periodic faces and forwards
their axes through physical/copy initialization, subtree redistribution and
final collection. Existing stage halos, periodic drift, source transfer and
physical source correction now operate with those axes. Mismatched face pairs
reject collectively before publication.

The MPI EB spray application enables this path for particle-only recursive SDC.
Periodic film-copy placement remains unimplemented: supplying film state to the
periodic recursive library, or enabling periodic wall impacts in the application,
still rejects explicitly.

The application checker appends its requested periodic/domain-wall face override
after the default SDC outflow setting. Previously its ordering let that default
replace the requested condition. The new recursive periodic application test
therefore executes actual periodic faces.

## Independently completed checks

Artifacts: `build/spray-recursive-periodic-checks/`. Application cases:
`build/spray-amr-reactive-debug/spray-sdc-recursive-periodic-app/eb_spray_yoeq52d0/`.
Existing fixed-mechanism Ninja Debug build; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_periodic_2$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_periodic_app$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both selected CTests passed: the two-rank library test took **7.21 s**, and the
application checker **100.98 s** (101.03 s including CTest overhead). Provenance,
build, Python syntax, project contract and whitespace checks passed. No failures
or skips; existing compiler warnings remain. Outputs are `ctest.log`,
`ctest-app.log`, `LastTest-library.log`, and `LastTest-app.log`.

The library fixture advances a coarse physical particle across the right
periodic face into the left fine patch, with chemistry, transport, heating and
evaporation. Its combined inventory residual is **1.6251053807333461e-15**
against the unchanged `3e-12` criterion. The existing early-rejection rollback
and successful retry checks run in the same fixture.

The application runs three root intervals on both two and four ranks with
periodic x/z faces, EB reflection, gas forcing and particle acceleration. It
checks wrapped coordinates, gas/particle rank comparison, component inventories
and combined budgets using the existing `3e-12` criteria, plus excessive-force
and missing-restart rejection. Application comparison is tolerance-based,
not a claim of bitwise equality.

No periodic film evolution, accepted periodic restart, selected mechanism,
physical validation or full external-solver parity was run. Other boundary
combinations affected by the checker override ordering were not rerun here.
