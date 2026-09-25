# Frozen sources through full splitting — 2026-09-26

Implementation `9312c538`, parent `df776f2d`, local main, version 0.246.0.

`advance_reactive_full_3d` accepts conserved source rates with four supplied
halo layers. Rates are held fixed during the step and routed to the single
PPM/PeleC hydro update between transport half steps. Periodic hydro receives
physical cells and wraps its own source halos; physical boundaries use the
supplied halos. Shape, finiteness and reconstruction compatibility are validated
before reaction or transport work. Existing candidate-state rollback and
boundary-only inventory semantics are retained.

## Verification

`unit_reactive_ppm_step_3d` passed **1/1**, 4.94 s test / 4.98 s total.
New periodic and outflow comparisons require exact state/temperature agreement
between the full dispatcher and explicit transport-half / source-hydro /
transport-half composition. A source failure after the first transport stage
checks complete state/temperature rollback and zero boundary change. Existing
inventory, source-step and hydro tests also ran. No failures or skips.
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

Artifacts: `build/ppm-source-split-checks/{provenance.log,build.log,ctest.log,unit.log}`.
Chemistry was disabled in the new source fixtures. This establishes source
routing through the library composition, not reacting-source qualification,
application source generation, MPI exchange, time-dependent quadrature/SDC or
physical validation.
