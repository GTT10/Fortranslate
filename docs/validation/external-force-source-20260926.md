# External force source construction — 2026-09-26

Implementation `27615093`, parent `968ffbd4`, local main, version 0.246.0.

`compute_reactive_external_source_3d` implements the regular-cell base forcing
in pinned `Source/External.cpp::fill_ext_source`. Input force is force per volume,
not acceleration: momentum sources equal its three components and total-energy
source is velocity dot force. Density/species sources are zero. The caller
supplies the state time level and any halo cells. Invalid dimensions, nonfinite
inputs or nonpositive density return failure and zero sources transactionally.

The pinned old/new source constructors evaluate this same function at different
state time levels. `Advance.cpp` reconstructs the candidate and recomputes new
sources before the final state update. Therefore this constructor alone does
not establish the required temporal work integration or SDC parity.

## Verification

`unit_reactive_ppm_step_3d` passed **1/1**, 5.00 s test / 5.05 s total.
New checks use nonuniform density, velocity (20,10,-5) and force (2,-3,4):
momentum sources are density-independent and work is -10. Mass/species sources
vanish; a bad density at the last cell leaves the entire output zero. Existing
source/split/CTU checks also ran. No failures or skips. Provenance, Python
syntax, project contract and whitespace passed. GNU Fortran 13.3.0, Python
3.12.3, existing tests-enabled Ninja Debug build; compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/external-force-source-checks/{provenance.log,build.log,ctest.log,unit.log}`.
Reference inspected via `git show` at
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` in
`/mnt/sn850/Projects/gtt2/external/PeleC`.

Application input, old/new work correction, SDC integration, MPI/EB and physical
validation remain open. No application forcing support is claimed here.
