# Conserved-to-primitive source conversion — 2026-09-26

Implementation `367004eb`, parent `f904758f`, local main, version 0.246.0. Implementation adds public
`reactive_pelec_source_to_primitive` in `reactive_1d_mod`.

The routine implements pinned `Hydro.H::pc_srctoprim`, including the reference's
full density-source times velocity-squared term in the internal-energy source.
Velocity and species sources subtract density dilution. Pressure conversion uses
ideal-mixture derivatives `dp/de|rho = rho R/cv` and `dp/drho|e = R T`, with NASA7
properties at the accepted primitive state. This preserves the reference policy;
it does not add composition derivatives absent from that policy. Invalid state,
shape or nonfinite source returns failure and zero outputs.

## Verification

`unit_reactive_pelec_riemann` passed **1/1**, 0.00 s test / 0.05 s total.
The new fixture uses rho=2, velocity=(3,-4,12), pressure=300000 and pure N2.
Conserved sources rho=0.2, momentum=(5,-6,8), total energy=1000 produce velocity
sources (2.2,-2.6,2.8) and pinned internal-energy source 898.8. Species dilution,
pressure derivatives, zero-source identity and NaN rollback are checked.
No failures or skips. Provenance, Python syntax, project contract and whitespace
checks passed. GNU Fortran 13.3.0, Python 3.12.3, existing tests-enabled Ninja
Debug build; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_pelec_riemann --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_pelec_riemann$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-source-conversion-checks/{provenance.log,build.log,ctest.log,unit.log}`.
Reference inspected with `git show` in `/mnt/sn850/Projects/gtt2/external/PeleC`
at `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

This is a source-grounded focused check, not a compiled external conversion
comparison. Grid-level conversion and application source generation/integration
remain open; physical validation and complete PeleC equivalence are unproven.
