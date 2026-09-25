# PPM double-transverse primitive sources — 2026-09-26

Implementation `5b6b2d27`, parent `34a65a62`, local main, version 0.246.0.

`reactive_ppm_transverse_pressure_3d` accepts optional primitive source rates,
an internal-energy-density source rate and a source interval. Pinned `pc_transdd`
adds these after conservative transverse conversion; pressure is floored after
its source. Velocity/species division uses pre-source density. Partial arguments,
invalid rates/intervals, single-transverse sources and nonpositive resulting
density fail with zero outputs. Existing callers retain zero-source behavior.

## Verification

Compiled unchanged `pc_transdo`/`pc_transdd` bodies passed **216 cases / 3,456
values**, maximum scaled error **1.4551886041586193e-16**, unchanged tolerance
`3e-12`. The original 144 cases remain; 72 double-transverse fixtures add finite
sources for every primitive component and internal energy, across three axes
and ordinary/signed-species/pressure-floor conditions. Tiny pressure floors
require exact agreement.

`unit_reactive_ppm_halo_3d` passed **1/1**, 6.40 s test / 6.45 s total, including
partial-source, wrong-stage and density-failure checks. No failures or skips.
Provenance, Python syntax, contract and whitespace passed. GNU Fortran/C++
13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug build; compiler warnings
remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_transverse.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_halo_3d --work-dir build/ppm-transverse-source-checks/oracle
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_transverse.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-transverse-source-checks/` contains provenance, build,
comparison, CTest and unit logs; `oracle/` contains pinned source, generated
driver, inputs, outputs and report with compiler and source/executable hashes.
Reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

This establishes the local correction kernel only. Full CTU source routing,
conserved-to-primitive source conversion, application integration and physical
validation remain open; no application or MPI source path is claimed.
