# Compiled pinned PPM normal projection — 2026-09-26

Comparison implementation `435b8dbe`, parent `e5d88e7f`, local main, version
0.246.0. No solver arithmetic, tolerances or golden data changed.

`tools/check_pinned_ppm_normal.py` reads PPM.H and PPM.cpp at manifest-pinned
PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` with `git show`.
It compiles unchanged scalar reconstruction/profile bodies and the complete
species/left/right characteristic projection block. A face-array shim and NASA7
SI adapter supply the environment; this is not a full AMReX/PeleC executable.
The Fortran stream calls the production cell reconstruction used by 3D normal
prediction. Thermodynamic coefficients come from the fixed seven-species database.

## Results and limits

**252 cases / 6,048 primitive-state values passed**, maximum scaled error
**1.1102230246251565e-16**, fixed limit `3e-12`. Each observable uses
`abs(reference-actual)/max(1,abs(reference),abs(actual))`.
Cases combine seven velocities (-2500, -100, -1e-8, 0, 1e-8, 100, 2500 m/s),
three supplied flattening factors (0, .37, 1), constant/linear/material-contact/
random stencils, and three dt/dx values (0, 1e-5, 8e-5). Random seed: 20260926.
Both face density, three velocities, pressure and seven species weights are
compared. Internal energy is not an observable; its profile scratch values do
not enter the pinned density/velocity/pressure result, and the reference body
recovers final internal energy separately through EOS.

The original normal/CTU unit passed **1/1**, 6.27 s test / 6.31 s total.
Provenance, Python syntax, project contract and whitespace checks passed.
The first oracle build failed because stencil index constants followed their
use; the shim declaration order was corrected, retaining the initial log.
No numerical comparison failures or skipped cases occurred.

GNU Fortran/C++ 13.3.0, Python 3.12.3, existing Ninja Debug tests-enabled
`build/spray-amr-reactive-debug`; artifacts `build/ppm-pinned-normal-checks/`:
`build.log`, initial `comparison.log`, successful `comparison-fixed.log`,
`unit.log`, and `oracle/` (pinned sources, driver, input/thermo/output tables,
compiler command, executable/source hashes and report).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_normal_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_normal.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_normal_3d --work-dir build/ppm-pinned-normal-checks/oracle
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_normal_3d$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_normal.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This closes a bounded normal-projection comparison, not full CTU parity.
WENO, multidirectional flattening calculation, sanitization/fallback extremes,
physical boundaries, transverse correction, coupled evolution, installed audit
and physical validation are outside this check. Complete PeleC parity is unproven.
