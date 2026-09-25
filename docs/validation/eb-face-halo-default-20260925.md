# Face-only EB halo callers — 2026-09-25

Code `9fd9154`, parent `eb3194e`. The legacy face-only halo fixtures explicitly
request `difmag=0`, preserving their raw Riemann donor comparison, unused NaN
slots, conservation and transactional rollback checks. The ghost fixture also
checks that omitted viscosity rejects an external halo without volume fractions
and returns zero flux. Production callers were inspected; integrated EB paths
already forward complete halos for positive viscosity. No solver code,
tolerances or golden hashes changed.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. **3/3 CTests passed**, no failures or skips:
EB hydro regression 0.93 s, ghost unit 0.19 s, patch halo unit 0.07 s
(1.23 s total). This is a focused caller compatibility check, not a full
regression matrix or external solver/physical validation. Provenance, Python
syntax, project contract and whitespace checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_reactive_ghost_3d test_eb_patch_halo_3d test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '(eb_reactive_ghost_3d|eb_patch_halo_3d|reactive_eb_hydro_3d)' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-face-halo-default-checks/{build.log,LastTest.log}`.
