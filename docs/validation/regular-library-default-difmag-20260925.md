# Regular hydro API default viscosity — 2026-09-25

Code `fb2820f`, parent `1d9561c`. PPM periodic and external-halo flux APIs now
resolve omitted `difmag` to 0.1 independently of pressure tracing or solver.
The regular boundary hydro and full-composition dispatchers likewise resolve
omission to 0.1. Thus periodic PCM/PLM full composition selects the complete
boundary-halo correction route when the coefficient is omitted. Explicit zero
retains the uncorrected path; negative library coefficients remain invalid.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. **2/2 CTests passed**: PPM halo 6.36 s and
PCM/PLM viscosity 3.36 s (9.76 s total), no failures or skips. Compiler warnings
remain. Provenance, Python syntax, contract and whitespace checks passed.

PCM/PLM periodic and no-slip/outflow fixtures check omitted versus explicit
0.1 exact state/temperature equality, zero producing a different result,
and omitted full-composition/boundary-hydro equality. Existing positive 0.037
mean-flux reconstruction, inventory, species closure and invalid rollback
checks pass. The existing PPM periodic/external-halo comparisons, pressure
tracing policies and viscosity/EB-mask analytic checks pass unchanged.

This covers these regular entry points only. Direct low-level PCM/PLM step
APIs, AMR/MPI/spray wrappers and EB library defaults still require review;
application inputs already supply their resolved coefficients. No full
regression matrix, external full-solver parity or physical validation was run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_pcm_plm_viscosity_3d test_reactive_ppm_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_(pcm_plm_viscosity|ppm_halo)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/regular-library-default-difmag-checks/build.log`, `LastTest.log`.
