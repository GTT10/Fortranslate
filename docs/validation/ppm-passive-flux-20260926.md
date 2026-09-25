# Godunov passive upwinding and final normalization — 2026-09-26

Implementation `222ee031`, parent `98df1b7b`, local main, version 0.246.0.
Pinned PeleC `pc_cmpflx` discards the species flux returned by `riemann` and
calls `pc_cmpflx_passive` using contact speed and original face species weights.
PeleF previously used the inner Riemann species flux directly. The primitive
Riemann API now has optional `passive_upwind`; regular CTU, its inflow replacement
and EB PLM normal/transverse face solves select it. The raw kernel API default
remains available for direct Riemann comparison. The contact speed is separate
from the sampled Godunov velocity used for pressure work.

Pinned Hydro.cpp normalizes final species fluxes regardless of `difmag`.
Regular pressure-traced CTU now does this even at zero viscosity. The viscosity
library's zero-coefficient identity contract is retained. EB CTU already had
its own zero-coefficient normalization. Intermediate face species are not
normalized before transverse flux differences. Cell acceptance and transactional
rollback remain unchanged.

## Focused evidence

GNU Fortran/C++ 13.3.0, Python 3.12.3, existing Ninja Debug tests-enabled
`build/spray-amr-reactive-debug`. The expanded Riemann tool compiles the unchanged
pinned `riemann` and `pc_cmpflx_passive` bodies, using the existing NASA7 SI
adapter. **74 cases / 1,036 values passed**, maximum scaled error
**3.111830054382345e-14**, fixed tolerance `3e-12`. Cases include signed and
nonunit positive weights in both velocity directions and stationary interfaces.
Using the raw Riemann output for these same cases has maximum scaled difference
**0.2** from the passive-upwind reference, confirming the dispatch distinction.
This comparison excludes physical-boundary state replacement inside pc_cmpflx.

Four distinct relevant tests passed:

- EB hydro regression: 22.92 s.
- PPM multistep update: 3.31 s, now alternating zero/default viscosity.
- Raw Riemann unit: 0.01 s.
- Final halo unit: 6.58 s test / 6.63 s total.

The first combined run passed 3/4. Its halo assertion assumed that zero-viscosity
flux was unnormalized and could be normalized only after adding viscosity.
Since those operations do not commute for species, that obsolete assertion was
replaced with exact hydro-component composition and separate three-axis species
closure checks for zero viscosity (3e-14 scaled bound). The final halo run passed.
No solver tolerance was relaxed and no tests were skipped. Initial failure log
is retained; the combined run is not reported as passing.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_pelec_riemann test_reactive_ppm_step_3d test_reactive_ppm_halo_3d test_reactive_eb_hydro_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_riemann_kernel.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_pelec_riemann --work-dir build/ppm-passive-flux-checks/oracle --face-flux
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_pelec_riemann|unit_reactive_ppm_step_3d|unit_reactive_ppm_halo_3d|regression_reactive_eb_hydro_3d)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_pinned_riemann_kernel.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace passed. Build warnings remain.
Artifacts: `build/ppm-passive-flux-checks/` with build logs, `comparison.log`,
initial `ctest.log`, successful `ctest-halo.log`, `unit-halo.log`, pinned
`Hydro.cpp`, and `oracle/` with both headers, driver, inputs, raw/reference/final
outputs, dispatch difference and hash-bearing report.
Reference commit: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
No selected/MPI/application run, installed audit, full-CTU assembly comparison
or physical validation ran. Complete PeleC parity remains unproven.
