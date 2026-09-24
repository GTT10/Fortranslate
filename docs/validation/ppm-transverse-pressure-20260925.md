# Separate transverse pressure work — 2026-09-25

Source: base `ba909a7`, implementation `5895fe3`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-transverse-pressure-checks/`.

The periodic and four-halo CTU APIs accept optional `pressure_tracing=.true.`
with the `pelec` solver. Normal and corner Riemann solves retain interface
pressure/normal velocity. Single and double transverse corrections apply the
separate pressure-work equation, while density, momentum, species and total
energy receive conservative flux differences. The primitive Riemann entry then
consumes the corrected pressure directly. The correction helper also exposes
corrected internal energy independently. Invalid states leave outputs zero.
The default conservative/EOS correction and application dispatch remain unchanged.

Pinned source: PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/Godunov.H::pc_transdo`, `pc_transdd`, `pc_cmpflx`, read with `git show`
from `/mnt/sn850/Projects/gtt2/external/PeleC`. Header SHA256:
`1bf575907a03764f17d7fe1580c5f546755846d6c9ee4badb42d8499fb4a9f27`.
Single pressure floor is `epsilon*1e-100`; double floor is `tiny`.
This covers regular cells without forcing sources. Normal conserved predictors
still come from reconstructed primitive states through the EOS: separate normal
internal-energy tracing and full reference parity are not established.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure`:
  **1/1 passed**, 4.72 s test / 4.76 s total; no failures or skips.
- Analytic single/double pressure work, independent internal energy, conservative
  species correction and invalid-density rollback; nonuniform 3D periodic/halo
  and central-subdomain exact flux parity for both old and new policies;
  new-policy solver mismatch/CFL rollback. Existing six-face Mach-inlet inventory
  and rollback, unused-corner and stationary-wall checks remain in the test.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No selected-runtime execution, MPI, full matrix, compiled external oracle,
convergence or physical validation was run. Compiler warnings remain.
