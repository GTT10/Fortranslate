# CTU transverse ordering — 2026-09-25

Source: base `66f8cdc`, implementation `614d5e6`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-transverse-order-checks/`.

Source inspection of pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` showed that the final
`pc_transdd` operations use transverse axes Y,Z for X faces; X,Z for Y faces;
and X,Y for Z faces. PeleF's cyclic ordering had Z,X for Y faces. Although
flux differences are summed, the pressure-work terms are subtracted sequentially,
so exchanging their order can change roundoff. The pressure-traced path now
uses the reference order. The conservative/EOS audit path retains its ordering.

Inspected via `git show` in `/mnt/sn850/Projects/gtt2/external/PeleC`:

- `Source/Godunov.H`: SHA256
  `1bf575907a03764f17d7fe1580c5f546755846d6c9ee4badb42d8499fb4a9f27`.
- `Source/Godunov.cpp`: SHA256
  `0d6a5fa5ec900c0af2f97d51c1493da074c876d430d3955cfcd4148e429f7efe`.

Both inspected files are retained in the artifact directory. This is source
alignment, not a measured full-CTU external error reduction.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d test_reactive_ppm_halo_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_step_3d|unit_reactive_ppm_halo_3d)$' --output-on-failure`:
  **2/2 passed**, 3.06 s and 5.03 s respectively, 8.14 s total; no failures/skips.
  Existing tests cover multistep conservation, species closure, periodic/physical
  transport composition, exact halo/subdomain flux parity and rollback.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

## Newly located remaining boundary work

After final CTU solves, `pc_umeth_3D` replaces hard inflow/user-boundary fluxes
using `pc_low_order_boundary`: PLM slopes, `pc_plm_d`, then `pc_cmpflx`, without
CTU transverse corrections for those replacement states. PeleF currently uses
its full PPM CTU result at inflow faces. This reference boundary policy is not
yet connected and prevents claiming complete physical-boundary reference parity.
The reference comment says no half-time predictor, but its actual call passes
`dt` into `pc_plm_d`; implementation must follow the active arithmetic rather
than assume the comment means `dt=0`.

No MPI/application rerun, selected build, full matrix, compiled full-CTU oracle,
convergence or physical validation was performed. Compiler warnings remain.
