# Physical PPM Riemann faces — 2026-09-26

Implementation `0db12d8c`, parent `ca14f2b9`, local main, version 0.246.0.
Pinned PeleC Godunov.H `pc_cmpflx` corrects hydro face variables immediately
before the Riemann solve at physical walls/outflow. Ghost construction alone did
not connect that rule to PeleF's regular pressure-traced CTU face solves.

The prepared-halo API now accepts reflecting/outflow masks, corrects hydro
variables at every normal/corner/final face solve and preserves the original
species face weights. Wall exterior normal velocity is the negative of the
interior value; other hydro variables copy the interior. Ordinary outflow copies
all five interior hydro variables. Explicit pressure outlets keep their existing
EOS-adjusted exterior state and are excluded from the ordinary outflow override.
Overlapping wall/outflow masks reject with cleared outputs.

Regular boundary stepping, serial two-level fine stepping, serial/MPI patch-tree
common/recursive stepping, MPI sparse coarse/fine stepping and regular MPI spray
hydro forward the physical masks. Tree integer layouts and slab/rank endpoints
restrict them to actual domain faces, not internal patch boundaries. Existing
inflow PLM replacement and non-PeleC reconstruction behavior remain unchanged.

## Focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, existing Ninja Debug tests-enabled
`build/spray-amr-reactive-debug`. Required targets rebuilt, including the regular
MPI spray unit (built, not executed). **5/5 CTests passed**, 46.27 s total:

- PPM step: 3.24 s.
- PPM halo: 6.72 s.
- Sparse two-level PPM hierarchy, 8 ranks: 19.21 s.
- MPI tree PPM common, 2 ranks: 4.14 s.
- MPI tree PPM recursive, 2 ranks: 12.91 s.

New halo checks supply deliberately unreflected donors with nonzero velocity.
All six wall faces have exactly zero mass, total-energy and species fluxes at
zero viscosity. A mismatched exterior velocity at ordinary outflow produces the
analytic interior mass flux; the explicit pressure-outlet variant retains the
exterior donor. Conflicting masks reject. Existing update conservation,
serial/MPI parity, physical-boundary/pressure-outlet cases, idle owners and
rollback assertions execute unchanged. No failed or skipped checks.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d test_reactive_ppm_step_3d test_mpi_amr_patch_tree_hydro_3d test_mpi_amr_patch_tree_subcycle_3d test_mpi_amr_reactive_3d test_spray_mpi --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_halo_3d|unit_reactive_ppm_step_3d|unit_mpi_amr_patch_tree_ppm_3d_2|unit_mpi_amr_patch_tree_ppm_subcycle_3d_2|mpi_amr_ppm_hierarchy_3d_np8)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. No tolerances or
golden data changed. Build warnings remain. Artifacts under
`build/ppm-physical-riemann-checks/`: `build.log`, `build-final.log`, final
`build-pressure.log`, `ctest.log`, `unit.log` and pinned `Godunov.H`.
Reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Godunov.H SHA256: `1bf575907a03764f17d7fe1580c5f546755846d6c9ee4badb42d8499fb4a9f27`.

No selected build, spray runtime/application rerun, installed audit, full-CTU
external comparison or physical validation ran. Complete PeleC parity remains
unproven.
