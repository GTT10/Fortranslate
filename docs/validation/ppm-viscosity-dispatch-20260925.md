# PPM artificial-viscosity dispatch — 2026-09-25

Source: base `26c5ed1`, implementation `146d1e9`. GNU Fortran 13.3.0,
OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/ppm-viscosity-dispatch-checks/`.

Periodic and supplied-halo CTU APIs now accept an optional `difmag`. The default
is 0.1 with pressure tracing (`pelec` by default), zero for the conservative/EOS
policy. An explicit zero disables viscosity; negative/nonfinite coefficients
reject with unpublished fluxes. Positive explicit coefficients also work with
the conservative/EOS policy through the shared halo implementation.

Viscosity runs after final CTU and inflow PLM replacement, before publishing
fluxes. Regular physical, serial/MPI two-level AMR, serial/MPI patch-tree common
and recursive, and regular MPI spray callers now forward physical no-slip masks.
Both no-slip and isothermal walls suppress tangential artificial momentum flux;
internal patch/rank faces do not. Tree masks reuse integer domain geometry;
sparse/slab masks use physical extents and rank ownership. Periodic calls use
no physical walls. Flux correction is still intensive and callers retain area,
reflux, boundary-budget and transactional state handling.

Applications selecting `characteristic_ppm` with `pelec` therefore now include
reference-default artificial viscosity. This changes results relative to older
executables; same-input cross-version bitwise continuation is not established.
The coefficient is configurable at CTU API level only. Namelist input and
restart binding of a user coefficient remain unfinished.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d test_mpi_amr_reactive_3d test_mpi_amr_patch_tree_hydro_3d test_mpi_amr_patch_tree_subcycle_3d test_spray_mpi_ppm --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_halo_3d|mpi_amr_ppm_hierarchy_3d_np8|unit_mpi_amr_patch_tree_ppm(_subcycle)?_3d_2|spray_mpi_ppm_4)$' --output-on-failure`:
  **5/5 passed**, 48.12 s total, no failures/skips.
  - Halo unit, 5.62 s: default flux changes for compression; zero-coefficient
    CTU plus explicit viscosity matches the default composition exactly;
    negative coefficient leaves outputs zero. Existing analytic/wall tests remain.
  - Eight-rank two-level hierarchy, 19.72 s: serial parity including inflow,
    reacting transport, boundary/reflux budgets, empty fine ranks and rollback.
  - Two-rank tree common/recursive, 4.21 / 13.23 s: exact serial fields,
    boundary/reflux diagnostics and recursive step counts.
  - Four-rank regular spray hydro, 5.30 s: exact serial fields for both solvers,
    six Mach-inflow orientations and collective rollback on short slabs.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No selected build/execution, application restart, full matrix, compiled complete
CTU reference, convergence or physical validation was run. Compiler warnings
remain. PCM/PLM and EB artificial-viscosity dispatch are not added by this change.
