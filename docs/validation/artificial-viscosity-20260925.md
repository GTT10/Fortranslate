# Regular-cell artificial-viscosity kernel — 2026-09-25

Source: base `5eff6af`, implementation `2228208`. GNU Fortran 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/artificial-viscosity-checks/`.

`adjust_reactive_artificial_viscosity_3d` accepts conserved states with one ghost
cell per side, spacing, nonnegative coefficient and six physical no-slip masks.
It averages velocity derivatives at nodes, averages four node divergences at
each face, and adds `dx*difmag*min(divu,0)*(U_right-U_left)` to conserved fluxes.
Tangential momentum diffusion is suppressed at marked no-slip faces. Species
fluxes are rescaled by mass flux divided by their sum when that sum is nonzero.
Fluxes stay intensive; callers own area factors. A zero coefficient is an exact
identity after input validation. Invalid inputs or nonfinite candidates preserve
all three supplied flux arrays transactionally.

This is a callable kernel, not yet dispatched by applications. The reference
uses `difmag=0.1` by default; coefficient input and local physical-wall masks
must be connected before enabling it across production paths. No numerical
results from existing applications changed in this increment.

## Reference

Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, read with
`git show` from `/mnt/sn850/Projects/gtt2/external/PeleC`:

- `Hydro.H::pc_divu`, `pc_norm_spec_flx`: SHA256
  `a3f60b710e7bbf7048222f681b110188ff1f8f9b488e49175cb920df7c5eb6b1`.
- `Godunov.H::pc_artif_visc`: SHA256
  `1bf575907a03764f17d7fe1580c5f546755846d6c9ee4badb42d8499fb4a9f27`.
- `Hydro.cpp::pc_adjust_fluxes` applies viscosity after CTU, then species
  normalization and face-area scaling. `Params/_cpp_parameters` sets `difmag=.1`.

The first two headers are retained with artifacts. PeleF has no transported
UTEMP/auxiliary slots in this conserved vector. The reference explicitly skips
those slots; all represented PeleF conserved components receive the correction.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure`:
  **1/1 passed**, 5.64 s test / 5.69 s total; no failures/skips.
- Added analytic linear velocity with divergence -6 and directional density
  gradients: exact expected mass diffusion, conserved-component jumps, species
  closure, no-slip tangential exclusion, zero correction in expansion and
  invalid-density rollback. Existing halo/PLM/inflow checks remain.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No MPI, selected build, application execution, compiled viscosity reference,
full matrix, convergence or physical validation was run. Compiler warnings remain.
