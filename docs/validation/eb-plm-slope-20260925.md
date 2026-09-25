# Covered-cell PLM slope policy — 2026-09-25

Code `1b30c6a`, parent `144ca41`. `compute_pelec_plm_slopes` ports
`plm_slope` / `plm_slope_eb` from pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, `Source/PLM.H:21-115`.
The local reference is `/mnt/sn850/Projects/gtt2/external/PeleC`; source was
read with `git show` at that revision.

The order-four stencil needs five uncovered cells. If an outer cell is covered,
a three-cell order-two stencil is used; a covered center or immediate neighbor
forces zero slope. Order one always returns zero. The order-two reference
formula retains 4/3 times the centered difference before limiting. The new
routine avoids reading covered/inactive stencil entries and validates orders
1/2/4, shapes and flattening. Invalid inputs return zero slopes.

The existing `reconstruct_pelec_plm_cell` now calls this shared implementation
and accepts optional order/covered controls. Its omitted policy remains order
four. Normal wave tracing and composition normalization remain unchanged.
This adds the slope prerequisite for EB PLM; the EB hydro face loops still use
cell-centered Riemann states. Aperture-aware EB normal prediction, transverse
corrections, centroid flux interpolation and application dispatch are not
implemented by this increment. No complete EB Godunov parity is claimed.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. **1/1 CTest passed**, 6.20 s (6.24 s total),
no failures/skips. The existing PPM halo test adds all 32 covered-cell masks
for three orders (96 slope cases), with covered values filled with NaNs.
Linear-field analytic slopes distinguish order four, the reference order-two
factor and zero fallback. Normal reconstructed edges exercise order-two/zero
fallback, and invalid order rejects. Existing PPM halo, boundary, flattening,
PLM predictor and composition checks also pass. No tolerances/goldens changed.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed. Compiler warnings remain. No compiled external slope oracle, full
matrix, external full-solver comparison or physical validation was run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-plm-slope-checks/{build.log,LastTest.log}`.
