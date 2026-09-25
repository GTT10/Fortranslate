# Aperture-aware EB PLM normal prediction — 2026-09-25

Code `fe19ee1`, parent `2c18ce0`. New
`predict_reactive_eb_plm_cell_3d` ports the normal `pc_plm_d_eb` prediction
from PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/PLM.H:314-485`, inspected with `git show` in
`/mnt/sn850/Projects/gtt2/external/PeleC`.

Inputs are a five-cell primitive stencil, a separate internal-energy-density
stencil, covered-cell flags, sound speed, direction, spacing, interval,
flattening and two face apertures. The predictor rotates velocities into the
normal direction, uses the pinned covered-cell slope order policy, and returns
primitive states plus independently traced internal-energy density. The latter
retains the enthalpy/acoustic/contact terms needed for EB transverse work.
Only open faces receive states; closed outputs are zero. A fully closed cell
can have poisoned covered data without being read. Invalid controls/CFL leave
all outputs zero.

The shared regular PLM predictor accepts optional open-face masks, defaulting
to both open. It checks/normalizes composition only on active faces, preserving
the existing normalized-composition EOS boundary. The EB predictor accepts
orders 1/2/4 and does not multiply states by aperture; aperture weighting belongs
to conservative flux integration.

This is a normal-prediction API. Existing EB hydro face loops still use
cell-centered states. Wide halo preparation, Riemann/CTU integration,
transverse corrections, centroid flux interpolation and application dispatch
remain. This is not a claim of complete EB Godunov or compiled external parity.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. **2/2 CTests passed**, no failures/skips:
EB hydro regression 0.92 s and PPM halo unit 6.36 s (7.32 s total).
The EB regression adds all three directions and four open/closed face masks,
linear spatial density/velocity/internal-energy predictions, second-order
fallback with poisoned covered values, finite-time fully flattened states,
analytic internal-energy contact-wave integration, CFL rejection and fully
covered data avoidance. Existing hydro conservation/rollback and regular
PPM/PLM checks remain. No tolerance or golden changed.

Provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. No application selection, MPI predictor dispatch,
full configuration matrix, compiled external oracle or physical validation was
run for this API.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d test_reactive_ppm_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(regression_reactive_eb_hydro_3d|unit_reactive_ppm_halo_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-plm-normal-checks/{build.log,rebuild.log,LastTest.log}`.
