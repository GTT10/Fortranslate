# EB CTU physical boundary corrections — 2026-09-25

## Implementation

Code `feb883c`, parent `144d8b0`, introduces `eb_plm_boundary_3d` with
physical cell bounds in working-geometry coordinates and six face types.
Normal, single-transverse and final CTU Riemann solves now apply the pinned
`pc_cmpflx` physical-boundary correction to local primitive arguments:
wall/symmetry faces mirror the normal velocity and copy interior density,
pressure and tangential velocity; outflow copies all five components.
Reconstructed states, independently traced internal energy and species
fractions remain unchanged. Omitted boundary data retain the previous path.
Interior/periodic/inflow labels require caller-prepared donors and apply no
additional Riemann correction. Invalid labels or reversed bounds fail without
publishing output.

Reference: PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/Godunov.H:303-445`, especially the boundary corrections at 384-419.
The local `isothermal_wall` label uses the reference no-slip hydro treatment;
thermal wall-state preparation remains the caller's responsibility.

`advance_reactive_eb_plm_ctu_3d` forwards the boundary descriptor. Its viscosity
stage converts physical bounds to the cropped margin's face indices:
`lower-3`, `upper-2`. The shared artificial-viscosity API accepts optional
physical face positions, so tangential momentum correction is excluded at
actual no-slip/isothermal walls inside the extended region. Existing callers
default to faces zero and N as before. The optional pair is validated together.

## Checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6 / Python 3.12.3, existing Ninja Debug build
`build/spray-amr-reactive-debug`; no compiler or backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_pcm_plm_viscosity_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_pcm_plm_viscosity_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2 distinct tests passed**. EB regression: 10.15 s / 10.19 s total.
Regular PCM/PLM viscosity: 7.80 s / 7.84 s total. Before complete-update and
interior-wall-mask integration, the EB flux-only increment also passed
(9.20 s / 9.24 s); this is not counted as an additional test.

The new regression exercises five boundary labels on six physical faces,
normal and final CTU solves, and unaffected interior faces. Independently
corrected primitive arguments feed the existing Riemann solver for expected
flux, pressure and velocity; spatially varying species expose accidental
composition copying. It also checks invalid-label rollback, five complete
prepared updates and their inventory balance, and tangential viscosity
exclusion on six interior physical faces while ordinary faces and working
region edges retain that correction. The regular test verifies existing
PCM/PLM default/explicit viscosity and rollback behavior.

Provenance matched SHA
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Python syntax, project contract and diff checks passed. No build/test failures,
no skips and no tolerance changes in this increment. Compiler warnings remain.
Artifacts are under `build/eb-physical-flux-checks/`: `build.log`,
`integration-build.log`, `regular-build.log`, `LastTest-flux.log`,
`LastTest-eb.log` and `LastTest-regular.log`.

## Boundary of the claim

This is prepared-input library evidence. The boundary fixture checks the
correction algebra and inventory, not a physical wall-flow solution. The
application still needs extended geometry/state halo preparation and input
selection, followed by coupled MPI/AMR dispatch. No compiled external oracle,
application matrix, installed-binary audit or physical validation was run.
