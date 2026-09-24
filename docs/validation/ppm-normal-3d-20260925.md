# Three-dimensional characteristic PPM normal fluxes — 2026-09-25

Implementation `cfbabe2`, base `44a1dc8`.
The new `compute_reactive_ppm_normal_fluxes_3d` API constructs
periodic time-centered normal face fluxes in x/y/z from NASA7 reactive states.
It rotates each seven-cell primitive stencil to the existing normal PPM
predictor, carries both tangential velocities, restores physical momentum
ordering and invokes the selected directional Riemann solver. Shock flattening
and bounded contact steepening retain the existing predictor semantics.
Each periodic boundary shares one flux value. Candidate fluxes publish only
after all directions succeed; invalid input/CFL leaves zero output fluxes.

This is the normal predictor needed by a CTU update. It is not yet a complete
3D PPM solver: transverse/corner corrections, conservative time stepping,
physical/AMR halos, MPI and application selection remain to be connected.
The current 3D PCM/PLM application behavior is unchanged.

## Provenance and scope

The reused normal characteristic reconstruction is described by
[decision 0012](../design_decisions/0012-reactive-characteristic-ppm-steepening-flattening.md).
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` sources were extracted
with `git show` into `build/ppm-normal-3d-checks/`:

- `PPM.H`: SHA256 `294cb223a95bf569b57700cbd9083ae77e562b7124f9cdda6aa338eada890e93`.
- `PPM.cpp`: SHA256 `312dd0ae66bd56b62ec0c3de8d7690ce8181ac8098124188d6576e8a79eab703`.

The source separates parabolic reconstruction/profile integration from tracing.
This increment reuses existing PeleF equations rather than claiming complete
reference tracing or multidimensional CTU parity.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed elementary H2/O2 mechanism.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_normal_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_normal_3d$' --output-on-failure --output-log build/ppm-normal-3d-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, test **0.71 s**, total **0.75 s**. An 8-cubed periodic N2
fixture checks analytic uniform Euler flux and a sinusoidal density wave,
with normal velocity 20 m/s and tangential velocities 5 and -7 m/s. Rotation
among all three axes agrees within 2e-13 relative scaling. Periodic face
identity is exact; the summed flux divergence closes within the same scaling.
Changing the trace interval changes the nonuniform prediction. An excessive
interval rejects and clears all face outputs. These are predictor/conservation
checks, not a time-evolved convergence study.

Build, provenance, Python syntax, project contract and whitespace passed;
no failures or skips. Compiler warnings remain. Artifacts:
`build/ppm-normal-3d-checks/{build.log,ctest.log,LastTest.log,PPM.H,PPM.cpp}`.
No selected build/execution, MPI run, full matrix, reference solver run or
physical validation was performed for this increment.
