# Static LES face flux extraction — 2026-09-18

Implementation `67c3494f60be4e76a7a03434eddbca879976573d`, based on `3f2f69c`.

The ghosted static LES RHS API returns optional x/y/z face flux arrays, all three
or none. Exact shapes are required, with one extra face in the normal direction.
The shared face loop writes private candidate fluxes and preserves its existing
RHS arithmetic. The wrapper publishes fluxes only after successful evaluation;
invalid inputs preserve caller face buffers. Disabled LES returns zero fluxes.
No output allocation occurs when callers omit the face arrays. Fluxes exclude
area/time factors and retain the positive coordinate orientation. Boundary masks
select ledger contributions only. This is an AMR synchronization building block,
not an implemented LES-AMR integrator or integrated SSPRK2 flux interface.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, existing configurations unchanged.
Fixed tests-enabled `test_spray_core` and selected H2O2 tests-disabled
`pelef_selected_spray_runtime` built. Selected runtime was not run.
Focused CTest **3/3 passed**, 0.40 seconds, zero failures/skips:
`spray_les_boundary`, `spray_les_face_fluxes`, `spray_les_isotropic`.

The new test covers Smagorinsky/WALE/Vreman with CI 0.1, periodic and outflow
halos. Requesting faces preserves RHS and diffusivity exactly. Their divergence
recovers every cell RHS within `1e-12` relative with unit floor; area-weighted
boundary faces reproduce the ledger to the same tolerance. Invalid-model and
partial-output rejection preserve sentinels; disabled LES returns zero faces.
Existing boundary/isotropic checks also pass. No tolerances changed.
Provenance, Python syntax, project contract and whitespace checks passed.
No builds/tests failed. No AMR integration, MPI execution, full matrix, external
solver parity or physical validation ran.

## Commands and artifacts

Logs: `build/les-face-flux-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(face_fluxes|boundary|isotropic)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
