# Periodic dynamic LES time integration — 2026-09-18

Implementation `7f657397c403a6a4c1a60fc9e090e7440dbad656`, based on `8101ba6`.
`advance_dynamic_les_ssprk2_3d` connects two-stage SSPRK2 with dynamic coefficients
recomputed from each stage. EOS recovery follows each candidate update. State
and temperature publish only after full success. Zero dt still validates the
model/state; invalid or oversized steps preserve both arrays.

The RHS optionally returns a frozen-coefficient diffusion estimate. At each
face it uses arithmetic mean `mu0 = rho Delta^2 |S|`, minimum adjacent density,
and the maximum of `4 (Cs^2 + CI)` and `(cp/cv) Cs^2/PrT`. Each stage must satisfy
`dt Dmax sum(1/dx^2) <= 0.2`. This is an implementation timestep heuristic that
includes isotropic stress, not a rigorous nonlinear stability bound for the
nonlocal dynamic coefficients. This increment is periodic library integration;
application/coupled-spray selection, MPI, physical boundaries and AMR remain open.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `test_spray_core`;
selected H2O2 tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`. Selected runtime was compiled, not exercised.

Final CTest **2/2 passed**, 0.17 seconds, zero failures/skips: new
`spray_les_dynamic_step` and existing `spray_les_dynamic_flux`. The step test uses
variable density, pressure and three velocity components; checks component
conservation, exact density/species preservation, zero step, excessive step and
invalid zero-step configuration rollback. One/two/four steps over the same
interval satisfy the second-order refinement ratio gate of 3 to 5. Existing
spatial flux conservation and crossed-shear dissipation checks remain.
Provenance, Python syntax for the contract checker, project contract and
whitespace checks passed. No full matrix, external solver parity or physical
validation was run.

## Initial failures and diagnosis

The first 1/2 run failed an inappropriate monotonic kinetic-energy assertion:
for the divergence-free initial shear fixture, dt was `4.27520125e-7`, kinetic
energy changed from `3000` to `3000.0080814524`, while initial instantaneous
production was `-1.578`. Isotropic SGS stress can accelerate flow and exchange
kinetic/internal energy; monotonic kinetic energy is not a general invariant.
The test was replaced by total/component conservation and temporal refinement.

The initial symmetric shear refinement then failed the 3-to-5 ratio gate:
coarse/fine differences were `6.9005182e-8` and `9.7788870e-9` (ratio about 7.06).
That fixture does not exhibit the expected generic leading temporal error at
this interval. The step fixture was changed to nonsymmetric variable density,
pressure and three-component velocity; the unchanged refinement gate passed.
The original spatial shear test remains separate. No numerical implementation
or conservation tolerance was changed to address these test failures.

## Commands and artifacts

Logs under `build/les-dynamic-step-checks/`: `provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, diagnostic build/run logs,
`build-final.log`, `ctest-final.log` (the failed symmetric refinement),
`build-refined-fixture.log`, and final `ctest-refined-fixture.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_dynamic_(step|flux)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
