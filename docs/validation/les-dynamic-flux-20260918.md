# Periodic dynamic LES fluxes — 2026-09-18

Implementation `0b9874b12dd971faa9882a95f63cbe6d4a8b935a`, based on `a3129b8`.
`compute_dynamic_les_rhs_3d` connects smoothed dynamic coefficients to periodic
momentum and energy flux divergence. It arithmetically interpolates cell-centered
stress quantities, isotropic strain quantity, thermal flux and coefficients to
faces, then forms `sigma = Cs^2 alphaij - CI alpha I/3` and heat flux
`-cp (Cs^2/PrT) flux_T`. Face cp is evaluated at arithmetic face temperature and
composition. One face flux is shared by its neighboring cells. Mass and species
RHS are zero, following the pinned dynamic model rather than adding a species
closure. Failure preserves caller RHS and all inputs remain unchanged.

The reference is pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`LES.H::pc_dynamic_smagorinsky_sfs_term` and arithmetic interpolation in
`LES.cpp` dynamic LES stage 6; downloaded files are in `build/les-reference/`.
Harmonic interpolation, nonperiodic/MPI/AMR paths, time integration, timestep
selection and application dispatch remain open. Existing regularization/unit
caveats from the algebra record apply. This is a spatial RHS implementation,
not a completed dynamic-LES solver or physical validation.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled
`build/spray-amr-reactive-debug` built `test_spray_core`; selected H2O2
tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`. No configuration change.

Final focused CTest **2/2 passed**, 0.08 seconds, zero failures/skips:
`spray_les_dynamic_flux` and existing `spray_les_dynamic_field`.
The new test covers zero uniform RHS, active crossed-shear fluxes, exactly zero
mass/species RHS, periodic component conservation, negative kinetic-energy
production for the shear fixture, periodic translation, input preservation and
invalid-state rollback. It is not an independent external flux oracle or a
thermal-convergence qualification.

Initial CTest passed 1/2; the flux case failed its unchanged `1e-12` relative
conservation gate because energy RHS accumulation left a `-1.4900853e-11`
residual against absolute RHS sum `1.2810000000149`. A direct diagnostic rerun
reproduced it. Storing each face flux and then differencing left/right faces
before adding directional contributions removed this cancellation problem.
No tolerance was changed. Both builds and focused tests were rerun successfully.
Provenance, Python syntax for the contract checker, project contract and
whitespace checks passed. No full matrix or selected runtime test was run.

Logs: `build/les-dynamic-flux-checks/` contains `provenance.log`, `build.log`,
`build-diagnostic.log`, `build-final.log`, `build-selected.log`,
`build-selected-final.log`, `ctest.log` and `ctest-final.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_dynamic_(flux|field)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
