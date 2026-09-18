# WALE and Vreman deviatoric stress kernels — 2026-09-18

Implementation `4d442cf`, base `5048cf0`.

## Reference and scope

Pinned [PeleC LES.H](https://github.com/Pele-Suite/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/LES.H)
provides WALE and Vreman closures. Source identity and SHA-256 are recorded in
`references/les_algebraic.json`; the downloaded reference is
`build/les-reference/LES.H` (SHA-256
`d9d24bc5dc24206408121fee64a9d616c78ac3edac99e16857bfa4aa93b69a17`).

`wale_stress` and `vreman_stress` accept G(i,j)=du_i/dx_j, density, spacing
and a coefficient, returning dynamic eddy viscosity and deviatoric stress.
Both use Delta=(dx*dy*dz)^(1/3), retaining the pinned isotropic filter convention
on anisotropic cells. WALE uses Cw squared; Vreman uses 2.5*Cs squared.
The returned stress is 2*mu*(S-tr(S)*I/3), corresponding to CI=0.
Grid fluxes, timestep estimates, model selection and application/restart controls
are not connected yet; the existing Smagorinsky path is unchanged.

The implementation normalizes G by its largest absolute entry before forming
homogeneous invariants. Vreman's Gram invariant is evaluated as squared minors,
which are nonnegative without cancellation. Zero gradient/coefficient returns
zero rather than dividing by zero. Input bounds match the existing Smagorinsky
kernel; invalid input returns zero outputs and false status. These numerical
choices are explicit deviations from direct evaluation of the pinned expressions.

## Focused checks

GNU Fortran 13.3.0, Debug/Ninja fixed tests-enabled and selected H2O2
tests-disabled configurations. Algebra tests run with invalid/zero/overflow traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(algebra|alternative_algebra)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Cases cover zero gradient, rank-one simple shear, diagonal strain on anisotropic
spacing, deviatoric stress and trace, WALE large-gradient homogeneity, and invalid
Vreman density. For diag(1,2,3), WALE uses S:S=14 and D:D=98/3;
Vreman uses B/(G:G)=49/14. These are analytic checks, not a C++ runtime oracle,
LES flow qualification or physical validation.

Initial build missed the test's IEEE finite-check import; corrected.
Initial tests were 1/2: the new trace check compared with zero through a relative
helper whose floor was 1e-25, rejecting roundoff of 1.7e-18. It now compares
sigma11+sigma22 with -sigma33 using the same 1e-13 relative tolerance.
Final tests **2/2 passed**, zero skips, 0.05 s total. Fixed/selected builds,
provenance, Python syntax, project contract and diff checks passed.

Logs: `build/les-alternative-checks/{provenance,build,build-final,build-test-final,build-selected,ctest,ctest-final}.log`.
