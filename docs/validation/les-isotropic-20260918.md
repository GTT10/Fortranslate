# Static LES isotropic stress — 2026-09-18

Implementation `546a826ad8e645c137c9cb7026d69fb90e8a66e9`, based on `da3484e`.

Smagorinsky, WALE and Vreman stress APIs now accept trailing optional `ci`.
The shared serial periodic/physical-boundary RHS, supplied-halo RHS and SSPRK2
APIs forward it. CI must be finite and in `[0,1]`; its default is zero.
The reference term is `-(2/3)*CI*rho*delta^2*|S|^2 I`, where
`delta=(dx*dy*dz)^(1/3)` and `|S|^2=2*sum(S_ij^2)`.
It changes normal stress and pressure work but leaves eddy viscosity and the
existing turbulent species/heat diffusion unchanged. Positive CI works when the
eddy coefficient is zero. Zero CI preserves the previous disabled-step identity.

The pinned reference is PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/LES.H`; source SHA and routines are in `references/les_isotropic.json`.
Local cached source SHA was checked in this increment. The existing species
flux is retained, so this is not full C++ flux parity.

Each stage adds `4*CI*rho_face*delta^2*|S|/min(rho_L,rho_R)` to the previous
maximum-diffusivity estimate. This is a frozen-gradient interval heuristic,
not a proof of stability for the nonlinear model. EOS failure or an excessive
interval leaves state and temperature unchanged. MPI wrappers and coupled/
application CI selection remain to be connected.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `test_spray_core` and selected H2O2 tests-disabled
`pelef_selected_spray_runtime` built successfully; selected runtime was not run.

Focused CTest **5/5 passed**, 0.25 seconds, zero failures/skips:
`spray_les_isotropic`, `spray_les_algebra`, `spray_les_alternative_algebra`,
`spray_les_zero`, `spray_les_boundary`.
The new case checks all three models: anisotropic-spacing analytical diagonal
stress, unchanged viscosity/shear, exact zero-CI RHS/diffusivity equivalence,
nonzero RHS and SSPRK2 update with zero eddy coefficient, unchanged mass flux,
periodic conservation, outflow boundary ledger, excessive-dt rollback and NaN
CI rejection even at zero dt. Scalar checks cover rigid rotation and invalid
CI bounds. No test tolerances were relaxed. A preliminary 5/5 run passed before
adding the outflow ledger check; the final 5/5 above includes it.
Provenance, Python syntax, contract and whitespace checks passed. No build/test
failed. No MPI execution, full matrix, external parity or physical validation ran.

## Commands and artifacts

Logs: `build/les-isotropic-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(isotropic|algebra|alternative_algebra|zero|boundary)$' --output-on-failure
sha256sum build/les-reference/LES.H
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
