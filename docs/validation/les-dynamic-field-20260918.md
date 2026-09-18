# Dynamic LES coefficient fields — 2026-09-18

Implementation `518b9bbf83462dbd9f5f2745dcc2f053980d2cd2`, based on `5a04809`.
`compute_dynamic_les_coefficients_3d` evaluates periodic regular-grid coefficients
in the ordering `Cs^2`, `CI`, `Cs^2/PrT`. It recovers primitive variables from
original and test-filtered conserved state, forms centered velocity/temperature
derivatives, evaluates grid and test quantities, filters grid quantities,
contracts Germano tensors and smooths coefficients with the pinned ratio-six
box kernel. Filtering conserved state before EOS recovery follows pinned PeleC
`LES.cpp` dynamic LES stages 2–5, rather than filtering primitive temperature.
Default grid filter ratio is one; callers select test filter type and ratio.
The API preserves caller coefficients on failure and does not modify inputs.

References remain pinned to PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`
and its PelePhysics dependency; local source is `build/les-reference/LES.cpp`.
The algebra and filter manifests are `references/les_dynamic.json` and
`references/les_filter.json`. Existing pinned regularization/unit caveats apply.
This increment does not connect dynamic coefficients to face fluxes, time
integration or applications. Physical boundaries, MPI and AMR coefficient fields
remain open. No full dynamic-LES solver or physical validation is claimed.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled build
`build/spray-amr-reactive-debug` built `test_spray_core`; selected H2O2
tests-disabled build `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`. Configurations were unchanged.

Final focused CTest **2/2 passed**, 0.06 seconds, zero failures/skips:
`spray_les_dynamic_field` and existing `spray_les_dynamic_algebra`.
The field case checks the uniform-state floor, nontrivial coefficients on a
variable-density/temperature/velocity field, positivity and CI upper bound,
input preservation, exact periodic translation, box/Gaussian distinction and
transactional rejection of invalid state/filter settings. The initial run also
passed 2/2 in 0.06 seconds; primitive-range guards were then added before gradient
formation and both builds/tests rerun. No failed builds or tests occurred.
Mechanism provenance, Python syntax for the contract checker, project contract
and diff whitespace checks passed. Selected runtime was compiled, not exercised;
no full matrix, external solver parity or physical validation ran.

Logs: `build/les-dynamic-field-checks/{provenance,build,build-final,build-selected,build-selected-final,ctest,ctest-final}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_dynamic_(field|algebra)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
