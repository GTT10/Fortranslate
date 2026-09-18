# LES spatial filters — 2026-09-18

Implementation `582d71ca566075b78d856a30eaa8662cf6d3bd81`, based on `69dccc1`.
`les_filter_3d_mod` supplies standard box, truncated normalized Gaussian and
identity weights, plus component-first tensor-product filtering. With no origin
argument the field wraps periodically; `first_cell(3)` selects a crop of a
caller-supplied halo field. The output shape defines the crop. Failure preserves
output, including insufficient halos, invalid weights and nonfinite input.
Three separable passes replace the reference triple stencil sum; operation order
and floating-point roundoff differ. Repeated periodic samples are retained when
the filter is wider than the domain.

## Reference and boundaries

The pinned PelePhysics dependency is `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
Reference sources are `Source/Utility/Filter/Filter.cpp` (standard weight builders)
and `Filter.H` (`run_filter`). SHA256 values are in `references/les_filter.json`;
local downloaded sources are under `build/les-reference/`.
Box ratios are one or positive even integers; Gaussian ratios are positive even
integers. Ratios are bounded to 64. Box endpoint weights represent half cells;
ratio six supports the smoother used by pinned PeleC dynamic LES. Gaussian's
common prefactor cancels in normalization. Polynomial/optimized filter variants,
physical halo generation, MPI exchange, AMR filtering and dynamic LES field/flux
integration remain outside this increment.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged build configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `test_spray_core`;
selected H2O2 tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`.

Final CTest **2/2 passed**, 0.05 seconds, zero failures/skips: new `spray_les_filter`
and existing `spray_les_dynamic_algebra`. Filter checks cover analytic box weights,
impulse response, periodic total/constant preservation, Fourier attenuation,
width-six coefficients, exact supplied-periodic-halo agreement, Gaussian weight
ratios/normalization/Fourier response, identity, invalid ratios, insufficient halo,
nonfinite halo and oversized weight rollback. Initial 2/2 passed in 0.04 seconds;
the oversized-weight guard/test was added and the final build/tests rerun.
No build or test failures occurred. Provenance, Python syntax for the contract
checker, project contract and whitespace checks passed. No full configuration
matrix, MPI filter runtime, external solver parity or physical validation ran.

Logs: `build/les-filter-checks/{provenance,build,build-final,build-selected,build-selected-final,ctest,ctest-final}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(les_filter|les_dynamic_algebra)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
