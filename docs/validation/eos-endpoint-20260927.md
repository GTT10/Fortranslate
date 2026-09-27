# Primitive conversion at exact EOS endpoints, 2026-09-27

Implementation: `2fb2931efdf233a3c08b7f8a7abc1368bd80a898`, based on
`8a138647`, local `main`.

## Defect and repair

The HIT input work exposed a failure before velocity loading: pure N2 at
300 K and 101325 Pa was rejected by the uniform initial-state constructor.
The focused reproducer established that forward EOS density `p/(R*T)` inverted
through `p/(rho*R)` to `299.99999999999994 K`. The common full-H2/O2 interval
is [300,3500] K, so strict property evaluation rejected the result. The upper
endpoint similarly inverted to `3500.0000000000005 K` at 300000 Pa.

`reactive_primitive_to_conserved` now attempts endpoint recovery only after
ordinary property evaluation fails. It first requires the computed temperature
to lie outside the common species interval. It then evaluates the forward EOS
at the corresponding endpoint and requires **exact equality** with the input
density. Only that endpoint is passed to the unchanged property evaluator.
Interior successful evaluations keep their previous path. Density, momentum,
composition, NASA data, temperature limits and numerical tolerances are not
changed. A distinct out-of-range density remains rejected; no broad clipping
or tolerance-based acceptance was added.

## Completed checks

Existing `build/amr-cvode-app-debug`, GNU Fortran 13.3.0 Debug, MPI/spray/
SUNDIALS enabled, Open MPI 4.1.6, SUNDIALS 7.2.0, selected full H2/O2.
Compiler/backend configuration was unchanged. Existing real-comparison and
executable-stack warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug \
  --target test_reactive_general_eos pelef_reactive_3d pelef_reactive_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/amr-cvode-app-debug \
  -R '^unit_reactive_general_eos$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The added endpoint regression failed on the unmodified solver at 300 K.
After repair, **1/1 CTest passed**, 0 unexpected failures and 0 skips.
It tests both common endpoints at 101325, 300000 and 1000000 Pa, exact density
and species-density preservation, and rejection after perturbing density by
32 representable spacings toward the out-of-range side. Existing transformation
round trips, characteristic operations, overflow/nonfinite rejection and
zero-output failure contracts are part of the same test and passed.

Two additional direct application runs use periodic 8×4×4 cells, pure N2 at
300 K and 101325 Pa, mean velocity [3,2,1] m/s, inert SDC with molecular
transport, and final time 1e-10 s. Fixed/selected outputs are byte-identical.
Inputs and logs are `build/eos-endpoint-checks/{ascii,selected}/`; binaries
are `build/amr-cvode-app-debug/pelef_reactive_3d{,_selected}`. Each was invoked
with its directory's `input.nml`, working in that directory with
`LD_LIBRARY_PATH`, `LIBRARY_PATH` and `DISPLAY` unset. These runs isolate the
uniform initializer and do not include HIT velocity perturbations, which can
physically cool some cells below the mechanism's lower bound.

Provenance, builds, syntax, project contract and whitespace checks passed.
The pinned YAML SHA-256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eos-endpoint-checks/{red-build.log,red-ctest.log,build.log,ctest.log,LastTest.log}`.
This is a bounded floating-point conversion repair, not a new thermodynamic
validity range or a full flow/reference qualification.
