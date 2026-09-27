# HIT Reynolds/Prandtl transport construction, 2026-09-27

Implementation `02c2946680e8c11789bcb1a7d5857776a4d91ab9`, based on
`5958c904`, local `main`.

## Behavior and reference

`reactive_constant_transport` now accepts `mode='hit'`, `lambda0`,
`reynolds_lambda0` and `prandtl`. The HIT initializer returns the uniform
initial density, target component RMS velocity and mixture cp. After successful
initialization, the application computes:

- `urms = mach_t0 * sound_speed / sqrt(3)`
- `viscosity = rho * urms * lambda0 / reynolds_lambda0`
- `conductivity = viscosity * cp / prandtl`

These are the equations in `Exec/RegTests/HIT/prob.cpp` at pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, inspected with `git show`.
PeleF quantities are SI. The velocity file normalization `uin_norm` affects
sample scaling, not the target RMS used in the Reynolds definition.
`mode='specified'` remains the default for explicit coefficients. Mixed explicit
and derived coefficients, nonpositive parameters, missing/zero initialized
HIT scales and overflowing intermediate products/quotients reject before output.

The initializer publishes the optional scales only with a successfully accepted
velocity field. No alternate sound-speed or Mach setting was added to transport.
Species diffusion flags remain caller-controlled; disable them for the pinned
HIT reference's zero species diffusion. Gamma-law EOS, forced-turbulence
statistics and external evolved-field parity remain outside this evidence.

## Completed checks

Existing `build/amr-cvode-app-debug`: GNU Fortran 13.3.0 Debug, Open MPI 4.1.6,
SUNDIALS 7.2.0, MPI/spray/SUNDIALS enabled, selected full H2/O2. No compiler or
backend configuration changed. Existing real-comparison and executable-stack
warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug \
  --target pelef_reactive_3d pelef_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest \
  --test-dir build/amr-cvode-app-debug -R '^application_hit_transport_scales$' --output-on-failure
python3 -m py_compile tests/application/check_hit_transport_scales.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All listed checks passed. CTest **1/1 passed**, 0 failures, 0 skips, 1.77 s
(total 1.83 s). The checker contains five successful solver runs and four
expected rejections. It uses periodic 4³ cells, pure N2 at 400 K/101325 Pa,
Mach 0.1, lambda 0.25 m, Re 80, Pr 0.72, WALE, SDC to 1e-8 s, and disabled
species diffusion.

An independent Python evaluation of NASA7 N2 cp and ideal-gas density/sound
speed matches both derived coefficients with **relative error 0** (gate 3e-14).
Explicitly specifying those coefficients and using the matched selected
mechanism produce byte-identical CSVs. Doubling Re halves both coefficients
exactly. Doubling `uin_norm` leaves both coefficients unchanged. The baseline
normalized Euler inventory residual is `8.0965343391571416e-16`, below 3e-12.
Re=0, simultaneous explicit coefficients, absent HIT velocity input and an
overflowing lambda value are rejected without the final CSV. No unexpected
failures occurred and no tolerances or golden data changed.

Logs: `build/hit-transport-scales-checks/{build.log,ctest.log,LastTest.log}`.
Inputs/CSVs/logs:
`build/amr-cvode-app-debug/hit-transport-scales/hit_scales_vi63wdv6/`.
The mechanism provenance check retains YAML SHA-256
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
