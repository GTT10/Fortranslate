# Dynamic HIT reference and SI regularizers — 2026-09-27

Numerical change: `0a6ecf63`; comparison extension: `98a07af4`.
The full-reference setup and pinned source identities are unchanged from
[the hydro/transport/WALE comparison](hit-full-reference-20260927.md).

The checker now accepts `--transport --les dynamic`. PeleF selects the pinned
PeleC defaults: optimized three-point box test filter, test ratio 2, grid ratio
1, harmonic face interpolation. Dynamic coefficients determine the isotropic
term; the explicit PeleF `ci` control stays zero as required by its input API.

## SI correction

The reference's denominator regularizers have dimensions, unlike the
dimensionless coefficient floor. For SI grid callers, stress-squared,
stress and heat-flux-squared denominators now receive `1e-10`, `1e-9` and
`1e-6`, respectively. This follows stress_CGS = 10 stress_SI and
(rho u T)_CGS = 0.1 (rho u T)_SI. The coefficient floor remains `1e-8`.
The raw coefficient kernel retains its CGS default for existing reference
callers; the SI grid path opts into the converted denominators.

An independent manufactured CGS/SI coefficient equivalence check passes with
maximum absolute difference `6.66e-16` (gate `1e-13`). The first test invocation
incorrectly used a relative-error helper against zero and failed; it was
replaced by the intended absolute comparison without changing the gate.

## Completed checks and limits

GNU Fortran 13.3.0 Debug in the unchanged MPI/spray/SUNDIALS build:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_reactive_3d test_spray_core --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/amr-cvode-app-debug -R '^spray_les_dynamic_algebra$' --output-on-failure
python3 tools/check_hit_evolution_reference.py --pelec build/pelec-pinned-hit/Exec/RegTests/HIT/PeleC-HIT --fortran build/amr-cvode-app-debug/pelef_reactive_3d --exporter build/pelec-plotfile-export/pelef_plotfile_csv --velocity build/pelec-pinned-hit-source/Exec/RegTests/HIT/hit_ic_4_32.dat --work build/hit-full-reference --transport --les dynamic
```

Build targets were run separately; all completed successfully. Final CTest:
1/1 passed, zero skips, 0.05 s total. Provenance preflight, Python syntax,
project contract and diff whitespace checks also passed. Build/test logs are
`build/hit-dynamic-si-build.log`, `build/hit-dynamic-si-test-build.log`,
`build/hit-dynamic-si-test-failure.log` and `build/hit-dynamic-si-test-pass.log`.

Both pre-correction and post-correction full-reference runs passed the existing
`1e-10` gate. Maximum normalized difference is `7.67090568842832e-11`, in
temperature; momentum differences are at most `1.141878905816737e-13`.
The regularizer correction did **not** change the reported HIT residual.
Its cause remains unresolved; this result must not be described as the
roundoff-level agreement observed for the algebraic models.

Before: `build/hit-full-reference/evolution_5h3ig6qd/`.
After: `build/hit-full-reference/evolution_9m4m8ugc/`.
Both contain logs, inputs, reference plotfiles, CSVs and binary hashes.
Scope remains 4³ cells, one `1e-8` s step. Multiple-step dynamic parity and
long-time turbulence validation are not established.
