# Saved runtime-stencil extremum-slope comparison

The isolated scalar EP slope is not an adequate standalone replacement for
the current quartic candidate: it reduces sampled velocity-extremum errors,
but leaves the fine y=3/4 residual band unchanged. No evolution run or
production change is justified by this comparison alone. C2-C5 remain open.

`tools/check_eb_extremum_runtime.py` rebuilds and rechecks the 284 scalar
controls, reads the completed admission diagnostic log, and evaluates all
components, three directions and 25 lines per direction for each recorded
family. The 64 families yield 72,000 lines. All results are finite and child
pairs retain their parent averages within the checker's tolerance. This is
an isolated line reconstruction, not a full tensor/volume-weighted family,
EOS check, stage update, or fixed-PeleC comparison.

For the initial central y lines, recorded density is constant and recorded
momentum/density matches the declared 16-cell analytic shear. Comparison to
exact fine-cell velocity averages gives:

| Fine y | Conventional MC absolute error (m/s) | EP error (m/s) |
|---|---:|---:|
| 1/2 | 0.02748608 | 0.02748608 |
| 3/4 | 0.02330158 | 0.02330158 |
| 5/6 | 0.01556962 | 0.01556962 |
| 7/8 | 0.19091493 | 0.005467322 |

The sampled profile is fine y=1..16, parent y=1..8, at z=6 and 8;
it is not the full periodic domain. Its maximum EP error is 0.02748608 m/s,
versus MC 0.19091493 m/s. Do not compare this sampled maximum directly to
the previous full-domain scalar maximum. The original quartic local target
had error 0.000873185 m/s: replacing it with this linear EP reconstruction
would lose that improvement. These are interpolation estimates only.

Across all recorded lines, changes relative to componentwise MC appear
only along y. The final-fill x-direction momentum perturbations receive
exactly the same scalar slope as MC; this implementation has no whole-family
smoothness gate. That does not prove coupled EOS or multidimensional
continuity, and adding this slope behind the existing binary screen would
not remove that screen's switching.

Artifacts: `build/eb-extremum-runtime-20261023/report.json` (scalar controls)
and `runtime-report.json` (input SHA, binary SHA, changed-line counts,
central line outputs and analytic profile). Additional representation audits
compare momentum-to-velocity normalization and independently evaluated sine
inputs: no resolved slope discrepancy was found in the sampled central
initial lines; the normalization maximum is 4.44e-16.

```sh
python3 tools/check_eb_extremum_runtime.py \
  --input build/eb-inflection-admission-20261023/grid16.log \
  --work-dir build/eb-extremum-runtime-20261023
```

GNU Fortran 13.3.0 with bounds/FPE checks completed the replay. Python syntax,
mechanism provenance, project contract and whitespace checks pass. Next work
must retain higher-order interpolation away from extrema while constraining
its correction; merely replacing the full family with a linear EP slope
cannot address the measured residual band. No calculation is active.
