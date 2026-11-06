# Bounded quartic outer-sign roundoff candidate

The recorded second-step x-momentum stencil reproduces the limiter jump in
compiled Fortran: the frozen implementation gives slope 0.5977067070313679,
while the guarded candidate gives 0.8241920564945773. The upper child value
changes from 3.2831887317466024 to 3.339810069112405. This is scalar evidence,
not a flow improvement or fixed-PeleC comparison.

Only the experimental bounded-quartic path requests the new sign guard.
`extremum_slope` retains the original strict sign-product branch. For the
candidate, opposite signs count only when both differences exceed
8 * epsilon * max(abs(stencil)) in magnitude. This is a local floating-point
resolution guard, with an experimental safety factor of eight; it is not a
published limiter formula or a calibrated physical parameter. The magnitude
bound, family mean correction, EOS contraction and rollback are unchanged.
A large constant offset can change which differences are numerically resolved;
arbitrary-offset invariance is not claimed.

GNU Fortran 13.3.0, O0, runtime checks:

- `python3 tools/check_eb_bounded_roundoff.py --work-dir build/eb-bounded-roundoff-20261023/regression`:
  165 controls pass (three scales, three offsets, sign reversal, +/-4 ULP
  outer perturbations, and resolved opposite outer slopes), plus 284 scalar
  polynomial, step, tanh and sine controls.
- The original EP path separately passes its 284 scalar controls.
- `check_eb_bounded_tensor.py`: 960 signed-axis/mean families pass; maximum
  mean and covariance residuals both 2.220446049250313e-16.
- `check_eb_guarded_quartic.py` with the existing fixed-core registration
  build: 190 checks pass, including weighted conservation and EOS guards.
- Mechanism provenance, project contract, Python syntax and diff checks pass.

The bounded candidate's full-sine child errors at 8/16/32/64 cells are
0.0024483304, 8.7318536e-5, 2.8192714e-6 and 8.8822302e-8. These are isolated
interpolation errors, not temporal/spatial convergence of the flow solver.
Artifacts and compiler/source hashes are under
`build/eb-bounded-roundoff-20261023/`; `compiled-reproduction.json` retains
the frozen-versus-candidate scalar comparison.

Next: one matched 16-cell/100-ns flow comparison against the saved bounded
candidate, inspecting second-step pressure and transverse momentum as well
as all ordinary norms and conservation. A matched flow run is now active as described below; no
production path has changed. C2-C5 remain open.

The candidate flow executable built successfully by replacing only the frozen
slope object in the energy-observer link. The tensor, guarded family, tree
and observer objects are reused. Source hash and exact compile/link commands
are recorded in `build-record.json`. The run uses identical input bytes to
the energy observer, retains all diagnostics, and checks complete finite
output and initial-state equality. All 102 noninitial norm comparisons are
recorded, including any worsening; conservation must stay below the existing
1e-12 normalized budget limit. Runtime is bounded to 240 seconds and no
automatic retry is enabled. Results are pending. `handoff.json` records the
active handle and the downstream closure/energy checks.

## Completed 100-ns intervention

The run completed in 119.041907 seconds, exit zero. All six tables are
complete/finite and all 102 noninitial norm comparisons show no increase
against the unguarded bounded candidate. Normalized budget residual is
5.929347223590117e-14. Directional closure, 66 synchronization pairs and
eight fine-cell momentum trajectories pass.

| 100-ns Linf | Before | After | Ratio |
| --- | --- | --- | --- |
| y momentum | 4.772151124e-5 | 1.032886318e-6 | 0.0216440 |
| z momentum | 4.910872257e-5 | 1.088027005e-6 | 0.0221555 |
| pressure (Pa) | 1.191953541 | 0.013597119 | 0.0114074 |

At the second macro-step predictor, cut/full probe pressure offsets change
from -0.622653/-0.622555 Pa to -0.00581287/-0.00570520 Pa. Their held-energy
changes change from approximately -0.8925 to -0.00275741/-0.00270044.
`causal-comparison.json` preserves this comparison. With only the scalar
sign guard changed, this supports the outer-sign branch as a major cause
of the observed second-step error. It is not a proof that all remaining
interface error has the same cause.

A 16-cell/200-ns duration control now uses the same frozen binary, with a
360-second bound. It compares all 204 norms against the completed unguarded
bounded history. Its result is pending; `long-record.json`, `long.log` and
`handoff.json` identify this job. The first three output times must also
match the completed 100-ns candidate. Refined-grid duration convergence and
production/MPI integration remain unverified; C2-C5 remain open.

## Completed 200-ns duration control: remaining growth

The duration control completed in 253.603621 seconds, exit zero, with all
204 norms nonincreasing against the unguarded bounded candidate. Normalized
budget residual is 6.025620858862583e-14. `audit_long_prefix.py` passes exact
100-ns prefix equality for all four ordinary tables and identical binary
identity. This establishes duration/observer consistency, not acceptance of
the numerical candidate.

At 200 ns, y/z momentum Linf are 4.814352967336979e-5 and
5.088280331685412e-5 (ratios 0.125002 and 0.122229 to unguarded). Pressure
Linf is 1.1782779003697215 Pa (ratio 0.356120). The pressure history at
50/100/150/200 ns is 0.00681394/0.01359712/0.02048031/1.17827790 Pa:
the late jump remains and precludes claiming a durable resolution.

At 200 ns, pressure and y momentum maxima lie at fine-node (1,11,15), and
z momentum at (1,11,16); all are full cells adjacent to the coarse-fine
interface. The recorded j=6 donor's central line never crosses the guarded
negative outer-slope threshold (`long-limiter-lines.json`, ten distinct
lines). That observation does not diagnose j=11, which was not logged.
The epsilon guard fixes the demonstrated j=6 roundoff trigger but has not
established robust behavior across the entire trajectory. Next investigate
the new donor location and the limiter's continuity, without merely tuning
the epsilon multiplier upward or launching a refined grid prematurely.
No job remains active. Production adoption, C2 reference/convergence and
C3-C5 remain open.

## Compiled continuity audit: threshold guard is insufficient

`check_continuity.py` verifies the frozen scalar binary hash and sweeps the
outer sample through -32 to +32 ULP, both in the recorded orientation and
with the stencil reversed (130 samples). Between -15 and -14 ULP the input
changes by 4.440892098500626e-16 while the output slope jumps by
0.22648534946320997 in both orientations. `continuity-audit.json` preserves
all inputs/outputs. This independently demonstrates that the epsilon guard
moves the finite discontinuity instead of eliminating it; the earlier
+/-4 ULP controls were too narrow to establish continuity. Their reported
passes remain valid for their limited scope.

This scalar result is not a recorded j=11 evolved stencil and therefore
does not prove the cause of its 200-ns maximum. It does provide a reason to
reject the threshold guard as a complete robustness fix without another
flow run. Next evaluate a continuous combination of monotone and curvature
bounds against exact polynomial averages, resolved nonsmooth extrema,
step/tanh controls and perturbation sweeps. Preserve conservation/EOS guards
and avoid promoting a candidate solely for passing the first flow steps.
No new calculation was launched.
