# Conservative mesh comparison, 2026-10-21

`tools/compare_tgreact_eb_meshes.py` now compares completed mesh diagnostics
at a common time step and endpoint. It separates differences in initial
sampling, final fields and evolution increments. The 8/16 pair passes.
The 32-cubed calculation also finished, but fails the unchanged same-grid
field comparison. The analyzer therefore refuses the three-grid order.

## Projection and independent check

For each coarse composite fluid cell, the analyzer collects the corresponding
eight cells on the same AMR level in the next resolution. Solid children are
excluded. Their fluid volumes must sum to the coarse fluid volume; projected
values use those volumes as weights. The projection checks its full-domain
integrals using absolute-integral scales, including fields with signed or
nearly cancelling integrals. Covered coarse cells are excluded at both sizes.
Stored temperature is volume averaged as an output diagnostic, not recovered
through an EOS from restricted conserved values.

The focused analytic check uses a constant density and a linear momentum
field evaluated at the exact fluid centroids of planar cut cells. Projection
reproduces the coarse fluid means and the exact domain fluid volume. This
would fail with equal weights in the cut row. Missing data are rejected.
The second focused check refuses mixed time steps and mechanism identities.
Both checks pass in 0.019 seconds. These are two focused checks of the new
analysis, not a repeated solver regression suite.

## Completed common-time-step pair

The additional 8-cubed run and the reused 16-cubed pilot both take two 50 ns
root steps to 100 ns, with native activation-basis chemistry and the native
held-EOS comparator. All same-grid comparison limits pass. The analyzer
loads complete initial and endpoint fields from each solver and independently
checks their physical times, level geometry and EB fractions.

| PeleF volume-weighted RMS difference, 8 versus restricted 16 | Initial | 100 ns | Change from initial |
| --- | ---: | ---: | ---: |
| Density (kg/m3) | 5.31227640e-5 | 5.34683765e-5 | 6.57232900e-6 |
| x-momentum (kg/m2/s) | 2.61810255e-2 | 2.48613322e-2 | 3.25503703e-3 |
| Total energy (J/m3) | 2.87652809e1 | 3.62544854e1 | 2.21046582e1 |
| Stored temperature (K) | 6.62986272e-1 | 6.63351802e-1 | 4.00978645e-2 |

The evolution-increment difference is the difference between each solution's
final-minus-initial field, after restriction. It is not an estimate of the
exact solution error. The initial temperature sampling difference dominates
the endpoint temperature difference, showing why an endpoint norm alone can
misrepresent the evolution discrepancy. Native held EOS gives essentially
the same pairwise differences. The largest scaled projection-integral
residual is `2.064e-16` across both solvers and both sampled times.

Two grids do not provide an observed order. Even a three-grid ratio at this
fixed time step will need interpretation for initial sampling, EB geometry
and temporal error; it cannot establish physical validation. The original
pinned spatial-EOS algorithm remains a separate failing reference gate.

## Completed 32-cubed comparison failure

Both solvers reached 100 ns: native held-EOS took 77.019 seconds and PeleF
1060.606 seconds. The checker reports twelve field exceedances: density,
three momenta, total energy and held energy at each of 50 and 100 ns.
Temperature, species, inventories and numerical energy budgets remain within
their original limits. This is a failed diagnostic comparison, not a timeout.

Saved initial composite fields match exactly. At 50 ns the largest total
energy difference is 0.224462655 J/m3 at fine cell `(2,1,28)`; only six
composite cells exceed one tenth of that maximum, all in fine indices
`i=1:2, j=1, k=26:29`. The y-momentum maximum is 7.17431058e-5 kg/m2/s at
`(2,1,27)`. This localizes the leading discrepancy near the lower periodic
corner, above the cut-cell row, but does not yet establish its cause.
`failure-localization.json` records both levels, both evolved snapshots,
the largest differences and source hashes. No solver rerun was needed for
this localization. A separately instrumented, intentionally truncated
first-spatial-stage capture follows; no tolerance is changed.

The completed four-evaluation capture in
`build/eb-mesh-stage-diagnosis-20261021/` covers three coarse evaluations and
the first fine evaluation. In the first fine evaluation, valid conserved
fields match exactly; the largest halo total-energy difference is
8.732e-11 J/m3. The largest combined hydro/transport energy-rate difference
is 0.003154 J/m3/s, corresponding to less than 8e-11 J/m3 over its 25 ns
step. Thus the large saved-field discrepancy is absent from the initial
fine spatial evaluation. Coarse correctors already show much smaller
differences (maximum input energy difference 1.674e-5 J/m3 in the third
evaluation); their later influence remains to be determined.

The first native instrumentation attempt incorrectly counted tiles and
then aborted during cleanup. It remains recorded as failed. The corrected
capture preserves original native tiling, counts complete function calls,
retains each tile's valid-rate ownership and exits intentionally with code
86 after flushing its fourth capture. The existing Fortran captures were
reused. Temperatures in the comparison use the primitive EOS temperature
on both sides, rather than comparing that value with native stored UTEMP.
These are truncated diagnostic captures, not completed evolution tests.
The first completion notification was accepted and resumed the analysis;
notification success did not override the failed calculation status.

The nine-evaluation continuation in `build/eb-mesh-stage-nine-20261021/`
completed: native took 51.926 seconds and PeleF 554.148 seconds. Both
instrumented executables exited with the expected code 86 and nine capture
sets. Its completion notification resumed analysis automatically.

The first large discrepancy appears in the halo entering the first fine
corrector (zero-based evaluation 4). Valid conserved fields still agree
to roundoff, but the halo y-momentum difference reaches 0.00101133 kg/m2/s.
The combined energy-rate difference then reaches 8.27216e6 J/m3/s at fine
`(1,1,28)`, creating the later interior-field discrepancy.

Complete eight-child ghost blocks recover each coarse mean and its limited
linear slope independently of the interpolation implementation. At the
periodic coarse donor `(1,32,14)` (unwrapped `j=0`), water density is about
4.4624e-9 kg/m3 and its centered z variation is only about 5e-17 kg/m3.
The native and Fortran water means differ by 5.70e-14 kg/m3. The water
profile sets the common slope factor to 0.314295 and 0.0877557 respectively,
which also multiplies the much larger momentum slope. Recovered means of
total energy differ by only 1.1e-8 J/m3. A second donor `(2,11,17)` likewise
has water-controlled factors 0.0298588 versus 0.0778989. A third examined
donor is limited by z-momentum. `ghost-blocks.json` and
`recovered-limiter.json` retain the recovered means, slopes and residuals.
The source limiter formulas agree; this evidence identifies amplification
of small input differences, rather than a different limiter formula.

The historical probe passes a mass-fraction absolute tolerance of 1e-14;
`react_sdc_state_3d` multiplies it by local density before calling CVODE.
Native uses 1e-18 g/cm3, or 1e-15 kg/m3, with typical-value scaling off.
At density near 0.1 kg/m3 the two effective tolerances are close, although
PeleF's varies with density. An earlier interpretation of the probe argument
as a density tolerance, and the resulting factor-ten claim, were incorrect.
Both effective tolerances exceed the first donor's tiny water variation.
This identifies sensitivity, not proof that a tolerance mismatch causes all
of the discrepancy.

V3 now passes the explicitly requested density tolerance through a callback
that bypasses the production driver's mass-fraction conversion. The first
V3 implementation omitted that bypass; its saved builds and controls are
historical intermediate diagnostics, not the final density-control semantics.
The final callback reports the effective CVODE tolerance range over both
initialization and evolution. V1/V2 keep their previous conversion.

The five-evaluation attempt in `build/eb-reaction-strict-prefix-20261021/`
failed: native captured all five stages in 992.929 seconds, but PeleF timed
out after 1800 seconds before its first spatial capture. That attempt used
the intermediate V3 behavior: the requested 1e-20 was multiplied by density,
whereas native used 1e-23 g/cm3. No paired strict result is claimed. Saved
native data remain useful: both water-controlled donor factors become zero
under its tighter settings, showing that even the original native factor
was sensitive to chemical integration accuracy.

A bounded isolated-cell pilot (`build/eb-reaction-cell-pilot-20261021/`)
then exercised the 81 exact distinct initial reaction inputs representing
33,792 fluid cells. With relative tolerance 1e-12 and a direct density
absolute tolerance 1e-18 kg/m3, all 81 solves completed within 5000 steps.
At 1e-20 none completed within that pilot step limit. Those are isolated
initialization solves, not a forced spatial trajectory or a claim of failure
at every possible step limit. All 324 attempts took 31.468 seconds.

The diagnostic probe now caches successful identical initialization inputs
by exact floating-point bits, with a 128-entry bound and no evolution cache.
An intermediate cache build reduced initialization to 81 solves and 33,711
hits: the 32-cubed default initial gas, held energy, temperature, reaction
history and budget matched saved bytes in 6.244 seconds. That build preceded
the density-control bypass; its moderate V3 timing is not a constant-density
tolerance result. Final callback/cache checks are recorded separately in
`build/eb-density-atol-20261021/`. Acceptance limits, production limiter,
canonical mechanism and historical V1/V2 defaults are unchanged. No
three-grid convergence claim is made.

The final callback/cache build passes the two V1/V2 one-step byte comparisons
and a V3 run that verifies its actual CVODE density tolerance is exactly
1e-15 kg/m3 throughout initialization and evolution. At 32 cubed, default
initialization takes 6.236 seconds, with 81 solves and 33,711 exact cache hits;
initial gas, held energy, temperature, reaction history and budget match the
previous complete capture byte-for-byte. The reported historical effective
density tolerance spans 9.58769e-16 to 1.03838e-15 kg/m3. Correct V3 controls
with rtol 1e-12 and constant density atol 1e-18 kg/m3 initialize in 6.437
seconds with the same cache counts. These timings cover initialization and
initial output only, not a complete timestep.

Four invalid probe control inputs and four invalid Python density controls
are refused. The two projection tests pass, including mixed-tolerance
rejection. Python syntax, diff and project-contract checks pass. Canonical
and activation-variant provenance/regeneration passed in the pinned Cantera
environment; the earlier system-Python attempt lacked Cantera and failed.

The bounded five-evaluation comparison in
`build/eb-reaction-moderate-prefix-20261021/` also timed out on PeleF after
900 seconds, retaining two coarse spatial evaluations. Native completed its
five evaluations in 291.783 seconds. Both used the verified 1e-18 kg/m3
density tolerance (native 1e-21 g/cm3), rtol 1e-12. Replaying those saved
operators isolated a forced-CVODE internal constraint stall; its fix,
analytic test and coupled follow-up are owned by the
[forced-CVODE record](eb-forced-cvode-20261021.md). Both earlier timeouts
remain failed attempts and are not counted as passing checks.

## Reproduction

Artifacts: `build/eb-mesh-convergence-20261021/`. `comparison8_16.json`
contains all field norms, projection residuals and input report/snapshot
hashes. The 16-cubed source and completed checks are in the
[pilot record](eb-mesh-pilot-20261021.md). `grid8-job.json` records the new
coarse run; `grid32-job.json`, `campaign.json` and `resume.json` identify the
completed calculation and failed comparison. The original 16-cubed run
was not repeated. Source provenance/regeneration was checked before launching.

Commands for the completed analysis:

```sh
python3 tests/python/test_compare_tgreact_eb_meshes.py
python3 tools/compare_tgreact_eb_meshes.py --runs <grid8-run> <grid16-run> --output comparison8_16.json
```

Python syntax, project contract and diff checks are recorded separately in
the campaign artifact directory. No full matrix, installed-binary audit,
MPI/restart or physical-validation campaign is part of this increment.
