# Continuous scalar envelope experiment

The new isolated `continuous_quartic_slope` removes the outer-sign threshold
switch. It takes the maximum of a monotone magnitude bound (zero when the
two central differences disagree in sign) and the existing curvature bound.
The quartic slope is then capped by that envelope and twice the centered
slope magnitude, with sign agreement required. The latter cap makes the
output vanish as the centered slope approaches zero. This is an experimental
construction, not the published EP algorithm or fixed-PeleC interpolation.
No production application selects it; the isolated tensor path is described below.

Independent compiled evidence (GNU Fortran 13.3.0, O0, runtime checks):

- `python3 tools/check_eb_continuous_quartic.py --work-dir build/eb-continuous-quartic-20261023/continuity`
  passes 284 polynomial, translated step, tanh and sine controls, plus 2,560
  perturbation inputs / 1,343 pairs. The latter include all 243 five-sample
  stencils with entries -1,0,1, perturbing each coordinate by +/-1e-10, and
  the measured donor's +/-32 ULP sweep in both orientations. Maximum observed
  slope-change/input-change ratio is 1.875. This finite sampling is not a
  global continuity or Lipschitz proof.
- Existing EP 284 outputs and guarded bounded 284 outputs remain byte-identical
  to their previous snapshots; the guarded variant's 165 controls also pass.
- Mechanism provenance, project contract, Python syntax and diff checks pass.

The maximum sine child errors at 8/16/32/64 cells are 0.0024483304,
8.7318536e-5, 2.8192714e-6 and 8.8822302e-8. These scalar results do not prove
multidimensional bounds, EOS admissibility, conservation under integration,
long-time convergence or physical validity. The monotone/curvature maximum
can differ from the former minimum bound near smooth extrema; nonsmooth
controls and family guards remain necessary before any flow experiment.

Exact commands, inputs, outputs, compiler and source/binary hashes are in
`build/eb-continuous-quartic-20261023/`. C2-C5 remain open.

## Isolated tensor/family integration

The isolated tensor now accepts an optional `continuous` selection, and
`guarded_children` an optional `continuous_correction` selection. Both default
false; the latter selects the bounded tensor path with the new scalar slope.
The mean recentering, weighted family handling and transactional EOS guards
are shared. No production application option has been introduced.

Both continuous and default configurations pass 960 tensor families and
190 family assertions. For continuous tensor, maximum mean residual is
2.220446049250313e-16 and signed-axis covariance residual is
4.440892098500626e-16. The family suite contains shared guard checks and
26 bounded-path-specific assertions; 190 is the suite total, not 190 new
continuous-specific cases. Commands use `check_eb_bounded_tensor.py` and
`check_eb_guarded_quartic.py` with `--continuous`; generated probe copies
and their hashes distinguish the tested path. Python syntax, project
contract and diff checks pass.

## Completed 16-cell duration experiment

The frozen isolated tree selects `continuous_correction`; the executable SHA256
is `3057dfabbcfea82391f174fa590040ae79aab44795e0bc44c75dbc8bbfd223ab`.
`python3 build/eb-continuous-quartic-20261023/run_long.py` completed in
247.470274 seconds, return code zero. Four connected-only tables contain
5/70/85/85 rows, with maximum budget residual 6.02560672052441e-14
(limit 1e-12). `summarize_long.py` verifies hashes and matches all 85 analytic
norm/location rows. These are not stage-flux tables.

At 200 ns, pressure Linf is 0.02740320390148554 Pa, versus
1.1782779003697215 Pa for the threshold-guarded candidate. Transverse
momentum y/z Linf are 3.881314363847576e-6 / 4.301895949030444e-6,
respectively 0.08062 / 0.08455 times the threshold-guarded errors.
The pressure Linf sequence at 50/100/150/200 ns is
0.00681393924 / 0.01359662256 / 0.02047952096 / 0.02740320390 Pa.

None of 204 noninitial norm comparisons worsens against the unguarded
bounded candidate. Against the threshold-guarded candidate, 36 comparisons
at 100/150 ns worsen slightly (largest ratio 1.0004056288973835).
The late jump is suppressed through this duration; this does not establish
long-time stability, production integration or fixed-PeleC parity.
`long-record.json` and `long-summary.json` retain every comparison.

## Completed 32-cell duration experiment

`python3 build/eb-continuous-quartic-20261023/run32.py` uses the same frozen
executable, root 32/fine 64, four 50-ns macrosteps to 200 ns. The 2400-second
bound permits one run, without automatic retry. Input SHA256 is
`cd62ef123f1c2a5bc0fed1a463debc32d118c0b0aec6e8db1aea5175be15421a`.
The run completed normally in 1790.8195885941386 seconds (29.85 minutes),
return code zero, within the 2400-second bound. All four tables are finite
and complete (5/70/85/85 rows). Maximum budget residual is
1.7663889050734788e-13, below the unchanged 1e-12 limit. Exact receipts,
commands and hashes are in `grid32-record.json`.

The exit observer `await_grid32.py` completed successfully and ran
`compare_grids.py`. The latter checks input hashes, common header controls,
same executable, matched times, table hashes, and analytic/location
agreement. All 204 noninitial norm comparisons are nonincreasing on
refinement from root 16 to 32 (zero-valued components may remain equal).
At 200 ns:

| Linf quantity | Root 16 | Root 32 | 32/16 |
| --- | ---: | ---: | ---: |
| y momentum | 3.881314363847576e-6 | 1.2852253214017234e-6 | 0.33113 |
| z momentum | 4.301895949030444e-6 | 1.1383091571196026e-6 | 0.26461 |
| pressure (Pa) | 0.02740320390148554 | 0.005044550460297614 | 0.18409 |

For y momentum, the ratio at 50/100/150/200 ns is
0.12185 / 0.16205 / 0.23809 / 0.33113; the gain weakens with duration.
At 200 ns the fine-grid maxima are node 2 cells (1,32,31), (1,2,32),
and (1,2,24), respectively. All three have volume fraction one and are
coarse-fine adjacent, not cut cells or z-wall-adjacent cells. The boundary
localization therefore persists despite the smaller errors.

`grid-comparison.json` retains all comparisons and selected locations.
The separately saved 50-ns snapshot agrees with all 51 norms of the old
32-grid bounded candidate; this early equality does not cover later times.
A two-grid comparison at fixed timestep is a refinement diagnostic, not
an asymptotic spatial-order proof. Neither candidate-to-candidate comparison
is fixed-PeleC or physical validation. The next discriminating check is
matched timestep refinement, starting with the cheaper 16-grid case, before
extending the spatial study or adopting the isolated candidate in production.

## Scalar branch review (exact arithmetic)

Source inspection of `continuous_quartic_slope` gives a stronger, narrowly
scoped result than finite perturbation sampling. For finite real inputs,
the scalar map is continuous and piecewise linear:

- The monotone bound switches only at `dm=0` or `dp=0`, where its minimum
  magnitude tends to zero.
- The curvature sign changes at `d2c=0`. Its magnitude is bounded by
  `1.875*abs(d2c)`, so both sides tend to zero.
- With curvature sign fixed, the choice of `2*abs(dm)` versus `2*abs(dp)`
  switches at `dc=0`. Exact arithmetic gives `dm+dp=2*dc`, hence the two
  caps agree at this boundary. At the intersection with `d2c=0`, the
  preceding zero bound also applies.
- The final sign agreement switches at `high=0` or `dc=0`. The output
  magnitude is bounded by both `abs(high)` and `2*abs(dc)`, so it vanishes
  at either boundary. Minima and maxima preserve continuity elsewhere.

A conservative global Lipschitz upper bound is 7.5 in the input infinity
norm, in exact arithmetic. Each linear piece selects zero or a signed
linear form: the quartic difference has coefficient absolute sum 50/32,
the centered cap 2, the one-sided caps 4, and the scaled second differences
at most `1.875*4=7.5`. A line segment crosses finitely many polyhedral
pieces; continuity permits summing their bounds along the segment.
The observed ratio 1.875 from the compiled sample is not the global bound.

This is an analytic source review, not a formal machine-checked proof or
a floating-point error bound. Rounded differences need not satisfy the
exact identities; overflow, tensor composition, EOS contraction and
long-time flow stability are outside this claim. The reviewed source hash
is saved in `scalar-continuity-review.json`. No production change follows
from this scalar argument alone.

## Completed temporal control

`python3 build/eb-continuous-quartic-20261023/run_halfdt.py` launches one
root-16 run to 200 ns with eight 25-ns macrosteps. Only the macrostep count
and timestep change from `long.dat`; the same frozen executable is used.
The run has a 720-second bound and no automatic retry. Expected table
counts are 9/126/153/153. The terminal audit verifies finite tables,
identical initial rows, unchanged budget tolerance, matching analytic and
location errors, and 204 norm comparisons at the common 50-ns output times.
Every increased norm is retained rather than rejected or hidden.

The purpose is to determine how much the remaining error depends on timestep
before interpreting the two-grid result as spatial convergence. Two timesteps
cannot establish temporal order. The run completed in 502.9038311149925 seconds, return code zero.
All four tables passed their count, finite-value and initial-state checks;
maximum budget residual is 5.923288273398557e-14. The full input comparison
confirms only line 5 (macrostep count and timestep) changed. Output hashes
were independently rechecked after completion.

At 200 ns, reducing dt from 50 to 25 ns changes Linf as follows:

| Quantity | dt 50 ns | dt 25 ns | 25/50 |
| --- | ---: | ---: | ---: |
| y momentum | 3.881314363847576e-6 | 3.950949164707118e-6 | 1.017941 |
| z momentum | 4.301895949030444e-6 | 4.367278391879482e-6 | 1.015199 |
| pressure (Pa) | 0.02740320390148554 | 0.027919867119635455 | 1.018854 |

83 of 204 matched-time norms increase; the complete list is retained in
`halfdt-record.json`. This is a completed diagnostic, not a successful
error-reduction intervention. Final maxima remain full coarse-fine-adjacent
cells: y (1,1,15), z (16,16,16), pressure (16,1,15), all on node 2.
The y maximum changes sign/location relative to the 50-ns run, so scalar
maximum ratios alone are not pointwise trajectory differences.

This control does not support timestep reduction as a remedy for the
remaining 16-grid boundary error. It does not prove time independence at
32 cells, exclude cancellation, or establish a temporal convergence order.
The next source diagnostic should target the continuous candidate's actual
maximum locations with stage pressure and directional flux capture;
previous stage captures concern the earlier bounded candidate and cannot
be used as direct evidence for the residual in this candidate.

## Completed stage audit at the remaining maxima

The isolated continuous-candidate stage run completed in 250.8000466004014 s
(return code zero), at root 16 and 200 ns with 50-ns macrosteps. Source base
is dirty main at `0a0e00ebe9138145fbb0e0e745665d0d768df009`; frozen source,
compiler/link commands and executable identity are in
`build/eb-continuous-stage-20261023/build-record.json`. Executable SHA256 is
`52365f13e240d536358e99fd430606b79076ff604eca0e5e741eb061906957a6`.
This changes probe coverage, not the continuous numerical candidate.

Re-running `audit.py` and `summarize.py` in that directory verifies all six
CSV hashes, exact parsed-row equality of the four ordinary tables with the
saved continuous run, 1152 stage face rows, 180 pre/post closure probe pairs
(2160 face rows), and 36 fine-probe momentum trajectories. All checks pass.
`maximum-origin.json` retains the normal-face reconstructed left/right
pressures at every sampled stage, alongside the integrated directional fluxes.

At the final full-volume maxima, the measured changes are:

| Probe / component | Observed change | x-face contribution | Normal-face contribution | Direct closure |
| --- | ---: | ---: | ---: | ---: |
| 12, (1,16,15), y momentum | -3.881314363847576e-6 | 8.080880274387145e-10 | -3.88212245e-6 (y) | 0 |
| 13, (16,1,16), z momentum | -4.301895949030444e-6 | 8.89428642327029e-10 | -4.30278544e-6 (z) | 0 |

The remaining cross-direction contributions are 3.01e-17 and 6.20e-14.
Integrating the recorded stage rates and applied closure increments gives
exactly the recorded changes at both probes in the Python audit. Summing
individual large pressure-containing face terms in a different order differs
by about 2.3e-15; this is not a physical extra contribution.

Already at the first predictor (25 ns before its corrector), probe 12 has
cell pressure offset +0.00354680694 Pa. Its lower y face has left/right
pressure offsets -0.00356206733/+0.00354680694 Pa; its upper face has
+0.00354680694/0 Pa. The corresponding normal momentum flux offsets are
-7.630224e-6 and +0.00177340348, giving a y-directed rate of -5.69930786.
Probe 13 has cell offset -0.00353123358 Pa; lower z-face pressures both
have that offset, while the upper face has -0.00353123358/0 Pa. Its normal
flux offsets are -0.00353123358 and -0.00176561679, giving rate -5.64997369.
The sampled initial transverse rates are only about -5.96e-8.

Thus the remaining transverse maxima accumulate predominantly through
normal spatial face fluxes after predictor pressure differences appear,
not through a directly applied closure increment at these cells. Normal
momentum flux includes advective and Riemann terms: reconstructed pressures
alone are not an exact pressure/advective decomposition. Prior synchronization
may still affect donor states. The next diagnostic must trace the first
predictor density/internal-energy balance and the interface-side halo donor
at these exact locations before selecting another interpolation change.
No production numerical policy is changed, and fixed-PeleC C2 parity,
MPI qualification of this candidate, and physical validation remain open.

## First-predictor density and internal-energy balance

The targeted diagnostic in `build/eb-continuous-predictor-20261023` completed
one 50-ns macrostep in 54.97695751301944 seconds, return code zero. It uses
the same numerical candidate with additional trace locations 12/13/14.
Executable SHA256:
`1c81d3a511d3e3dc8fe122f09170a93f87b73f517d4746e740684bc479d3b681`.
Build source and commands are retained in `build-record.json`; run input,
command and all six output hashes are in `record.json`. The four ordinary
CSV tables match the saved continuous-run prefix exactly (2/28/34/34 rows).
There are 288 stage-face and 540 synchronization-face rows.

`python3 build/eb-continuous-predictor-20261023/analyze.py` verifies all six
hashes, 12 energy-stage records, and six old-to-predictor density/held-energy
balances. The compiled trace records time-integrated mass and held-energy
face fluxes; division by fine spacing gives their update contributions.
All three sampled cells have full volume. The mass-flux balance passes
without a pressure-work source. The first fine interval is 25 ns:

| Probe | Density change | Held internal-energy density change | Predictor pressure offset (Pa) | Density-only pressure change (Pa) |
| --- | ---: | ---: | ---: | ---: |
| 12 | +1.327055809685618e-8 | +0.005093732092064 | +0.003546806940 | +0.003544792007 |
| 13 | -1.322734294317840e-8 | -0.005065619829111 | -0.003531233582 | -0.003533248504 |
| 14 | -1.322734294317840e-8 | -0.005065619829111 | -0.003531233582 | -0.003533248504 |

At this first interval both updates arise entirely from x-directed face
flux differences; y/z contributions are zero. Density state-minus-rate
residuals are at most 2.14e-17 and held-energy residuals at most 7.02e-12.
The held-rate minus recorded flux update is below 9.2e-16. Specific internal
energy increases by about 1.52157e-5 J/kg at all three probes.

The density-only column evaluates `p_old * delta_rho / rho_old`, holding
temperature and composition fixed in the ideal-mixture EOS. The observed
pressure increment differs by about +2.015e-6 Pa. This is a conditional
comparison, not a complete nonlinear EOS decomposition or proof that
composition is unchanged. It supports tracing the x-face mass flux and
its reconstructed halo states next; changing transverse reflux would not
address the measured first-predictor origin. Production numerics and the
fixed-PeleC qualification verdict remain unchanged.

## Initial x-face mass flux and halo velocity

An initial-operator-only probe completed in 13.558205898851156 s (exit zero)
in `build/eb-continuous-initial-20261023`. It uses the frozen continuous tree
and changes only the probe to evaluate the paired/connected policy once.
`build-record.json` retains source hash, compiler/link commands and binary
hash; `record.json` retains input/output hashes. No time evolution is run.
`analyze.py` validates the three table hashes and independently reproduces
the saved first-stage density rates at probes 12/13/14, then verifies them
against the two x-face mass fluxes at fine spacing 0.0003125 m.

Probe 12 has lower/upper x-face mass fluxes
0.3713754855275657 / 0.37120960355162164 kg m^-2 s^-1, giving density rate
+0.5308223230210628 kg m^-3 s^-1. Probes 13/14 have
0.3712096035516208 / 0.3713749453383526, giving -0.5290937175416843.
The frozen initial halo has x velocity 0.979470818621288 m/s outside the
fine patch versus 0.978597633264189 m/s inside: a difference of about
8.73185357e-4 m/s. Both outer halo cells along x have the same value;
the three sampled interior cells also agree. Initial pressure is 101325 Pa
on both sides to about 6e-11 Pa. The high-x side has the same velocity
mismatch, with reversed location of the ghost and interior states.

This ties the residual first-predictor density change to an initial
interpolation-induced streamwise velocity mismatch at the coarse/fine
interface. The recorded states are cell primitives, not reconstructed
Riemann left/right states; the face flux is measured, not inferred from
averaging these velocities. `initial-origin.json` retains the five-cell
lines for a subsequent reconstruction/donor replay. A numerical correction
still requires checking the coarse-to-fine cell-average interpolation and
its conservation/EOS constraints. Fixed-PeleC qualification remains open.

## Saved-donor replay: the residual is the unlimited quartic error

`python3 build/eb-continuous-initial-20261023/replay.py` compiles the scalar
continuous interpolant with GNU Fortran, O0 and runtime bounds checks. It
first requires byte identity with the frozen scalar source and validates
the saved halo table hash. It extracts the coarse y-lines at donors
(16,8,8) and (9,1,8), then applies the compiled half-cell interpolant.

For probes 12/13, replayed velocities are 0.9794708186212876 and
0.9794708186212832 m/s. Each agrees with the recorded halo within 7e-16.
The scalar slope equals the unlimited quartic expression in both cases:
-3.8356838751434235 and +3.835683875143426. The child pair preserves the
parent mean exactly in this replay. The analytic fine-cell average is
0.9785976332641891 / 0.9785976332641870 m/s, leaving errors
8.731853570984605e-4 / 8.731853570961290e-4 m/s.

Thus these initial residuals require neither limiter clipping nor a net
EOS/family adjustment to reproduce them. This scalar replay does not prove
that every family operation was inactive, or describe subsequent evolved
states. It shows that the unlimited quartic cell-average interpolant itself
accounts for the saved mismatch. Raising its order simply to force this
analytic fixture closer is not yet justified by fixed-PeleC evidence.

An independent scalar sine sweep at 8/16/32/64 coarse cells yields errors
0.02448330405850907, 0.0008731853570970172, 2.819271428672332e-5,
8.882230193851814e-7 m/s at the corresponding first child. These approach
fifth-order local interpolation error; they are not flow convergence orders.
`interpolation-replay.json` retains the donor lines, compile command, source
hash, replay and scalar results. Next compare the admissible numerical
policy and trajectory against the pinned reference before selecting another
numerical change. C2 remains unqualified.

## Fixed-reference reacting comparison: interrupted, not qualified

The isolated candidate was connected to the existing paired-reacting probe
with connected directions enabled, using the frozen continuous tree and
unchanged pinned TGReact capture. This is the existing root-8, 300-ns
reacting/transport contract, not a native reference of the analytic shear.
Artifacts are in `build/eb-continuous-reacting-20261023`.

The outer 600-second bound did not fire. The existing checker terminated its
probe at its internal 240-second limit; the checker returned 1 after
240.91349376365542 seconds. The initial wrapper label `completed_comparison`
only meant checks.json existed; inspection found TimeoutExpired and no
comparison steps. The receipt is corrected to `timed_out_in_checker`.
No 300-ns trajectory or H2-consumption verdict is available. No automatic
retry or increased limit was applied.

A complete 100-ns field table was retained. `partial.py` checks finite values,
1024 unique matching cells per table, time, geometry and composite-cell
selection, then converts the native CGS fields to SI before comparison.
It retains source hashes and historical same-contract 100-ns comparisons
in `partial-comparison.json`. At 100 ns the candidate differs from the fixed
reference by 0.161659017736838 K in temperature, 2.603295199538773e-4 in
relative density, 3.694910420879888e-4 in relative total-energy density, and
0.02823423414791973 kg m^-2 s^-1 in x momentum. These already exceed the
unchanged limits 1e-3 K, 1e-8, 1e-8 and 1e-6 respectively. This partial
field evidence rules out a passing full trajectory, but is not a completed
comparison or a physical-accuracy verdict. Candidate adoption is not
supported by the analytic-shear improvements. The next C2 task should
return to the much smaller default-policy discrepancy against the fixed
reference rather than extend this candidate run solely to obtain its endpoint.
