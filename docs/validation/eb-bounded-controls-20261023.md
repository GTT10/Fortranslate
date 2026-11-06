# Bounded tensor regular and duration controls

The regular-geometry 16-cell/50-ns solver completes in 98.668 seconds, exit
zero, with normalized budget residual 6.45868e-14. Compared with the recovered
vector-roundoff quartic regular candidate, all 51 final norms do not increase.
The y/z momentum and pressure Linf ratios are 0.0317692, 0.0317673 and
0.0317673; final maxima are 2.70956e-7, 2.77099e-7 and 0.00681334 Pa.
All four ordinary tables are complete/finite and the initial rows agree.

The original postprocessor selected the historical interrupted paired tables
and stopped on their incomplete row counts. That failed receipt is retained
as `regular-record.json`; no solver rerun was needed. `audit_regular.py`
uses the completed recovery in `build/eb-quartic-paired-only-20261023`, checks
identical input bytes, and records table hashes and full comparisons in
`regular-audit.json`. The runner now validates complete baseline tables
before launching a new calculation. This fixes evidence routing, without
changing numerical tolerances or turning the failed receipt into success.

The cut-geometry 16-cell/200-ns run completed in 252.734937 seconds,
exit zero, using the same frozen executable (SHA256
`1977e650311f9ffb10fb2b4dd272ed0fe29f1a4d7110b4382131244940a426f4`).
The comparator is the completed connected-direction paired-policy history,
not fixed PeleC and not the quartic 50-ns baseline. All 204 comparisons
(4 steps, 17 components, 3 norms) do not increase against that comparator.
The normalized budget residual is 5.929347223590117e-14. Complete finite
ordinary tables contain 5 native, 70 inventory, 85 analytic and 85 location
rows. `audit_long_prefix.py` independently passes exact CSV-row equality
for the initial and 50-ns outputs of all four tables against the completed
50-ns bounded candidate; duration extension preserves that prefix.

| Time (ns) | y momentum Linf | z momentum Linf | Pressure Linf (Pa) |
| --- | --- | --- | --- |
| 50 | 2.709559593e-7 | 2.770992566e-7 | 0.00681393924 |
| 100 | 4.772151124e-5 | 4.910872257e-5 | 1.191953541 |
| 150 | 1.799609028e-4 | 1.903617999e-4 | 2.290467329 |
| 200 | 3.851432217e-4 | 4.162921281e-4 | 3.308649543 |

At 200 ns the corresponding ratios to the connected-direction comparator
are 0.189242, 0.449047 and 0.581960. Improvement persists, but the very
small first-step error does not persist through later steps. This duration
control therefore does not justify extending the two-grid/50-ns convergence
claim to 200 ns or promoting the candidate into production. The next causal
check is the 50-to-100-ns transition: identify maximum locations and inspect
pressure reconstruction, limiter admission and directional flux increments
at those locations before changing the method. No larger-grid run is needed
merely to reproduce the favorable first-step slope.

Artifacts: `build/eb-bounded-controls-20261023/`; `long16-record.json`
records the completed run and all comparisons; `long16-prefix-audit.json`
records the independent prefix check. `handoff.json` has no active job.
No production option or MPI path changes. C2-C5 remain open.

## Duration maximum localization

`localize_duration.py` verifies the saved table hashes and matches all 85
location rows to the analytic Linf values and times. All 12 selected maxima
(y/z momentum and pressure at four noninitial times) are coarse-fine adjacent.
At 50 ns they have fine-node i=1, j=1. At 100, 150 and 200 ns, the y momentum
and pressure maxima are at node 2, (1,6,11), volume fraction 0.44; the z
momentum maximum is at node 2, (1,6,16), volume fraction 1. The pressure
maximum changes sign between 50 and 100 ns. These are maxima at different
locations, not a tracked single-cell trajectory.

The existing stage observer samples j=1, 8 and 9 on this grid and therefore
cannot explain the new j=6 maxima. The next observer run should cover the
two identified j=6 cells through 100 ns, retain original diagnostic points,
and verify the ordinary-output prefix against the completed duration run.
The saved results identify where to measure, but do not yet distinguish
limiter switching, reconstructed pressure and synchronization as causes.
Detailed localization: `duration-localization.json` in the control directory.

## Targeted second-step observer (completed)

The isolated paired probe now adds IDs 10/11 for fine-node cells (1,6,11)
and (1,6,16), retaining IDs 1:9. Only its stage/sync observer is changed;
the bounded evolution solver and interpolation objects are reused. The
frozen probe source, compile/link commands and hashes are recorded in
`build/eb-bounded-transition-20261023/build-record.json`. Binary SHA256:
`2ff35f77ddce4ab4a5a4896ac0dd73b1b365c8551d305f7812cfe8e9b1c8d68f`. The probe builds successfully.

The matched 16-cell input changes only macro-step count from four to two.
It completed in 115.869214 seconds, exit zero. All four ordinary tables agree
exactly with the saved duration history through 100 ns. Stage/sync tables
contain 384/792 finite rows. `audit.py` passes directional momentum closure,
66 synchronization pairs, and eight fine-cell momentum trajectories.
Mechanism provenance, Python syntax and project contract checks also pass.
The diagnostic has not changed production behavior.

At probe 10 (cut cell), y momentum changes by 4.77215112410843e-5;
integrated spatial rates give 4.77215112399002e-5 and measured synchronization
adds only 1.18405e-15. At probe 11 (full cell), z momentum changes by
-4.91087225670573e-5, exactly reconstructed from the sampled spatial rates,
with zero measured synchronization increment. This excludes a direct
synchronization increment as the dominant source at these two locations
for this trajectory; it does not exclude effects of earlier synchronized
states on later interpolation.

The second macro-step old-state cell pressure offsets are -0.00390267 and
-0.00383475 Pa at probes 10/11. Their predictor offsets already reach
-0.62265330 and -0.62255536 Pa. The corresponding y-direction contributions
to y momentum rate increase from about 3.65/3.58 to 989.95/989.86. Probe 11's
z-direction contribution to z momentum rate changes from -6.20 to -996.30.
Thus the large transverse forcing is present in the predictor before the
subsequent synchronization. Reconstructed face pressures and aperture flags
are retained in `transition-summary.json`; zero-aperture face placeholders
must not be interpreted as physical zero-pressure faces.

Next: trace the second-step predictor's energy and halo interpolation at
j=6 to distinguish a changed interpolation/limiter decision from the
subsequent pressure reconstruction. The current trace establishes when
and through which flux directions momentum changes, but does not yet
identify the upstream cause. There is no active calculation. C2-C5 remain
open; this is not fixed-PeleC agreement or physical validation.

## Energy observer follow-up (completed)

An isolated tree copy now logs held internal energy, its spatial rate, gas
state/rate and six integrated scalar face fluxes at the two target cells.
It also relocates the existing interpolation value/family log to j=6.
No numerical expressions or limiter choices were changed. The diagnostic
build completed; a two-step/100-ns run has started with a 240-second bound.
The runner requires all six CSV tables to be byte-identical to the completed
transition observer, in addition to the ordinary duration-prefix check.
`analyze.py` is prepared to compare predictor energy changes to their old-state
rates; those checks are reported below. Artifacts and the active handle are under
`build/eb-bounded-energy-20261023/`, with source/command hashes in
`build-record.json` and resume instructions in `handoff.json`.

The energy-observer run completed in 116.317453 seconds, exit zero; all six
CSV tables are byte-identical to the transition observer. Directional closure,
66 sync pairs and eight momentum trajectories pass again. Sixteen energy
records are complete. For the first fine predictor of macro-step 2, held
energy changes by -0.892518770823 (cut) and -0.892475827946 (full), matching
the measured old-state held rates times dt to absolute residuals below
1.01e-11. Across all eight predictor comparisons the largest residual is
1.25e-11. Total-energy changes differ from the logged raw gas rate prediction;
they must not be conflated with the held-energy update or claimed as a new
conservation failure without auditing the existing energy reset terms.

The logged interpolation families remain admitted with theta=1. The central
y line of x momentum reveals an upstream limiter discontinuity: initially
its outer difference is exactly zero, while after evolution it is
-1.3322676295501878e-15. This activates the strict sign-product extremum
branch in `slope_bound`. The quartic slope changes from 0.8241815310 to a
bound of 0.5977067070, despite the unbounded slope remaining 0.8241920565.
The scalar child value changes from 3.3398163986 to 3.2831887317; the logged
full-cell tensor halo similarly changes to 3.2831888919. The scalar replay
is not a full tensor replay and does not alone prove the complete causal
chain, but identifies a specific discontinuous branch to test before
another flow intervention. Artifacts: `energy-analysis.json` and
`limiter-analysis.json`. Next: compiled scalar reproduction and a principled
roundoff-robust branch criterion, with scale/translation and discontinuity
controls. No active calculation; production behavior remains unchanged.
