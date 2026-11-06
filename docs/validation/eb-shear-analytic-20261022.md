# Analytic variable-velocity EB-AMR shear test, 2026-10-22

This test evaluates the halo policies against a steady inviscid Euler
solution, independently of the pinned PeleC interpolation policy. The
reacting reference comparison and its failed field gates remain unchanged.

## Analytic solution and finite-volume initialization

Pure N2 at 900 K and 101325 Pa has constant density, composition and held
internal-energy density H. Velocity is u(y)=10 sin(2*pi*y/L), v=w=0,
where L=0.01 m. All fields are independent of x; x/y are periodic and the
z slip walls are tangent to the flow. Chemistry, molecular transport and
forcing are disabled. Euler advection and pressure divergence vanish, so
all fields retain their initial values at every physical time.

For a cell centered at phase theta with phase width delta:

- mean(u) = 10 sin(theta) sinc(delta/2)
- mean(u^2) = 50 [1 - cos(2 theta) sinc(delta)]
- mean(rho E) = H + rho mean(u^2)/2

Here sinc(x)=sin(x)/x. Momentum uses mean(u), while total energy uses
mean(u^2). The z cut plane does not alter these y averages. A bottom-up
restriction check verifies the analytic coarse averages to 1e-12 scaled
precision before any time integration. The solution does not import a
computed fine-grid or native-solver trajectory as its truth.

The distinction between mean(u^2) and mean(u)^2 matters: the existing held
energy closure can remove unresolved cell kinetic variance. This test
retains the true analytic total-energy average as its target and records
returned numerical-energy corrections separately. A corrected global budget
does not turn a total-energy field discrepancy into zero analytic error.

## Numerical comparison

The native common probe accepts `--audit-shear`, initializes the analytic
state, and runs the ordinary recursive held-MOL integrator for four 50 ns
root steps to 200 ns, with two MOL iterations. Four policies are compared:
default, separate momentum, plus consistent velocity, plus paired pressure.
Each resolution has two nested ratio-two patches with 8, 16 or 32 cells per
axis. The physical z plane remains at 0.0033 m; cut fractions vary with mesh.

`tools/check_eb_shear_reference.py` drives the comparison using V4 geometry
imports, freezes the binary, records commands and identities, checks every
trajectory and the unchanged 1e-12 normalized component-budget threshold,
and writes the observed orders when all three meshes finish. The 32-grid
run is admitted only if the 16-grid measured time times twelve (8x cells
and 50% margin) is no more than 3600 seconds. Each grid has the same hard
timeout; completed data are retained on failure or cost refusal.

Per-step analytic CSVs contain absolute volume-weighted L1, RMS L2 and Linf
errors over exposed fluid cells, excluding covered coarse and solid cells.
Components 1–14 are gas conserved fields in the LiDryer order, 15 is held
energy, 16 is EOS temperature, and 17 is EOS pressure. Conserved and held
errors use the exact analytic initial averages. Initial, intermediate and
endpoint errors are retained, including zero-error and roundoff fields;
orders of roundoff-sized errors must not be interpreted as accuracy orders.

## Completed results

All 12 trajectories (three meshes times four policies) completed four root
steps to 200 ns with zero execution failures, skipped grids or timeouts.
Grid runtimes were 67.8, 203.9 and 927.7 seconds, about 20 minutes total.
The process-exit watcher successfully queued this thread's completion message
once; `completion-notification.json` records the successful dispatch.
No job remains to resume, and no calculation was relaunched for this audit.

The table below gives selected RMS L2 errors. All 17 components and all
three norms for every policy/resolution are retained in `endpoint-norms.csv`;
`all-field-results.md` lists all nonzero L2 fields and their observed orders.

| Field | Policy | 8 | 16 | 32 | Order 16→32 |
| --- | --- | ---: | ---: | ---: | ---: |
| Density (kg/m3) | Default | 9.48580e-6 | 5.78728e-6 | 2.98532e-6 | 0.955 |
| Density (kg/m3) | Paired | 5.41318e-6 | 3.14959e-6 | 8.31902e-7 | 1.921 |
| x momentum (kg/m2/s) | Default | 6.11541e-3 | 3.95399e-3 | 2.23070e-3 | 0.826 |
| x momentum (kg/m2/s) | Paired | 3.50630e-3 | 2.21530e-3 | 7.00157e-4 | 1.662 |
| y momentum (kg/m2/s) | Default | 2.09927e-4 | 1.70732e-4 | 1.26812e-4 | 0.429 |
| y momentum (kg/m2/s) | Paired | 1.42420e-4 | 1.00705e-4 | 4.02822e-5 | 1.322 |
| z momentum (kg/m2/s) | Default | 1.26795e-4 | 1.10086e-4 | 8.62490e-5 | 0.352 |
| z momentum (kg/m2/s) | Paired | 1.46263e-4 | 1.65505e-4 | 1.30338e-4 | 0.345 |
| Total energy (J/m3) | Default | 3.68267 | 2.22420 | 1.14646 | 0.956 |
| Total energy (J/m3) | Paired | 2.15259 | 1.21752 | 0.321613 | 1.921 |
| Held energy (J/m3) | Default | 3.63854 | 2.21984 | 1.14579 | 0.954 |
| Held energy (J/m3) | Paired | 2.07674 | 1.20864 | 0.319476 | 1.920 |
| Temperature (K) | Default | 7.00822e-5 | 4.29147e-5 | 2.39160e-5 | 0.843 |
| Temperature (K) | Paired | 5.39396e-5 | 2.74774e-5 | 8.67941e-6 | 1.663 |
| Pressure (Pa) | Default | 2.53439 | 1.54622 | 0.797848 | 0.955 |
| Pressure (Pa) | Paired | 1.44640 | 0.841682 | 0.222396 | 1.920 |

N2 density error equals total density error. The other eight species remain
exactly zero at all reported times, with zero norms; their order is undefined.
Separate momentum alone, plus consistent velocity, and plus paired pressure
have nearly identical errors in this constant-density control. For example,
at 32 the pressure L2 errors are 0.22239503, 0.22239503 and 0.22239558 Pa.
Most of the improvement over default comes from separate momentum; this
case does not show an additional pressure-policy benefit.

## Limits revealed by the full norm audit

At 32 the paired policy reduces x-momentum L2 error to 31.4% of default,
pressure to 27.9%, and temperature to 36.3%. However, z-momentum L2 error
is 151.1% of default and is nonmonotone from 8 to 16. L1 decreases for every
nonzero field under all policies, but that is not sufficient for all-field
convergence.

Maximum-norm transverse momentum errors increase on both refinements. For
the paired policy, y-momentum Linf is 0.00195654 → 0.00383839 → 0.00491369
kg/m2/s; z-momentum Linf is 0.00222812 → 0.00440807 → 0.00753186 kg/m2/s.
The default shows the same increasing transverse-Linf trend. Several other
Linf fields are nonmonotone, including paired pressure (22.8926 → 26.8917
→ 18.2824 Pa). Do not claim uniform convergence or general second-order
accuracy from the improved integrated norms.

This is a fixed-dt, short-time three-grid experiment; cut fractions change
with resolution. It does not separate temporal, boundary, interface and
small-cell contributions to the remaining local error. The next useful
investigation is spatial localization of the transverse maxima followed by
a matched boundary/interface or timestep discriminator. These aggregate
norm files alone cannot identify the worst cell or its causal operator.

## Conservation and energy closure

All cases pass the unchanged 1e-12 normalized component-budget threshold.
The largest residual is 2.027e-13, in the default 32 case; paired maxima
are 1.604e-14, 5.954e-14 and 1.427e-13 at 8/16/32. The budget includes
returned signed numerical-energy corrections and is not strict raw-energy
conservation. Paired corrections are -3.05106e-7, -7.74218e-8 and
-1.94203e-8 J, approximately quartering at each refinement.

Independent integration of the analytic initial subcell kinetic variance
shows that the magnitude of the final correction is 1.00016–1.00221 times
that variance across all 12 runs. This is consistent with the existing
held-energy closure removing unresolved kinetic variance. It is not a
cellwise proof that every correction comes from that operation, and the
analytic total-energy error remains reported without subtracting it.

## Verification and retained evidence

`build/eb-shear-analytic-20261022/runs/record.json` contains complete commands,
input/executable hashes, elapsed times and all endpoint norms. `audit.json`
records an independent recomputation of 840 component-inventory residuals,
checks of all 1020 norm rows (complete component/time coverage, finite and
ordered norms, initial errors and endpoint identity), and recomputation of
the observed L2 orders. Frozen binary and input hashes match. Norm ordering
and report checks do not independently recalculate cellwise norms: full
cellwise fields were not exported by this run.

The prior native pressure diagnostic was rerun for one full root step with
both existing policies before this batch. All four summary/inventory CSVs
matched the preceding frozen executable byte for byte (`parity.json`).
Mechanism provenance, Python syntax, project contract and whitespace checks
passed again during this final audit. Debug array-temporary warnings in the
solver logs are retained; no trajectory was interrupted.

The completed verdict is **bounded integrated-error improvement, with
unresolved transverse/local-error convergence**. Original reacting
fixed-reference failures remain open. No production defaults, tolerances
or golden hashes were changed and no further calculation was launched.

## Follow-up localization discriminator

The 16-cubed follow-up in `build/eb-shear-localization-20261022/record.json`
completed cut50 and cut25 (four policies each, 200 ns endpoint), then failed
regular50 during geometry initialization. The original failed record and logs
are preserved. Its regular50 entry retains a stale `running` label, but
returncode=1 and the top-level failed status are definitive; no regular
trajectory completed in that attempt.

The diagnostic adds maximum-error coordinates without changing the numerical
method. Ties report the first maximum in traversal order. All 12 preexisting
cut50 summary, inventory and norm CSVs are byte-identical to the earlier
16-grid result. `audit_completed.py` independently recomputed 784 component
inventory residuals and checked 952 norm/location pairs, including complete
time/component coverage, finite ordered norms and signed-error/Linf identity.
This does not independently recompute cellwise norms from full fields.
`audit-completed.json` records the result and source/artifact hashes.

At 200 ns the paired-policy results are:

| Observable | cut50 Linf | cut25 Linf | Change on halving dt | cut50 maximum node / (i,j,k) |
| --- | ---: | ---: | ---: | --- |
| y momentum (kg/m2/s) | 0.0038383873 | 0.0038284471 | -0.259% | 2 / (1,1,11) |
| z momentum (kg/m2/s) | 0.0044080717 | 0.0044045550 | -0.0798% | 2 / (1,16,12) |
| pressure (Pa) | 26.8917212 | 26.9126576 | +0.0779% | 2 / (1,1,11) |

All these maxima adjoin a coarse/fine boundary. The y-momentum and pressure
maxima occur in a cut cell (volume fraction 0.44); the z-momentum maximum
occurs in a full-volume cell. At cut25 the z maximum is at (1,1,12), with the
same geometric classification. Default transverse maxima are full-volume
cells (1,1,15) and (1,1,16); its pressure maximum is a cut cell. Momentum-only
and consistent-velocity policies have the same classification as paired.
No reported maximum above adjoins a domain-z wall. These flags may overlap
and identify geometry, not a causal operator. They do not directly locate
the earlier 32-grid maxima. The small response to dt halving does not support
time discretization as the dominant source of these 16-grid maxima.

All eight completed trajectories pass the unchanged 1e-12 normalized
component-budget limit; the independently recomputed maximum is
6.12932096524e-14 (cut25 default). This includes the signed numerical-energy
correction; raw total energy is not strictly conserved. All 17 component
norms and maxima remain available in the original CSVs and record.

### Regular initialization failure and bounded recovery

`regular50.log` reports `ERROR STOP Invalid padded EB geometry` before any
time integration (0.42 s). The input placed the EB at z=0. Although this is
the physical lower domain boundary, it lies inside the padded geometry at an
exact grid face. `build_axis_plane_eb_geometry_3d` explicitly rejects an
interior face-aligned plane because no positive-volume cut cell can own its
wall metric. Thus this is an inadmissible diagnostic input, not observed
solver instability. The geometry guard and all solver tolerances remain
unchanged.

Only the failed regular condition was restarted in
`build/eb-shear-localization-20261022/regular-recovery/`, with plane=-0.01 m,
below the padded lower bound -0.00375 m. This selects fully regular geometry
through the existing geometry constructor; the domain slip wall stays at
z=0. The frozen executable SHA256 is
`9e04d48691d729e0f6ba88e3c5273bff02d5fcad3861d040619461d02986ad4c`.
The recovery completed four 50 ns root steps and the same four policies and
1200-second timeout. It does not rerun either completed cut condition.
Driver PID 1865417 and child PID 1865423 were recorded at launch;
`regular-recovery/record.json` owns subsequent status and exact commands.

The original failure notification was successfully queued to this conversation
(returncode 0). The recovery's separate `completion-notification.json` was
verified queued with returncode 0 after terminal status was saved. The
completion message triggered this final audit. No further calculation or
restart is needed for this comparison.

The regular case changes fluid volume and wall position relative to the cut
case (embedded wall at 0.0033 m). It is a geometry ablation, not an identical-
domain reference, and cannot alone isolate a single boundary operator.
Transverse/local convergence and physical qualification remain open.

Follow-up record checks: `python3 tools/check_mechanism_provenance.py --bundle
mechanisms/h2o2_full.json`, `python3 tools/check_project_contract.py --source .`,
Python compilation of the completed-case auditor and recovery runner, and
`git diff --check` passed. Two earlier provenance invocations were command
errors (missing required bundle; then a directory supplied as source); the
corrected invocation above checked the pinned YAML hash successfully. No
additional solver regression suite was run for this evidence-only update.

### Completed regular comparison

Recovery completed in 421.76 s with returncode 0 and all four policies reaching
200 ns. The frozen binary hash matches the cut runs. The retained log contains
Debug array-temporary warnings and IEEE underflow/denormal flags, but no fatal
error; all exported norms and inventories passed finite-value checks. The
analytic initial restriction residual is 8.881784197e-16.

| Observable | paired cut50 Linf | paired regular50 Linf | Change | regular maximum node / (i,j,k) |
| --- | ---: | ---: | ---: | --- |
| y momentum (kg/m2/s) | 0.0038383873 | 0.0020351907 | -46.98% | 2 / (1,8,15) |
| z momentum (kg/m2/s) | 0.0044080717 | 0.0009270558 | -78.97% | 2 / (1,8,16) |
| pressure (Pa) | 26.8917212 | 5.6852649 | -78.86% | 2 / (1,8,1) |

All three regular maxima are full-volume cells adjacent to a coarse/fine
boundary. The pressure maximum also adjoins the lower domain slip wall.
Momentum-only and consistent-velocity policies show nearly identical
reductions. Default maxima instead change by +0.00847%, +0.000540% and
-0.3049% respectively, and remain coarse/fine-adjacent full cells at node 2,
(1,1,14), (1,1,16), and (1,1,13). Thus removing cut geometry does not remove
the local error. The candidate's large reduction depends on geometry, but
the changed fluid volume/wall position prevents attributing it solely to
the cut-cell algorithm.

Paired regular L2 values for y momentum, z momentum and pressure are
1.005252878e-4, 1.412311017e-5 and 0.342132936 respectively (cut50:
1.007051001e-4, 1.655047134e-4 and 0.841682084). In particular, the y maximum
decreases substantially while its volume-weighted L2 barely changes;
maximum-error and integrated-error improvements must remain separate.
The completed evidence supports geometry-sensitive local errors and a
remaining coarse/fine contribution, not demonstrated local convergence or
a unique responsible operator. Time halving on cut geometry changes paired
maxima by less than 0.3%; no regular time-halving experiment was run.

`regular-recovery/audit_completed.py` independently recomputed 280 inventory
residuals and checked 340 norm/location pairs plus frozen binary/input hashes.
Together with the cut audit this covers 1,064 inventory rows and 1,292
norm/location pairs. All 12 trajectories pass the original 1e-12 corrected
budget threshold; the maximum is 8.14554310426e-14 (regular paired).
The signed numerical-energy ledger remains included, so this is not a claim
of uncorrected total-energy conservation. `comparison-all-fields.csv` records
all 17 endpoint components, all three norms and maximum locations for each
of the 12 trajectories; `final.json` closes the batch while preserving the
original failure history. No solver was started during this final audit.

### Initial spatial-operator follow-up

The diagnostic-only `--audit-shear-initial` probe mode initializes the same
finite-volume analytic shear and evaluates `compute_mol_eb_tree_hydro_3d`
once for each of default/momentum/velocity/paired. It exports prepared halo
pressure, temperature, held energy and all velocity components; face mass,
energy and all three momentum fluxes; and all conserved-component rates on
valid composite cells. Face fluxes are divided by dt. The rate output excludes
covered parent cells. This is a spatial-operator diagnostic before chemistry,
time integration, energy reset and reflux; its rates are not a full recursive
step or a conservation qualification. Existing audit output formats are
preserved outside this new mode. Production solver code is unchanged.

The matched cut/regular inputs from the completed 16-grid comparison are
complete in `build/eb-shear-initial-20261022/`; `record.json` retains exact
commands, frozen executable/source/input hashes and child PIDs. Driver PID
1897708 runs one condition at a time, with a 600-second bound per condition.
A finalization handler queues a completion message to this conversation on
success or failure, with dispatch evidence in `completion-notification.json`.
No time-evolution trajectory is repeated. The purpose is to test whether
transverse rates already arise in the initial spatial operator and to compare
halo/face observations before interpreting later maxima. No result is claimed
until the saved outputs have been audited.

The probe target built successfully in the existing GNU/OpenMPI registration
build. Pinned-mechanism provenance, project contract and whitespace checks
passed. This does not substitute for the pending diagnostic-output audit.

The existing `--audit-shear` mode was compared before/after this probe change
at zero root steps: all 16 summary/inventory/norm/location CSVs are byte-identical
(`build/eb-shear-initial-20261022/parity.json`). This checks initialization and
output compatibility only, not evolved-trajectory parity. Python syntax checks
for the driver and parity script passed. Completion notification was verified
armed after launch.

### Completed initial spatial-operator audit

Both initial-only conditions completed (cut 17.08 s; regular 31.70 s), each
with all four policies and returncode 0. The completion notification was
queued successfully. No simulation was launched during analysis.
`build/eb-shear-initial-20261022/analyze.py` verifies binary/input hashes,
finite values and unique row keys across 1,138,496 exported rows; results and
CSV hashes are in `analysis.json`. The frozen executable SHA256 is
`e0dea32753a44c6671e1f2ece65e6bcb02c28b131a99260e8025c815c98d648d`.
Logs retain array-temporary warnings and IEEE underflow/denormal flags; no
fatal marker was found. Restriction residual is 8.881784197e-16.

**The initial transverse operator is balanced in all eight cases.** All
prepared halo v and w values are exactly zero in the export. Maximum halo
pressure departure from 101325 Pa is 1.455191523e-11 Pa. On open faces,
y/z momentum fluxes equal the corresponding constant-pressure diagonal
flux or zero off-diagonal flux, to at most 1.455191523e-11 Pa. Both transverse
momentum rates are exactly zero on every exported valid composite cell.
Their zero maxima have no meaningful location; traversal-first cells must
not be interpreted as error hot spots. Area-weighted coarse/fine transverse
momentum flux mismatches also vanish in this initial-only comparison.

In contrast, the streamwise state and flux are not stationary:

| Geometry / policies | Maximum padded-halo u error (m/s) | Maximum density rate (kg/m3/s) | Maximum x-momentum rate (kg/m2/s2) | Maximum energy rate (J/m3/s) |
| --- | ---: | ---: | ---: | ---: |
| cut / all four | 0.959794154 | 584.869712 | 352795.355 | 3.81320789e8 |
| regular / default | 0.959794154 | 584.869712 | 352795.355 | 3.81320789e8 |
| regular / momentum, velocity, paired | 0.190914927 | 117.770434 | 71868.394 | 7.73134881e7 |

The u comparison uses the exact finite-volume sine average at the exported
cell's level and y index, including periodic halo coordinates. These are
maxima over the entire padded halo, not just the nearest active stencil.
The nonzero x-momentum rate maxima, however, are valid fine cells: node 2
(1,2,11) for all cut policies, (1,2,1) for regular default, and (1,8,1) for
regular candidates. They lie beside the coarse/fine interface. All cases
have nonzero density, x-momentum, total-energy and N2 rates; absent-species
rates remain zero. The analytic stationary shear has zero exact rates.

Maximum coarse/fine flux discrepancies below are positive-coordinate fine
aperture-weighted 2x2 averages minus the matching coarse aperture-weighted
flux, per coarse full-face area. The x/y low interfaces wrap to the root's
high periodic face; the high interfaces and z high interface use the
matching parent internal face. The physical z wall is excluded.

| Geometry / policy | Mass (kg/m2/s) | x momentum (Pa) | Energy (W/m2) |
| --- | ---: | ---: | ---: |
| cut / default | 1.45525064e-4 | 0.244064795 | 147.370223 |
| cut / momentum, velocity, paired | 1.04778046e-4 | 0.339535341 | 106.106561 |
| regular / default | 1.45525064e-4 | 0.244064795 | 147.370223 |
| regular / momentum, velocity, paired | 8.21559538e-6 | 0.339535341 | 9.29223558 |

Grouped policies agree at the displayed precision. Momentum separation
reduces mass/energy mismatch here but increases maximum x-momentum mismatch;
no single flux metric establishes overall accuracy. These unsynchronized
spatial outputs do not constitute a conservation-budget test. The previously
audited completed-trajectory budgets remain the applicable conservation
evidence; no new time-integrated budget is claimed.

This excludes a nonzero initial transverse rate or appreciable initial halo
pressure mismatch as the immediate explanation of the later transverse
maxima. It identifies nonstationary streamwise/density/energy behavior before
time integration. Whether subsequent pressure evolution arises through
predictor state, ancestor-time interpolation, reflux or held-energy closure
requires stage-local evidence. The simultaneous whole-tree initial call does
not reproduce all recursive substep clocks. Next locate the first nonzero
transverse rate/pressure departure during those stages, retaining their
energy correction and synchronization ledgers. No production policy or
threshold was changed and all-field convergence remains open.

Final initial-audit checks passed: pinned-mechanism provenance, project
contract, analyzer Python syntax and `git diff --check`. No solver regression
suite was rerun for this analysis/documentation update.

### Recursive-stage trace setup

`build/eb-shear-stages-20261022/build_trace.py` creates a diagnostic copy of
`eb_reactive_patch_tree_3d_mod.F90`, compiles it into an isolated object/module
directory, and links that object with the existing reference probe and runtime
libraries. The production module remains byte-identical. Build commands and
its original SHA256 are retained in `build-commands.json`; instrumented source
and build log are retained alongside it.

Read-only hooks report sequence, operator/closure stage, node/level, pressure
maximum and its cell, transverse momentum/rate maxima, maximum E-H-KE gap and
energy correction value where available. Hooks cover spatial entry/exit,
closure entry, after gas/scalar reflux, after leaf reset, and closure exit.
They include covered parent cells deliberately, since those feed subsequent
halos. Node maxima therefore are not composite-domain error norms. Entry
has depth=-1 when the optional depth selector has not yet been resolved;
exit reports the resolved target. The sequence is call order, not physical
time, and dt=0 on closure rows means no timestep argument is available.
Energy values are incremental after leaf reset and cumulative at closure
entry/exit; zeros at other rows are placeholders, not measured zero changes.

The bounded comparison uses cut and regular 16-grid inputs, all four policies,
and one 50 ns root step. It will compare all 16 existing summary/inventory/
norm/location CSV prefixes per geometry with saved completed trajectories;
no uninstrumented solver rerun is required. Numerical disagreement must be
resolved before interpreting the trace. The trace is not a substitute for
complete recursive clock identity or the existing conservation audit.
The runner records hashes/commands and sends a completion notification on
success or failure. No new scientific conclusion is claimed at setup.

The isolated build completed and driver PID 1908700 started the bounded trace.
`trace-provenance.json` hashes the original/instrumented module, object,
executable and runner scripts. Completion notification was verified armed;
the first trace row was written successfully. Pinned-mechanism provenance,
project contract, Python syntax and whitespace checks passed. The original
module hash was independently verified unchanged. Results and prefix parity
remain pending until the completion notification.

### Completed recursive-stage trace

Both conditions completed with all four policies (cut 61.19 s; regular
118.91 s). `build/eb-shear-stages-20261022/analyze.py` checked all 384 trace
rows for finite values, complete sequence/node coverage and the expected
24-event pattern per policy. Binary/input/instrumentation hashes match.
All 32 summary/inventory/norm/location CSV prefixes at steps 0 and 1 match
the saved uninstrumented trajectories exactly as parsed string-valued CSV
rows. This is output-prefix parity, not whole-trajectory or external-solver
parity. The maximum reported corrected component-budget residual is
6.058795169e-14, below the unchanged 1e-12 limit. `analysis.json` retains
per-policy event evidence and checks. Logs contain no fatal marker; successful
completion notification dispatch was verified. No simulation was rerun during
this analysis.

Within each policy's 24-event sequence, the first coarse spatial operation
(events 1–4, dt=50 ns) produces only roundoff-scale transverse rates:
up to about 6.47e-8 kg/m2/s2 and momentum about 3.23e-15 kg/m2/s. These
nonzero values differ from the earlier simultaneous whole-tree initial call;
the recursive depth-selected call is a different path. They must not be
reported as the resolved transverse defect.

The first fine spatial call (events 5–6, dt=25 ns) still sees zero pressure
offset and zero transverse state/rates on the fine node. At the next fine
spatial entry (event 7), the predictor already has the following pressure
error while both transverse momenta remain exactly zero:

| Geometry / policy | Fine predictor pressure maximum (Pa) | Node 2 cell | Next y/z momentum rate maxima (kg/m2/s2) |
| --- | ---: | --- | --- |
| cut / default | 3.911186951 | (1,2,11) | 6257.992511 / 6257.867959 |
| cut / momentum, velocity, paired | 3.911186951 | (1,2,11) | 6257.992511 / 6078.531347 |
| regular / default | 3.911186951 | (1,2,1) | 6257.992511 / 6257.867959 |
| regular / momentum, velocity, paired | 0.790890201 | (1,8,1) | 2530.312509 / 1265.423058 |

Event 8 is the exit of that spatial call. The pressure locations in this table
are pressure-maximizing cells; the trace did not export locations of rate or
momentum maxima. All pressure locations adjoin the coarse/fine boundary.
The cut location has volume fraction 0.44; regular locations also adjoin the
lower physical wall. At the first closure entry (event 9), before its reflux,
fine y/z momenta are already 7.82249e-5 / 7.82233e-5 for cut/default and
regular/default, 7.82249e-5 / 7.59816e-5 for cut candidates, and
3.16289e-5 / 1.58178e-5 kg/m2/s for regular candidates.

Thus the resolved pressure departure is bracketed between the first fine
spatial exit and the next fine spatial entry; the resolved transverse rate
appears in the following spatial evaluation, before the first observed
closure/reflux. This narrows the source to preparation of the fine predictor
and its thermodynamic state. It does not isolate every operation within
that interval, including any earlier energy reset. The pressure in event 7
is reconstructed from valid state and held energy, not a halo-only pressure
measurement. The preceding nonstationary density/energy rates are consistent
with this sequence, but are not by themselves a complete causal proof.

The first closure's reflux does not change these fine pressure/transverse
maxima (events 9–10), while it changes the coarse state. For example,
cut/default coarse pressure maximum grows from roundoff to 0.00494272 Pa
and then to 0.00933979 Pa at closure exit, which also includes average-down.
Covered parent cells are included in these node maxima, so these coarse
values cannot be read as composite-domain error norms or as a pure reflux
contribution. In every observed leaf-reset pair (events 10→11, 18→19 and
22→23), the pressure and both transverse-momentum maxima are unchanged.
The final reset does change the energy ledger: -1.74866e-11 J for cut/default,
-7.48588e-12 J for cut candidates, -4.51497e-11 J for regular/default and
-5.03982e-12 J for regular candidates. These increments are not to be summed
once per node: the trace duplicates the same scalar ledger on both node rows.
Earlier cumulative corrections are already nonzero on closure entry, so this
observation does not exclude a role for energy handling before the first hook.

At 50 ns, fine pressure maxima are 7.50774250 Pa for default in both geometries,
7.49927737 Pa for cut candidates, and 1.54444719 Pa for regular candidates.
The first-step candidate distinction is driven by momentum separation;
consistent-velocity and paired-pressure do not distinguish these reported
fine maxima at the displayed precision. A reset that restores E-H-KE does
not erase the observed pressure/transverse error.

The next discriminator is a cellwise density/held-energy/EOS balance across
the first fine predictor interval (events 6→7), including any reset and
ancestor-time boundary preparation inside it. The current trace establishes
the ordering before reflux but does not yet justify changing a production
numerical policy. All-field convergence and original-reference qualification
remain open.

### Fine predictor cell-balance follow-up

The next isolated trace build in `build/eb-shear-predictor-20261022/` adds
cellwise observations at the existing spatial/closure hooks and immediately
before `recover` and after its `reset_held_eb_tree_family_3d` call. Source
inspection shows that the predictor updates gas state, calls `held_update`,
and then calls `recover(valid,1)` before the next source evaluation. The
observation distinguishes that update from the following reset without
changing either operation.

Selected fine cells have i=1, j in {2,8}, k in {1,11}; covered cells are
skipped. These include the previously observed cut and regular predictor
pressure maxima and matched control cells. The new `.csv.cells` output
records rho, held H, total E, all momenta, kinetic energy, held-EOS pressure
and temperature, plus spatial rho/E/H rates where those are available.
Rate zeros outside spatial-exit hooks are placeholders. The event counter
is call order, not physical time; it is independent of the node-summary
counter. Node 2 is intentionally specific to this two-level diagnostic.

The runner again bounds the calculation to one 50 ns root step for each
geometry and four policies, with a 600-second per-geometry timeout and
same-thread completion notification. It checks the 32 saved output prefixes
against the completed uninstrumented runs before declaring success. This is
a new diagnostic measurement of the same short interval, not a duplicate
background job or a full rerun of the earlier 200 ns batch. Production module
source is unchanged. No predictor-balance conclusion is claimed at setup.

The isolated predictor trace build completed. Driver PID 1933200 launched
successfully, completion notification was verified armed and the first
cellwise row was observed. `trace-provenance.json` records executable/object/
source/runner hashes; the production module matches the retained original.
Mechanism provenance, project contract, Python syntax and whitespace checks
passed. `analyze_cells.py` is prepared to compare the first fine spatial rate
with the pre-recover predictor increment and the post-reset state; it has only
been syntax-checked and its numerical analysis awaits completed output.

### Completed fine predictor cell balance

Both instrumented conditions completed (cut 61.45 s; regular 118.73 s),
all four policies, returncode 0. `analyze_cells.py` and `verify.py` in
`build/eb-shear-predictor-20261022/` evaluated 24 selected cell/policy pairs:
two fluid cells per cut policy and four per regular policy. Binary/input/
instrumentation hashes match; all 32 saved CSV prefixes match as parsed
string-valued rows. Logs have no fatal marker and completion notification
was successfully queued. No calculation was launched in this final audit.

The first fine predictor obeys `U_pred = U_old + 25 ns * RHS` for density,
held internal-energy density H and total-energy density E. Maximum absolute
residuals over the sampled cells are 1.81012e-17 kg/m3, 2.12338e-12 J/m3,
and 1.26013e-11 J/m3 respectively, each within two ULPs of the resulting
state. Thus the sampled discrepancy is present in the spatial update, not
an unexplained mismatch between the exported rate and the predictor.

At node 2 (1,2,11), all cut policies have the same first predictor:

| Quantity | Before update | Before recover | After reset |
| --- | ---: | ---: | ---: |
| Density (kg/m3) | 0.379328123157 | 0.379313501414 | 0.379313501414 |
| Held H (J/m3) | 145435.213732 | 145429.592054 | 145429.592054 |
| Total E (J/m3) | 145436.806812 | 145427.273792 | 145431.159735 |
| Pressure (Pa) | 101325.000000 | 101321.088813 | 101321.088813 |
| Temperature (K) | 900.000000000 | 899.999951407 | 899.999951407 |

The spatial density rate is -584.869711618 kg/m3/s; H and E rates are
-2.24867104355e8 and -3.81320788603e8 J/m3/s. The reset adds
3.885943408 J/m3 to E, restoring E=H+KE at this cell. It changes neither
rho, H, momenta, held-EOS pressure nor temperature. Those equalities hold
exactly in the exported before/after data for all 24 sampled pairs; the
post-reset energy identity is within two ULPs. These are local energy-density
changes, not the global signed-energy ledger in joules.

For this fixed-composition ideal gas, the measured pressure change can be
split exactly as
`delta_p = p0*(rho1-rho0)/rho0 + p0*(rho1/rho0)*(T1-T0)/T0`.
This is an algebraic attribution at fixed initial temperature then final
density, not an independent experiment. Its largest residual against the
exported EOS pressure change is 1.723e-11 Pa.

| Geometry / policy / fine cell | Density contribution (Pa) | Temperature contribution (Pa) | Total pressure change before recover (Pa) |
| --- | ---: | ---: | ---: |
| cut / all / (1,2,11) | -3.905716444 | -0.005470508 | -3.911186951 |
| cut / all / (1,8,11) | -0.786462202 | -0.004427999 | -0.790890201 |
| regular / default / (1,2,1) and (1,2,11) | -3.905716444 | -0.005470508 | -3.911186951 |
| regular / candidates / (1,2,1) and (1,2,11) | -0.111936369 | -0.000186867 | -0.112123236 |
| regular / all / (1,8,1) and (1,8,11) | -0.786462202 | -0.004427999 | -0.790890201 |

Here candidates are momentum, velocity and paired. At the cut representative
cell, density change accounts for about 99.86% of the pressure drop; at the
regular candidate maximum sampled cell j=8, about 99.44%. The latter cell
has density rate -117.770434153 kg/m3/s. The smaller j=2 regular candidate
rate is -16.762146721 kg/m3/s. Thus the first pressure defect is dominated
by the local continuity update, with nearly unchanged temperature; it is
already present before `recover`. The observed energy reset does not create
or remove it. This does not exclude earlier halo construction or ancestor
preparation as the source of the nonzero mass flux divergence.

The maximum reported corrected component-budget residual remains
6.058795169e-14, and output parity retains the existing conservation evidence.
Local density error is compatible with global conservation; the next target
is the face-by-face mass flux divergence and reconstructed streamwise
velocity at these coarse/fine-adjacent cells. Neither an EOS patch nor an
energy-reset change is justified by these results. Production numerics and
qualification boundaries remain unchanged.

Verification artifacts: `cell-analysis.json` includes all 24 before/after
states, rate residuals and pressure decomposition; `verification.json`
records hashes/parity/numerical checks. Mechanism provenance, project
contract, Python syntax and whitespace checks were rerun for this update.

### Face balance and connected-direction candidate

`build/eb-shear-face-balance-20261022/analyze.py` reused saved initial face
fluxes and predictor states; it launched no solver. All 24 selected predictor
pairs satisfy the six-face finite-volume density balance to 1e-8 kg/m3/s.
Every y/z contribution is exactly zero. The x right-face mass flux equals
rho times the analytic fine-cell streamwise velocity within 1e-12 kg/m2/s;
the left coarse/fine face supplies too little mass. Aperture-weighted flux
imbalance divided by fluid volume reproduces the predictor density rate and
its pressure contribution. Results and input hashes are in `analysis.json`.

For cut/paired at fine (1,2,11), x-left/right intensive mass fluxes are
0.916591648887 and 1.099363433767 kg/m2/s. Both apertures and cell volume
fraction are 0.44, so that factor cancels; the density rate is
-584.869711618 kg/m3/s. The corresponding ghost/valid velocities are
1.938391787407 and 2.898185941550 m/s. For regular candidates at j=2,
the ghost velocity is 2.870699858852 m/s and density rate is reduced to
-16.762146721 kg/m3/s. At j=8 the ghost remains 9.744953584044 m/s versus
9.935868511442 m/s in the fine cell, leaving -117.770434153 kg/m3/s.
Thus small-cell volume amplification is not the explanation of these selected
cut-cell continuity errors. The mass-flux discrepancy is localized at the
coarse/fine face; no claim is made that its Riemann flux equals a simple
arithmetic velocity average.

Source inspection identifies a relevant distinction: `plan_eb_linear_halo_3d`
returns piecewise constant interpolation whenever any of the 3x3x3 donor
connectivity entries is false. The cut donor adjacent to the planar wall has
covered neighbors in z, so tangential y interpolation is also suppressed.
The saved cut ghosts equal the coarse cell averages for all four policies.
Regular candidates retain y slopes at j=2, while the limiter at the sine
extremum flattens the j=8 donor. The original connectivity rule is documented
as part of the pinned EB interpolater; changing it is a numerical policy
experiment, not a proven porting correction.

An isolated candidate in `build/eb-shear-tangent-candidate-20261022/` permits
linear reconstruction separately along axes whose two immediate neighbors
are connected. A blocked axis uses the donor for both slope neighbors, hence
zero slope. Existing missing-node/inflow checks, limiter arithmetic and EOS
admission are retained. This prototype changes only a copied transfer module;
the production source hash matches `module.before`. It is not an application
option or an adopted policy, and corner/curved-EB admissibility is unqualified.

The candidate built successfully and driver PID 2105981 started one 50 ns
root step for cut/regular and all four policies. Each condition is bounded
by 600 seconds. Cut outcomes will be compared against saved baseline errors;
regular output must match all 16 saved CSV prefixes exactly. The original
1e-12 corrected-budget threshold is enforced, with no golden or tolerance
changes. `provenance.json` and `build-commands.json` retain the candidate
source/object/executable hashes and commands. Its finalization handler queues
one same-thread completion notification on success or failure. Candidate
improvement is not claimed before evaluation.

### Completed connected-direction candidate comparison

The isolated candidate completed cut and regular (54.39 s and 103.64 s)
for all four policies, each to 50 ns. `analyze.py` under
`build/eb-shear-tangent-candidate-20261022/` independently recomputed 224
component-inventory residuals and checked 272 norm/location pairs, binary/
input/provenance hashes, finite values and unchanged regular outputs.
All 16 regular summary/inventory/norm/location CSV prefixes agree exactly
with the saved baseline as parsed string-valued rows. All eight corrected
budgets pass the original 1e-12 limit, maximum 6.058795169e-14. Logs have
no fatal marker and successful completion notification dispatch is recorded.
No simulation was launched during this audit.

At the cut endpoint, momentum/velocity/paired each improve all three norms
for every nonzero-error field (density, all momenta, total E, N2 density,
held H, temperature and pressure). The eight absent-species fields remain
exactly zero. Paired ratios relative to its original connectivity baseline
are below; momentum and velocity ratios agree at displayed precision.

| Field | L1 candidate/baseline | L2 candidate/baseline | Linf candidate/baseline |
| --- | ---: | ---: | ---: |
| density | 0.302413 | 0.273128 | 0.204922 |
| x momentum | 0.307030 | 0.266780 | 0.203919 |
| y momentum | 0.517671 | 0.620277 | 0.427963 |
| z momentum | 0.054401 | 0.103527 | 0.209603 |
| total E | 0.902866 | 0.451987 | 0.210827 |
| N2 density | 0.302413 | 0.273128 | 0.204922 |
| held H | 0.302421 | 0.273229 | 0.206964 |
| temperature | 0.427462 | 0.480048 | 0.389941 |
| pressure | 0.302406 | 0.273176 | 0.205946 |

Paired pressure Linf decreases from 7.499277372 to 1.544447190 Pa, and
pressure L2 to 0.0677830772 Pa. Its y/z momentum Linf becomes
1.265878162e-4 / 6.229499697e-5 kg/m2/s. These candidate cut maxima equal
the regular candidate maxima at this endpoint, but this is not whole-field
parity: geometry, volume and wall location still differ. Total-energy L1
only improves by about 9.7%, despite larger L2/Linf reductions; report all
norms rather than selecting only the strongest improvement.

**Default-policy regressions remain.** With the shared original limiter,
y-momentum L1 increases 4.6702% and L2 increases 5.4677%; x-momentum Linf
increases 0.0009946% (0.017412123607 → 0.017412296794 kg/m2/s).
These are the only increases among all 204 cut field/policy/norm comparisons;
other nonzero default norms decrease or remain equal. The candidate therefore
does not support a blanket replacement of the pinned connectivity guard.
The initial intuition that preserving tangential slopes helps is supported
for the separated-momentum policies in this bounded case, not universally.

`analysis.json` and `all-field-comparison.csv` retain all 17 fields, three
norms, four policies and exact baseline/candidate values. Output-prefix
parity and corrected budgets are distinct from physical accuracy. This
experiment covers one resolution and one root step with a planar wall;
200 ns evolution, refinement convergence, curved/corner connectivity,
conservative reconstruction, positivity/rollback, MPI and restart policy
identity have not been qualified for this new candidate. No application
option, production source, golden reference or tolerance was changed.

The next adoption gate is a bounded opt-in connected-direction design tied
to the separated-momentum policy, with geometry/stencil unit coverage and
longer matched shear/refinement evidence. Default/shared-limiter behavior
must remain pinned. The current result is a promising isolated numerical
candidate, not a completed production fix.

### Candidate refinement batch and opt-in kernel

The next batch uses the frozen isolated candidate at 8/16/32 and the original
four 50 ns root steps to 200 ns, with all four policies. Saved baseline runs
are reused. Driver PID 2155192 owns
`build/eb-shear-tangent-convergence-20261022/record.json`; the numerical runner
owns `runs/record.json`, per-grid logs, executable and input hashes. The
per-grid timeout is 3600 s, and the existing 16-grid runtime forecast can
refuse grid 32 if it exceeds that bound. Completed grids remain retained on
failure/cost refusal. The same-thread finalization notification was verified
armed. No convergence claim is made at launch.

In parallel, the planner's optional `connected_directions` selector was added
with false as default and no application caller enabling it. Focused unit
coverage checks the retained omitted/explicit-false fallback, tangential
donor indices, blocked-normal donor substitution, affine tangential values,
and symmetric-child mean preservation. The relevant target built successfully and `unit_eb_reactive_transfer_3d`
passed 1/1 in 0.07 s, with zero failures or skips. This kernel implementation and the
frozen numerical experiment are distinct artifacts.

`kernel-verification.json` records the dirty source hashes, HEAD, GNU compiler,
build log and exact CTest command. This is a focused kernel test, not a full
matrix or installed/application/MPI/restart qualification. Provenance,
project contract, Python syntax and whitespace checks passed. The copied
independent budget auditor and all-norm baseline/order analyzer are prepared
for batch completion; only syntax checks have run on them at this stage.

### Completed candidate 8/16/32 comparison at 200 ns

All three grids completed with all four policies: 70.83 / 213.67 / 942.04 s,
respectively. No grid was skipped, timed out or restarted. The successful
completion notification was queued. `audit_results.py` independently
recomputed 840 component-inventory residuals and checked 1020 norm rows,
initial/time coverage and frozen executable/input hashes. All 12 corrected
budgets pass the unchanged 1e-12 threshold; maximum residual is
2.070507750e-13. The numerical-energy correction remains part of that budget;
this is not uncorrected total-energy conservation. Its magnitude relative to
initial subcell kinetic variance ranges from 1.000019 to 1.002212. Logs retain
array-temporary warnings and IEEE underflow/denormal flags but no fatal marker.

`compare_baseline.py` produced 612 endpoint comparisons (three grids, four
policies, 17 fields, three norms) against the saved original experiment, and
204 all-norm grid sequences. No solver was launched during this audit.
All eight absent-species fields remain zero in every case; N2 density follows
total density. The paired candidate results are:

| Field | L2 at 8 / 16 / 32 | L2 orders 8→16 / 16→32 | L2 ratio to original paired at 32 |
| --- | --- | --- | ---: |
| density | 3.81516e-06 / 9.04969e-07 / 1.79508e-07 | 2.076 / 2.334 | 0.215780 |
| x momentum | 0.00242613 / 0.000593818 / 0.000127479 | 2.031 / 2.220 | 0.182073 |
| y momentum | 0.000133778 / 7.12531e-05 / 3.24318e-05 | 0.909 / 1.136 | 0.805114 |
| z momentum | 5.0757e-05 / 1.72542e-05 / 5.12039e-06 | 1.557 / 1.753 | 0.039286 |
| total E | 1.56968 / 0.376453 / 0.0779969 | 2.060 / 2.271 | 0.242518 |
| N2 density | 3.81516e-06 / 9.04969e-07 / 1.79508e-07 | 2.076 / 2.334 | 0.215780 |
| held H | 1.46412 / 0.347615 / 0.0690666 | 2.074 / 2.331 | 0.216187 |
| temperature | 4.90764e-05 / 1.32638e-05 / 3.02021e-06 | 1.888 / 2.135 | 0.347974 |
| pressure | 1.01956 / 0.241955 / 0.0480333 | 2.075 / 2.333 | 0.215981 |

All nonzero L1 and L2 errors decrease on both refinements for each of momentum,
velocity and paired. At 32, every nonzero norm improves against the matching
original policy. Paired pressure L2 is 0.0480332742 Pa (78.40% reduction), and
pressure Linf is 2.286939598 Pa versus 18.2824 Pa (87.49% reduction).
The paired z-momentum L2 reduction is 96.07%, whereas the y-momentum L2
reduction is only 19.49%. Total-energy L1 improves only 16.22%; it must not be
represented by its stronger L2 reduction. Momentum and velocity results are
close to paired; the three selectors do not constitute independent remedies.

**Local convergence remains incomplete.** Paired y-momentum Linf is
0.00195867206 / 0.00203519345 / 0.00184204799 kg/m2/s at 8/16/32, increasing
about 3.91% on the first refinement. Its observed orders are -0.055 / 0.144.
The z-momentum Linf decreases, but slowly: 0.000942255883 / 0.000927056870 /
0.000792909436, orders 0.023 / 0.226. Pressure Linf orders are 1.042 / 1.314.
These are observed three-grid endpoint orders at fixed dt, not proof of an
asymptotic spatial order; no new candidate timestep study was performed.
The only regression versus each original separated policy among all fields
and norms is y-momentum Linf on grid 8, about +0.109% (momentum, velocity,
paired). Therefore the 50 ns finding of universal improvement does not extend
unchanged to every grid at 200 ns.

The original shared/default limiter remains unsuitable for blanket adoption:
y-momentum L2 increases 1.4994% at 16 and 0.5939% at 32 relative to its
original policy, with L1 increases 0.3885% / 0.3304%. Its transverse Linf still
grows strongly with refinement. The exact comparison also records very small
positive differences, including z-momentum Linf +0.0004268% at 16 and
grid-8 relative increases no larger than 1.70e-8. These are retained without
loosening any threshold or claiming bitwise equality.

Artifacts are `audit.json`, `all-field-results.md`, `baseline-comparison.csv`,
`all-norm-orders.json` and `final.json` under
`build/eb-shear-tangent-convergence-20261022/`. The numerical experiment uses
the earlier frozen isolated candidate, not a newly enabled application path.
The separately added optional planner selector has a passing focused unit
test, but applications still do not enable it. These results support further
bounded opt-in development with momentum separation; they do not close
all-field/local convergence, curved/corner EB validity, MPI/restart lifecycle
or original-PeleC/physical qualification. No production default was changed.

### Connected-direction geometry admission hardening

The optional production planner kernel now requires connectivity of the full
retained subspace, not just the two axial neighbors. A covered donor center,
no usable direction, or disconnected diagonal within that subspace returns
PCM. Blocked directions still request the center twice; the default remains
false and no application caller enables this option. The planar-wall tangent
stencil is retained, while the earlier frozen numerical candidate lacks this
new diagonal guard and remains separately identified.

The extended `unit_eb_reactive_transfer_3d` passed 1/1 (0.07 s, no final failures
or skips). It retains the default/explicit-false tests and planar tangential
interpolation/mean checks, adds a synthetic sealed-diagonal topology that
must fall back, and enumerates 64 coarse donor positions in each of two
nodal level-set geometries. The sphere admits 24 and falls back at 40;
the corner admits 36 and falls back at 28. All admitted requests stay on the
available coarse level and refer to positive-volume fluid cells. The synthetic
aperture case tests connectivity logic only, not consistency of a physical EB
metric construction. Enumeration is not a proof for every curved geometry.

An initial fixture check failed because the corner at x=.13,z=.17 produced
no PCM donor cases in the sampled physical cells. Moving its walls to
x=.38,z=.42 exercises both branches; the test assertion and production
admission guard were retained. `build/eb-connected-admission-20261022/`
contains the final test log and source hashes/HEAD/command in `record.json`.
The application/MPI/restart selection route remains unimplemented; curved
flow accuracy, volume-weighted reconstruction conservation and local shear
convergence remain open. No long solver or duplicate batch was launched.
