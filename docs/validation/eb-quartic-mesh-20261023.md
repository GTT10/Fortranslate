# Quartic candidate mesh diagnostic — 2026-10-23

This study compares the same frozen roundoff-scale candidate at 8/16/32 root
cells and **50 ns**. It is not the earlier 200-ns study, fixed-PeleC trajectory
parity or physical qualification. At this record's creation, 8 and 16 are
complete; 32 is running with a 900-second limit.

## Matched controls and baseline selection

Artifacts: `build/eb-quartic-mesh-20261023/`. `run_grid.py` runs one requested
grid, refuses to overwrite any existing per-grid receipt, checks frozen binary
and input hashes, preserves the 1e-12 budget gate, checks complete output
tables and compares default output with the original baseline. The existing
16-cell result is reused from `eb-quartic-roundoff-20261023`; it is not rerun.
`summarize.py` checks binary equality before deriving all L1/L2/Linf ratios
and observed refinement slopes at the same physical time.

The historical connected-direction candidate also altered its default branch,
so its default results are **not** this study's default reference. Default
uses `build/eb-shear-analytic-20261022/runs/`; the separated policies use
`build/eb-shear-tangent-convergence-20261022/runs/`. Their first-step 16-cell
prefixes reproduce 15 available current baseline tables exactly. Historical
original-default location tables are missing, so no location-parity claim is
made for that branch. Input files differ from the earlier 200-ns inputs only
in the requested step count (one rather than four); 16-cell verification
found no other line difference beyond the equivalent numeric time spelling.

## Completed 8-cell pilot and interpretation

The 8-cell run completes in **30.156 seconds**. Default ordinary outputs agree
exactly with the original first-step baseline. All four budgets pass; the
maximum normalized residual is 1.05e-14. There are **zero quartic adoptions**,
and the separated-policy norms equal the connected baseline. Consequently,
the large 8-to-16 drop includes a change from complete fallback to partial
quartic admission. It must not be reported as formal high-order convergence.

An independent analytic-initial-stencil calculation,
`analyze_initial_screen.py`, reproduces the expected screening transition:
no coarse y positions pass at 8, 4 of 16 pass at 16, and 20 of 32 pass at 32.
It uses exact periodic shear averages for momentum and kinetic energy and
the unchanged curvature-sign rule. This calculation excludes geometry/EOS
admission and time evolution; its counts are not measured halo adoption counts.

At 8/16 the paired y-momentum Linf values are 1.22608e-4/8.52889e-6, and
z-momentum values are 6.09250e-5/8.72278e-6. The existing baseline maxima do
not decrease over that pair. These are useful error measurements, but the
32-cell result is needed before deciding whether to extend this candidate.

The 32-cell wall-time forecast is about 636 seconds, derived from the measured
144-second candidate run at 16 and the prior 16-to-32 timing ratio. This is a
forecast, not observed runtime. Results, failures and partial outputs remain
separate in the per-grid receipts. No production selector or checkpoint policy
has been changed.

## Partial 32-cell output audit

While the full sweep is still running, default and momentum settings have
emitted all four complete tables. `grid32-partial-audit.json` records their
hashes and checks original-default native, inventory and analytic equality.
This receipt explicitly marks the full run incomplete; velocity and paired
are not admitted from this partial audit.

For momentum separation, 16-to-32 y/z Linf decreases from
8.52889e-6/8.72278e-6 to 4.26088e-6/4.35778e-6. The corresponding observed
slopes are 1.00120 and 1.00120. Pressure Linf decreases from 0.214495 to
0.0527575 Pa (slope 2.02350). These first-order transverse maximum-error
reductions improve on the previous near-flat maxima, but do not establish
formal asymptotic order, all-policy convergence or longer-time behavior.


## Terminal sweep and bounded recovery

The original 32-cell sweep hit its 900-second bound. Its failed receipt is
preserved; no whole-sweep success is claimed. Terminal CSV validation in
`grid32-terminal-audit.json` admits default, momentum and velocity: each has
2 native, 28 inventory, 34 analytic and 34 location rows. Default matches its
original reference tables exactly. Momentum and velocity both give y/z Linf
4.26088e-6/4.35778e-6 and pressure Linf 0.0527575 Pa; their maximum normalized
budget residual is 8.48e-14. Paired emitted initial rows only.

`recover_paired.py` runs only the missing paired policy with a 600-second
bound, using the previously verified paired-only executable and unchanged
32-cell input. Its separate receipt is `grid32-recovery.json`. The combined
`audit_recovery.py` includes paired only after that receipt reports success;
otherwise its `all_policies_complete` remains false. The recovery completed in 356.388 seconds with return code zero. Production
interpolation and qualification are unchanged.


## Recovered complete 32-cell comparison

`grid32-terminal-audit.json` now validates all four policies and all 16 tables;
`all_policies_complete=true`. This combines three completed policies from the
interrupted sweep with the separately verified paired-only recovery. The original
sweep remains recorded as failed. All 17 components and L1/L2/Linf norms were
compared to the appropriate connected baseline; no candidate norm increased.
The default tables remain equal to the original default reference.

Momentum, velocity and paired give identical transverse and pressure maxima:
y-momentum 4.26088431777e-6, z-momentum 4.35777518824e-6, pressure
0.05275745498 Pa. Relative to their connected baseline these are 3.33376%,
7.05502% and 6.94838%, respectively. Their maximum normalized budget residual
is 8.47783e-14. The 16-to-32 slopes are 1.00120, 1.00120 and 2.02350.

At 32 cells, the y maximum is in fine cell (1,1,22), volume fraction 0.88;
the z maximum is in full-volume fine cell (1,1,32); the pressure maximum is
in full-volume fine cell (1,1,23). All are coarse/fine adjacent. Thus improved
amplitude does not remove the interface-localized residual or prove a general
high-order method. The next numerical question is the stage-level forcing at
these remaining maxima, including their donor admission and face pressure.
This is still a 50-ns isolated serial candidate, not fixed-PeleC trajectory
agreement, MPI/restart qualification, or physical validation.

Recovery tooling passes Python syntax checks; mechanism provenance, project
contract and `git diff --check` pass. No production numerical change was made
in this recovery.


## Remaining-extrema diagnostic extension

The reference probe now adds fine-cell IDs 6:9 at (1,1,zcut), (1,1,N),
(1,1,zcut+1), and (1,1,N-1), preserving historical IDs 1:5. Both stage
pressure/face flux and synchronization diagnostics use those IDs. This does
not change the numerical method. The isolated paired executable, exact build
commands, source hash and frozen candidate linkage are recorded in
`build/eb-quartic-extrema-20261023/build-record.json`.

The 16-cell diagnostic completed in 50.066 seconds. All four ordinary tables
are byte-identical to the saved candidate; all historical diagnostic rows
also match exactly. The expanded trace has 144 spatial rows and 324 sync
rows. `audit.py --grid 16` validates directional face-rate closure, 27
pre/post synchronization pairs, and all six fine-probe momentum trajectories.
At y-maximum probe 9, integrated spatial change is 8.52888765187e-6 and
applied closure is zero. At z-maximum probe 7, the corresponding values are
8.72277621234e-6 and zero. Cut-adjacent probes 6/8 instead receive a small
y closure of -2.95266107364e-8. Thus closure affects some nearby cells, but
not these 16-cell global transverse maxima. The same bounded diagnostic is
running at 32 cells; no 32-cell stage result is admitted yet.


The validated 16-cell forcing summary (`grid16-forcing-summary.json`) narrows
the remaining mechanism further. At fine probes 7 and 9, the initial y/z
rates vanish. At the first predictor evaluation, cell pressure has risen by
0.111649534 Pa. The y directional rate is 179.397254 (momentum-density per
second), entirely from y faces; probe 7 also has z rate 178.639279, entirely
from z faces. Raw and centroid fluxes agree at these probes. The first
predictor y-normal flux difference is approximately 0.0560616 Pa and the
z-normal difference at probe 7 is approximately 0.0558248 Pa, divided by
fine spacing 0.0003125 m. These measurements identify transverse forcing
following the predictor pressure disturbance, not a reflux or centroid
correction origin at these points. They do not yet identify which donor
state or interpolation operation creates that pressure disturbance.
`summarize_forcing.py --grid 32` will apply the same analysis only after
the 32-cell trajectory audit passes.


An analytic initial-state check (`analyze_boundary_shear.py`, saved as
`analytic-boundary-shear.json`) evaluates the actual MC momentum formula at
parent y=1, child offset -1/4, for constant-density cell-averaged sine shear.
Here MC retains the centered slope; this is not an extremum with a zero slope.
The predicted halo velocity exceeds the exact fine average by 0.0274861 m/s
at 16 and 0.00351997 m/s at 32 (observed analytic slope 2.96507).
The earlier analytic curvature screen rejects parent y=1 at both grids.
Thus a possible residual path is fallback interpolation near the shear
inflection, followed by pressure generation and a transverse pressure gradient.
Two spatial derivatives of an O(h^3) velocity mismatch could explain the
O(h) transverse error at fixed short time. This is a hypothesis, not a measured
runtime donor attribution: evolved/ancestor donors, geometry admission and
Riemann coupling are excluded from this analytic check. Do not loosen the
curvature predicate or promote the candidate on this evidence alone.


## Completed 32-cell stage attribution

The expanded paired diagnostic completed in 351.130 seconds with return code
zero. All four ordinary tables remain byte-identical to the recovered 32-cell
candidate. `audit.py --grid 32` passes finite/face completeness, directional
rate closure, 27 sync pairs and six fine-probe trajectories. The maximum
trajectory reconstruction residual is 3.50e-16 (16-cell: 9.31e-17).

At probes 7/9 the first predictor pressure offset decreases from 0.111649534
to 0.0285644640 Pa, a two-grid slope of 1.96668. The predictor y directional
rate decreases from 179.397254 to 91.6058918; probe 7 z rate decreases from
178.639279 to 91.4062880. This supports the measured pressure-disturbance to
transverse-gradient path, but does not prove its donor/interpolation origin.

The 32-cell global y maximum is probe 6, unlike the 16-cell probe 9 maximum.
At probe 6 its spatial contribution is 4.28195429173e-6, with applied closure
-2.10699739594e-8, giving 4.26088431777e-6. The 32-cell z maximum remains
probe 7: spatial contribution 4.35777518824e-6, closure zero. Therefore a
blanket statement that all maxima have zero closure would be incorrect.
The forcing summaries include wall terms explicitly: at cut probes the large
z pressure-face force cancels the wall pressure force; face-only sums must
not be interpreted as net acceleration.

Both grids and their audit hashes are saved under
`build/eb-quartic-extrema-20261023/`. No diagnostic job remains active.
Next evidence needed is runtime donor state and candidate admission at the
y=1 coarse/fine boundary before modifying the curvature predicate. C2 stays
open; no production interpolation change is included.
