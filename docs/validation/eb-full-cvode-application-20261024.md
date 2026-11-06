# Full-coordinate CVODE EB-tree application lifecycle, 2026-10-24

This is focused application lifecycle evidence, not C2 completion or physical
validation. The dirty `main` source identity is recorded in
`build/eb-native-atol-trajectory32-20261024/full-policy-integration.json`
(base `0a0e00ebe9138145fbb0e0e745665d0d768df009`, per-file SHA256).

The explicit `cvode_full` selection connects the full species-density and
temperature reactor through the shared serial/MPI EB-tree front end. Native
and reduced CVODE dispatch retain their existing defaults and behavior.
Checkpoint controls retain native/reduced encodings and reserve value 4 for
full-coordinate CVODE, independently of the chemistry-enabled bit.

## Completed checks

The existing GNU Debug/OpenMPI/SUNDIALS 7.2 build at
`build/eb-eos-application-20261021/isolated-openmpi/registration` rebuilt the
selected serial/MPI applications and full-H2/O2 fixture applications.
No compiler, MPI or optional-backend configuration was changed.

Five CTests passed, zero failed, in 93.66 s:

- `regression_selected_eb_tree_hot_cvode_restart` (39.61 s)
- `regression_selected_eb_tree_full_cvode_restart` (52.56 s)
- `unit_sundials_signed_trials` (1.16 s)
- `unit_sundials_forced_reactor_full` (0.15 s)
- `unit_sundials_forced_reactor` (0.13 s)

The full-coordinate application gate checks chemistry intervention, inventory
and EOS consistency, identical serial/two-rank/changed-rank final CSVs,
intermediate snapshots, and rejection of native/reduced reactor identities
when reading a full-coordinate checkpoint in serial or three ranks.
An additional saved-checkpoint check passed MPI-to-serial continuation with
byte-identical final CSV and reduced-to-full refusal in serial and three ranks.
Logs and exact commands are in `full-policy-tests.log` and
`full-interchange/result.json` under the integration artifact directory.

Lidryer and full-H2/O2 mechanism provenance, Python syntax, project contract
and `git diff --check` passed. Existing executable-stack linker warnings
remain; they did not prevent these builds or checks.

## Fidelity boundary and next work

The earlier 32-cell/100-ns fixed-reference pass used isolated diagnostic
objects. Replaying three captured inputs with the new source reactor and
those diagnostic dependencies reproduces their outputs exactly. Replacing
the diagnostic reactive/thermo/kinetics objects with production libraries
introduces small nonzero output differences; the native-activation mechanism
object is still retained in that replay. Neither is an unmodified application
trajectory result.

The replay explicitly supplies density-space absolute tolerance 1e-15 through
its adapter. The application caller normally multiplies the configured
mass-fraction tolerance by local density. This policy and the diagnostic
thermo/kinetics/native-activation arithmetic must be reconciled before
transferring the 32-cell result to a source-built trajectory. Mesh/time
convergence, broader fixed-reference agreement, C3-C5 and physical validation
remain unqualified.

## Extended permanent regression

The permanent full-coordinate restart CTest was extended to assert the
checkpoint-specific refusal message, reduced-to-full mismatch rejection in
serial/three ranks, unchanged checkpoint contents, and MPI-to-serial final
CSV identity. The extended test passed 1/1 in 67.52 s (CTest wall 67.56 s).
This is a later rerun of the same test, not an additional independent test
in the earlier five-test total. The exact source SHA and result path are in
`build/eb-native-atol-trajectory32-20261024/production-handoff.json`;
`full-policy-extended-tests.log` holds CTest output.

## Production-library trajectory and timing

The canonical-mechanism, production-library Debug 32-cell/100-ns probe reached
its wall-time limit at 1200.061 s. Only the initial state and budget were
written; no final-state comparison was possible. The probe has no checkpoint
restart path. The run and recovery records are in
`build/eb-native-atol-trajectory32-20261024/production_trajectory/`.
Its 775872 array-temporary warnings do not by themselves establish the cause
of the runtime cost.

The unchanged binary completed the bounded 8-cell/100-ns input in 45.357 s.
The original fixed-reference comparison reported `diagnostic_passed` with
zero violations. Exact input/binary identities and comparison evidence are in
`build/eb-native-atol-trajectory32-20261024/production_performance8/run.json`
and `diagnostic-comparison.json` in that directory. This is a source-built
reference probe using production libraries and canonical LiDryer, not an
application CLI trajectory, a 32-cell result, or convergence qualification.

The clean Release build with explicit system OpenMPI libraries completed.
All five selected lifecycle/unit CTests passed in 37.57 s, including both
full-coordinate and reduced-coordinate restart regressions. Exact commands
and timings are in `build/eb-full-cvode-release-isolated-20261024/run.json`;
`tests.log` in that directory contains the five results. The same-input
Release 8-cell comparison completed in 14.519 s (3.12 times faster than
Debug), with `diagnostic_passed` and zero violations. Its identities and
comparison are in `build/eb-native-atol-trajectory32-20261024/release_performance8/`.
The 32-cell/100-ns run reused that exact binary and completed in 712.416 s
within its 1800-s wall limit. The original comparison reported
`diagnostic_passed` with zero violations. Maximum absolute differences were
1.869e-7 and 1.711e-7 kg m^-2 s^-1 for y/z momentum and 4.923e-5 K for
temperature. Exact binary/input identities and per-field errors are in
`release_trajectory32/run.json` and `diagnostic-comparison.json` under the
same artifact root. This establishes the bounded production-library probe
comparison, not application CLI trajectory agreement or mesh/time convergence.
No acceptance tolerances or golden data were changed for this investigation.

## Application initial-state admission

The optional `initial_state_file` path now imports complete node states,
including covered coarse donors, into nonrotating gas applications with held
internal energy. MPI reads on rank zero and broadcasts the input; owner-local
candidate arrays are committed after admission. The density-sheet default
is retained. The format is documented in the EB-tree case README.

The permanent `regression_selected_eb_tree_initial_state` passed 1/1 in
24.41 s (23 executable invocations, not 23 independent CTests). It checks
single-node seed reproduction, multilevel serial/two-/three-rank byte-identical
initial and one-step CSVs, species ordering, missing/duplicate/trailing rows,
mass/temperature mismatch, nonfinite density and restart conflict refusal.
The two existing hot/full CVODE restart regressions passed 2/2 in 38.08 s.
Records, source/binary hashes and logs are linked from
`build/eb-native-atol-trajectory32-20261024/initial-state-handoff.json`.

Canonical LiDryer serial/MPI application fixtures build, but the standard
LiDryer bundle lacks the pinned Simple transport fits used by the reference
probe. Those fits must be imported with provenance before application
trajectory comparison. Neither initial-state admission nor these lifecycle
checks close C2 trajectory/convergence or physical qualification.

The canonical Simple-transport preparation is now available as
`tools/prepare_lidryer_simple_transport.py`. It checks canonical YAML bytes
against the pinned TGReact source and imports the header with SHA256, Git
blob and species-order validation. The bundle regeneration whitelist now
accepts that existing pinned manifest; no chemistry parameters are changed.
The command used was:

```sh
python3 tools/prepare_lidryer_simple_transport.py --reference-source build/pelec-pinned-hit-source --output build/eb-native-atol-trajectory32-20261024/lidryer_simple_bundle
```

Preparation and provenance checks passed. The artifact's `checks.json`
records exact source/bundle hashes, unchanged canonical chemistry, coefficient
round-trip identity, and rejection of modified header bytes, species order
and commit. The existing provenance Python suite passed 14/14; Python syntax,
project contract and `git diff --check` passed. These checks were run on dirty
main at `0a0e00ebe9138145fbb0e0e745665d0d768df009`. Building the application
against this prepared bundle and comparing its trajectory remain pending.

The separate Release serial/MPI build subsequently completed. The first
ordinary serial 8-cell application stopped during initial admission with
`Initial held energy mismatch` (0.266 s), before trajectory comparison.
All 32 discrepant fluid cells are covered coarse donors. Their conservative
and held fields reproduce fine-volume restriction to 5.83e-11 in absolute
SI components. The total-minus-held-minus-resolved-kinetic difference equals
the unresolved kinetic variance to 6.85e-11 J/m3. Evidence is in
`application8/covered-energy-audit.json` under the artifact root above.

Admission now permits this nonnegative covered-cell variance while retaining
the original uncovered-cell consistency tolerance and held-energy EOS test.
It rejects a covered energy deficit. No input values are reset and no
reference acceptance limits changed. Focused serial/MPI tests now include
opposite fine velocities and covered/leaf rejection cases. The incremental
build and expanded test are pending in `covered_admission/run.json`; they
are not yet reported as passed. Application trajectory comparison remains open.

The expanded energy-admission suite subsequently passed all 30 executable
invocations. The next actual-reference attempt passed energy admission but
stopped at `Initial temperature mismatch`. All 32 discrepant temperatures
are covered coarse donors, equal to fine-volume-averaged temperature within
2.28e-13 K (`application8_covered/temperature-audit.json`). Such independent
restriction does not commute with nonlinear EOS inversion. Covered donor
admission now preserves a positive temperature seed and still requires valid
held-energy EOS inversion; uncovered cells retain the original temperature
agreement test. The further expanded suite and rerun remain pending in
`covered_temperature_admission/` and `application8_restricted/` respectively.

The temperature-admission suite passed all 32 executable invocations in
32.27 s. The ordinary serial application then completed the pinned 8-cell,
100-ns input in 19.177 s. Initial and final composite-field gates passed for
544 fluid leaves; initial differences were zero. Maximum final temperature
and y-momentum differences were 8.064e-5 K and 1.435e-9 kg m^-2 s^-1.
The final-time sidecar confirms 100 ns. Correction-accounted energy residual
was -9.54e-18 relative; H2 consumption relative error was 6.575e-6 (limit
1e-3), and species mass closure was 3.34e-16. These selected checks are in
`application8_restricted/comparison_0.json`, `comparison_1.json` and
`budget-checks.json`. They do not certify omitted covered-cell/history output,
mesh convergence, or physical validation. MPI two-/three-rank comparisons also completed under `application8_mpi/`: both initial/final field gates passed, and both CSVs are byte-identical to serial at each rank count.

The expanded inventory checker `check_application_budget.py` in the same
artifact root passed for serial and MPI two/three ranks at 8 cells. It checks
output clock, nonnegative species, species mass closure, mass and H/O/N
molar inventory drift, correction-accounted energy, H2 consumption, and
independent agreement of CSV integrals with the application's component
ledger. Each run's `inventory-checks.json` preserves input/checker hashes.
The 32-cell ordinary serial application is now running with an 1800-s wall
limit (`application32/run.json`). Its initial comparison already passes for
32256 composite fluid cells. Final trajectory/budget results are pending;
this initial agreement does not establish time evolution or convergence.

The additional nonuniform 8-cell restart audit has not passed: serial 200 ns
completed in 32.415 s, but MPI2 checkpoint save at 100 ns was refused before
MPI3 continuation. A diagnostic-only MPI writer change reports covered-root
conservative restriction mismatch 1.204559657350755e-6. Independent analysis
of the saved same-input production-library probe reproduces that exact
maximum in rhoE only. Applying the configured held-energy closure E=H+KE
reduces rhoE mismatch to 2.446e-16; held energy and all other conserved fields
already match to roundoff. Records are in `checkpoint_diagnostic/`, including
`restriction-closure-audit.json`. This identifies a mismatch between the
checkpoint restriction invariant and post-reset coarse donors; no tolerance
was changed. Consistent serial/MPI read/write validation and the actual
restart rerun remain pending. This is separate from the successful MPI
initial/final field comparison without checkpointing.

The candidate checkpoint repair now validates each parent against its saved
immediate children, rather than cascading reconstructed grandchildren through
an already reset parent. Serial and MPI writers/readers share a comparison
that first accepts raw restriction, then permits the exact saved dual-energy
policy on a temporary copy of covered fluid cells. The original 2e-12 scaled
comparison tolerance is retained; checkpoint data are not modified by this
validation. Focused unit coverage has been added for raw and policy-reset
states, disabled reset, corrupted momentum/held energy, uncovered energy,
and a multilevel fixture independently checked against volume-weighted saved
child energy. The expanded unit test subsequently passed (see below). The first
application build hit a Fortran line-length error; the signature is corrected,
and `checkpoint-validation-supervisor.json` records the queued rebuild,
unit CTest, and fresh `application8_restart_policy` audit. This is candidate
implementation evidence only; the restart gate remains open until those
results are inspected.

The expanded `unit_eb_tree_checkpoint_3d` passed independently (1/1, 0.36 s
CTest wall time), recorded in `checkpoint-unit-corrected-result.json` with
test-source and executable hashes. Its first execution exposed a fixture
error: the cold density-sheet state had held energy -105146.92785536687 J/m3,
which exercises the negative-energy branch rather than the intended positive
eta2 reset branch. Explicit positive held energy fixes that test setup;
production implementation and comparison tolerances were unchanged. The test
covers the direct-child multilevel distinction and policy/corruption cases
above, plus its existing checkpoint checks. Actual nonuniform MPI restart
validation remains pending in the separate supervisor.

The ordinary serial 32-cell run completed 100 ns in 981.649 s, below its
1800-s limit, using the launch-time binary hash recorded in
`application32/run.json` (not the subsequently rebuilt checkpoint binary).
Both initial and final field gates pass for 32256 composite fluid cells.
Maximum temperature difference is 4.92347e-5 K; x/y/z momentum differences
are 3.46078e-7, 1.86897e-7, and 1.71114e-7 kg m^-2 s^-1, respectively
(limit 1e-6). Density and total-energy scaled differences are 1.622e-9 and
2.628e-9 (limit 1e-8). `application32/inventory-checks.json` also passes
clock, species closure/nonnegativity, mass/element inventory, corrected
energy and CSV-to-ledger checks. The normalized corrected-energy residual
is 2.670e-17; H2-consumption relative difference is 6.574e-6 (limit 1e-3).
This adds ordinary-application fixed-reference evidence at 32 cells for one
100-ns step; it does not establish multi-step trajectory agreement, formal
mesh convergence, MPI32 parity, or physical validation.

The repaired nonuniform 8-cell checkpoint audit subsequently passed:
`application8_restart_policy/run.json` records serial continuation to 200 ns
(32.365 s), MPI2 checkpoint at 100 ns (19.979 s), and MPI3 continuation to
200 ns (14.059 s). Prefix and final composite CSVs are byte-identical to the
uninterrupted serial outputs, and the checkpoint SHA is unchanged by reading.
`checkpoint-validation-supervisor.json` records the successful application
rebuild and rerun of the expanded unit test. This establishes the tested
MPI2-to-MPI3 restart equivalence; it is not a 200-ns fixed-PeleC comparison.
Serial-save/MPI-read and MPI-save/serial-read checks are separately underway
in `application8_checkpoint_cross/`.

The cross-implementation checkpoint audit also passed in
`application8_checkpoint_cross/run.json`: serial save at 100 ns (19.320 s),
serial continuation of the MPI2 checkpoint (13.252 s), and MPI3 continuation
of the serial checkpoint (14.110 s). All three CSVs match the corresponding
uninterrupted serial outputs byte for byte, and both checkpoint files remain
unchanged. A separate ordinary serial 8-cell, six-step 300-ns fixed-reference
comparison is now running against the saved original 50-ns-step reference;
its per-step and inventory results are pending in `application8_300ns/`.

The ordinary 8-cell, 50-ns-step trajectory now passes all seven output
comparisons from 0 to 300 ns against the saved original fixed reference
(`application8_300ns/run.json`, elapsed 56.596 s). At 300 ns maximum
T difference is 3.65327e-5 K, y momentum difference 1.49586e-8 kg m^-2 s^-1,
and scaled total-energy difference 4.702e-10. Final inventory gates pass:
corrected-energy residual 5.796e-17, species closure 2.221e-16, and relative
H2-consumption difference 2.172e-5. This transfers bounded multistep reference
agreement to the ordinary application; it does not prove time/mesh
convergence. A corresponding 25-ns-step, 12-step trajectory is underway in
`application8_300ns_dt25/` against its saved original reference.

All three ordinary 8-cell temporal trajectories now pass their fixed-reference
field and final inventory gates through 300 ns: dt=100 ns (4 output times,
47.663 s), dt=50 ns (7 times, 56.596 s), and dt=25 ns (13 times, 75.855 s).
The 25-ns final maximum temperature difference is 1.80590e-5 K. Records are
`application8_300ns_dt100/`, `application8_300ns/`, and
`application8_300ns_dt25/`; this is bounded trajectory agreement, not a
blanket convergence verdict.

`analyze_application_time_refinement.py` independently computes composite
fluid-volume L1/L2 and Linf differences between adjacent step sizes at the
same 300-ns endpoint, for both PeleF and original PeleC. It verifies matching
reference contracts apart from dt/step count, identical application initial
states, cell keys/volumes, final clocks, and finite values, and hashes the
input tables. `application-time-refinement.json` preserves the norms and
maximum locations. PeleF L2 observed orders are density 1.12535,
x/y momentum 1.18886/1.18889, z momentum 1.21603, and temperature 1.20696;
the corresponding original PeleC values are 1.12535, 1.18886/1.18889,
1.21603 and 1.20694. These three step sizes support similar roughly
first-order behavior in this bounded case, not asymptotic second order.

Species prevent an all-field monotone-convergence claim: H2 Linf differences
increase from 1.09671e-8 to 1.69489e-8 in mass fraction, and N2 from
7.06468e-9 to 1.20909e-8 (orders -0.6280 and -0.7752). Original PeleC shows
nearly the same differences/orders and the same maximum cells and signs.
All H2/O2/N2 sampled maxima are full-volume fine cells at k=8, the upper
coarse/fine interface; this is not a cut-cell-only effect. The matched
behavior does not identify its cause or close spatial/all-field convergence.
No numerical policy, reference, or acceptance tolerance was changed for this
analysis. All jobs in this bounded application/checkpoint sequence are done.

A saved-output audit now follows species temporal differences at common
100/200/300-ns clocks (`analyze_species_temporal_origin.py`,
`species-temporal-origin.json`). It decomposes each maximum exactly as
DeltaY = Delta(rhoY)/rho_coarser - Y_finer*Delta(rho)/rho_coarser;
all 36 decompositions close within 4.61e-17. The 18 solver-paired maxima
(three clocks, two step pairs, three species) share the same cell locations
in PeleF and PeleC; all are full-volume fine k=8 cells. H2/N2 maximum
step-pair differences decrease with refinement at 100 ns but invert by
200 ns. At 300 ns the H2 100-to-50-ns maximum combines -2.190e-7 and
+2.080e-7 contributions, giving -1.097e-8; its 50-to-25-ns maximum is at a
different cell and has contributions of the same sign. N2 similarly has
strong cancellation in the coarser pair. This identifies a changing balance
and moving maximum, not a proven interpolation, reaction, or flux-sync cause.
The upper fine layer contributes only 6.7%/16.9% of H2 volume-weighted L1
difference in the two endpoint step pairs, illustrating why a local Linf
failure coexists with decreasing integrated norms. A stage-level attribution
must use matched cells/clocks and flux/source data before changing numerics.

The ordinary serial application now exposes opt-in selected-cell MOL stage
and sync observers, documented in the case README. Both selected serial/MPI
applications rebuilt successfully; MPI trace requests are explicitly refused.
The new checker `tests/application/check_eb_tree_mol_cell_trace.py` passed
six application invocations, including diagnostic-off/on byte identity and
invalid count/node/cell and MPI refusal. Its 140 six-face
hydro balance groups close within 0.000e+00
scaled error (840 stage rows, 1512 sync rows). The
checker is registered as `regression_selected_eb_tree_mol_cell_trace`;
this run invoked it directly, not through CTest. Actual reacting dt100/300-ns
trace identity is pending under `application8_300ns_dt100_trace/`.

The ordinary dt100/300-ns reacting trace completed in 47.815 s with all six
ordinary CSV/budget files byte-identical to the saved uninstrumented
trajectory (`application8_300ns_dt100_trace/run.json`). Stage/sync outputs
occupy 2.68/1.75 MB. `audit_application_cell_trace.py` verifies 1512 six-face
hydro groups and 756 closure pairs: maximum normalized hydro-balance
residual is 8.34e-17, raw-flux closure residual zero, and nonenergy/energy
state-versus-applied-increment residuals are below 5.60e-17. The audit checks
all conserved components at the six selected full-volume fine-interface
cells. It does not reconstruct chemistry or the complete SDC update. Matched
50/25-ns diagnostic trajectories at those same cells are now running;
identity and local-balance results are pending.

Both remaining trace trajectories passed ordinary-output byte identity.
The dt50 trace passes 3024 hydro groups and 1512 closure pairs; dt25 passes
6048 groups and 3024 pairs. Maximum normalized hydro residual is below
9.68e-17 and state/applied closure residual below 5.95e-17. Each trace work
directory contains `run.json` and `trace-audit.json`.

`compare_application_final_sync.py --include-dt25` compares both temporal
pairs at the same six cells and 100/200/300-ns endpoints. It requires unique
root pre/post events and checks post-closure density/species densities against
accepted ordinary snapshots, avoiding ambiguous repeated attempts.
`application-final-sync-comparison.json` records 108 species comparisons.
At 300 ns, dt50-to-dt25 H2 maximum at node2 (3,3,8) is -1.69489e-8 both
before and after the last root closure; the N2 maximum at (3,4,8) similarly
remains +1.20909e-8. Both individual closure changes are zero. In contrast,
the dt100-to-dt50 H2 maximum at (4,1,8) changes from -1.23600e-8 before
closure to -1.09671e-8 afterward (difference +1.39285e-9), reducing its
magnitude. N2 there also decreases in magnitude. Thus the finer-pair maxima
are already present before the final root closure. This excludes that final
operation as their immediate origin, but does not exclude earlier closure,
coarse-boundary interpolation, fine-stage transport or reaction history.
The next attribution must examine those earlier stages at matching cells;
endpoint closure success is not a full SDC trajectory balance.

The source call sites establish that stage1 is the Euler-predicted endpoint
and stage2 is the first reacted endpoint used by the next spatial evaluation.
`compare_application_fine_stages.py` matches these endpoint states, fine
pre/post closure, and root post closure at 100/200/300 ns while explicitly
retaining their differing interval starts. It excludes stage0 from endpoint
comparison and refuses duplicate event keys. At the finer-pair H2 maximum
(node2,3,3,8), the 300-ns difference changes from +4.03694e-8 in the Euler
predictor to -1.69713e-8 after the first forced reaction pass and -1.69489e-8
before fine closure; both fine and root closure leave it unchanged. The same
sign reversal appears at 100 and 200 ns. At the N2 maximum (node2,3,4,8),
the 300-ns predictor difference is already +1.21178e-8 and the final value is
+1.20909e-8. `nitrogen-source-audit.json` checks zero net N2 stoichiometry in
all bundle reactions; this does not remove its third-body role.

These snapshots locate the change but do not isolate chemistry: the forced
reactor integrates nonchemical forcing too. The next discriminator is actual
hydro/transport forcing and chemical increments around the first pass at
these same cells. The existing stage observer runs inside `held_mol_source`
before diffusive flux merging, so its tree flux at this call site is hydro
only; do not interpret it as total hydro-plus-diffusive flux. No solver or
acceptance changes were made for this saved-data stage analysis.


The next serial diagnostic extension records the actual forced-reactor input
and returned state/source before subsequent energy recovery. It is restricted
to held MOL subcycling with no redistribution, and defaults off. The permanent
cell-trace checker now also runs reacting enabled/disabled comparisons and
checks each recorded conserved update against forcing plus returned source.
The energy source entry is a projection correction, not heat release.
Python syntax, project contract and diff checks pass; imported LiDryer source
provenance also passes. Runtime verification of this extension is pending in
`build/eb-native-atol-trajectory32-20261024/reaction-trace-supervisor.json`.
The supervisor awaits the existing build by pidfd, ensures selected serial/MPI
targets are current, then runs the focused checker and fresh dt50/dt25 cases,
requiring ordinary outputs to match the saved trajectories byte for byte.
`audit_application_reaction_trace.py` is prepared for subsequent balance and
same-cell endpoint decomposition; it has not yet produced runtime evidence.


The reaction-trace supervisor subsequently completed successfully: selected
serial/MPI targets built, the focused checker passed all eight app invocations
(140 reacting rows, maximum scaled update residual 2.26e-16), and fresh
50/25-ns runs retained byte-identical ordinary outputs. The serial executable
SHA256 is `a245deeeeca6bbf86691faa11205dd3b7098e19ebb996762637e718b6f8f2406`;
full source identities and commands are in `reaction-trace-supervisor.json`.
The direct checker is not a rerun of the CTest configuration matrix.

`application-reaction-trace-audit.json` verifies 2016/4032 update rows, with
maximum scaled residuals 6.33e-16/3.37e-16, and 324 endpoint comparisons.
At 300 ns, first-pass H2 at probe2 has dt50-minus-dt25 mass-fraction difference
-1.69713e-8. Its exact quotient decomposition is old-state contribution
+2.69450e-7, forced increment -1.00939e-8 and returned reaction increment
-2.76327e-7. N2 at probe3 has difference +1.20984e-8, decomposed as
-1.01139e-6, +1.02349e-6 and +1.42e-16 respectively. Interval starts differ
(275 versus 287.5 ns); these are local algebraic contributions, not a
reaction-off causal experiment or a decomposition of the complete history.

`application-momentum-forcing-audit.json` joins old/new stage face fluxes and
reconstructed pressures to actual momentum updates at the six fine probes.
All 1296 component/pass samples close with maximum scaled residual 5.49e-17;
the actual forcing equals the recorded mean spatial plus external rates
exactly for these samples. The report retains directional hydro increments,
wall, spatial-minus-hydro remainder and pressure reconstruction pairs.
These probes were selected for species localization, so this does not yet
cover all transverse-momentum maxima or establish spatial convergence.
No numerical policy or acceptance tolerance was changed.


The all-cell momentum localization was then recomputed from saved 100/200/300
ns snapshots (`application-momentum-maxima.json`, with source CSV hashes).
It keeps fixed-reference error separate from dt50-minus-dt25 self-difference.
Their union contains 14 root/fine cells. At 300 ns the y-momentum reference
maximum is fine (2,5,8), -1.49586e-8 for dt50 and -8.36572e-9 for dt25;
the temporal maximum is instead fine (8,1,8), +6.87639e-4 kg/m2/s.
This location mismatch explains why the earlier species-selected six probes
do not cover the momentum question. Fresh dt50/dt25 diagnostic runs at the
14 cells are pending in `momentum-maxima-driver.json`; they require identical
ordinary outputs and subsequently audit directional momentum increments.
This is still an 8-cell temporal/reference investigation, not the open
multi-mesh spatial convergence gate.


The 14-probe dt50/dt25 runs subsequently completed (57.35/78.83 s), with
ordinary outputs byte-identical to the saved trajectories. Their 2916
momentum component/pass balances close within 1.90e-15 scaled residual;
actual forcing minus recorded spatial/external mean is zero. Records are
`momentum-maxima-driver.json` and
`application-momentum-maxima-forcing-audit.json`.

The final-root-sync comparison verifies accepted output states at all
14 probes and three clocks before computing 126 comparisons. At 300 ns,
the y temporal maximum at fine (8,1,8) decreases from +7.12884517e-4 before
root sync to +6.87639093e-4 afterward; the sync difference is -2.52454240e-5.
The x maximum at fine (8,8,8) similarly decreases in magnitude, from
-7.09928063e-4 to -6.84476528e-4. The final sync therefore reduces these
maxima; this does not exclude accumulated earlier synchronization errors.

`application-momentum-maxima-stage-comparison.json` decomposes the last
forced pass at all nine temporal maxima. At the 300-ns y maximum, the old
state difference is -1.87052608e-3, the normal y hydro increment difference
is +2.55693164e-3, x hydro is +1.18125645e-5, z hydro is +1.38004349e-8,
and spatial-minus-hydro is +1.46525856e-5. The resulting +7.12884517e-4
is unchanged between the last reactor return and pre-root-sync state.
The two interval starts are 275 and 287.5 ns, so these unequal-duration
increments cannot identify interpolation or pressure reconstruction as the
cause. Endpoint pressure/flux rates must next be aligned with pinned native
stage evidence and the original spatial-refinement maxima. No numerical
correction or convergence qualification is inferred from this closure.


### Same-clock native momentum attribution (2026-10-26)

The isolated native capture under `build/pelec-momentum-maxima-20261024/`
completed with all seven exported 1024-row tables parsed-identical to the
saved dt50 reference. Its binary SHA256 is
`24d8203a71c2491a6ab78b1967b4e0eab0054b852b1b9cfaba1abc498c618780`.
The build receipt records the pinned Diffusion.cpp source identity and reused
object/archive hashes. `compare.py` matches 2808 faces at 468 fine-cell/stage
records. Native hydro flux is total minus transport, divided by face area
and converted from CGS to SI. All sampled fine faces are fully open.

`integrate.py` compares actual candidate forcing against the native mean
old/final spatial face divergence for 12 fine intervals. At the y-momentum
reference maximum, fine (2,5,8), the final candidate-minus-native difference
is -1.49585720477e-8 kg/(m2 s), and the integrated forcing difference is
-1.49585983501e-8; the unresolved remainder is +2.63024e-14.
This establishes cumulative local update agreement with the observed error,
not a demonstrated interpolation or Riemann-solver defect.

`integrate_split.py` verifies hashes of the saved evidence and decomposes
39 component/probe results. At that y-momentum maximum the integrated
x/y/z hydro contributions are -1.88818e-11, -1.50482263416e-8 and
-4.67502e-11 kg/(m2 s); the nonhydro difference is +1.55260e-10.
Candidate wall, external and mean-forcing remainder terms are all zero at
the sampled fine updates. The nonhydro term is obtained by subtracting
independently recorded hydro from total forcing; it is not a separate
transport-coefficient attribution. The y-normal hydro contribution is the
next reconstruction/Riemann-input target. Native reconstructed face states
are not in this capture, so cell-pressure differences alone cannot isolate
that cause. All of this remains the 8-cell temporal/reference case; the
spatial-convergence gate is still open. No numerical policy was changed.


The isolated native Riemann-input capture in
`build/pelec-momentum-face-state-20261026/` completed through the notification
runner with exit 0. Seven reference output tables remain parsed-identical,
and `native_flux.csv` remains byte-identical to the preceding capture.
`compare.py` matches 72 lower/upper y-face pairs at fine (2,5,8), and checks
the native pressure reconstruction against donor pressure, acoustic slopes
and sound speed. Initial reconstructed pressures agree exactly. At the final
stage (interval start 275 ns, stage2), the upper-face candidate-minus-native
left/right pressures are +4.69514e-5/+4.89999e-5 Pa and normal momentum flux
difference is +4.23257e-5 Pa. The lower-face pressure differences are
-1.39163e-5/+2.30964e-5 Pa and flux difference is +4.13801e-6 Pa.
Candidate reconstructed velocity/density are not yet captured, so a
pressure-only causal attribution or numerical fix is not established.
The isolated source/binary and table hashes are retained in its build/run
receipts and comparison JSON. No production numerical policy was changed.


The ordinary candidate face-kinematics extension completed on 2026-10-26:
selected serial/MPI build, eight direct checker invocations, and a 55.73-s
dt50 trajectory all passed; ordinary trajectory CSVs remain byte-identical.
Binary SHA256 is `482de22385120076d26d30ffdee9cbdbcdbe4371fb365fac727b8dc20dab9cc6`.
Commands/source hashes are in `build/pelef-momentum-face-state-20261026/record.json`.
The notification runner delivered completion and the prepared comparison
then matched 72 actual pressure/density/normal-velocity/donor-sound pairs.
This is not a full CTest or installed-binary matrix.

At the final upper y face, candidate-minus-native normal velocities are
-4.75185e-8 and +5.23883e-8 m/s; reconstructed densities differ by
-3.32124e-11 and -1.30402e-11 kg/m3. The normal flux difference is
+4.23257e-5 Pa, but pressure-only attribution is not justified because all
these inputs differ. Transverse velocities and species face compositions
are not included in the comparison.

`pressure_origin.py` verifies 36 paired target-cell reconstructions using
both sides of its y faces. Initially both the cell pressure and corrections
match. At the final stage, the cell-pressure difference is +3.50238988e-5 Pa;
the target-cell reconstruction correction differences on lower/upper faces
are -1.19275355e-5/+1.19275355e-5 Pa. Their sums reproduce the measured
face-pressure differences. Thus both the incoming cell state and its
reconstruction correction contribute. This is not a same-stencil limiter
test: differing stencil states can change the correction even for identical
reconstruction operators. No tolerance or numerical policy is changed.

The same-input reconstruction replay on 2026-10-26 is recorded in
`build/pelec-momentum-stencil-20261026/replay-record.json`. The isolated native
capture exited 0, reproduced all seven 1024-row output tables, and retained
byte-identical flux and face-state traces. It saved 288 stencil rows for
72 y faces (36 stage events, two sides) at fine (2,5,8), including all nine
mass fractions, velocity components, donor sound speed and y connectivity.
The source/binary hashes and exact commands are in the build/run/replay
receipts; no production numerical source was changed for this replay.

`replay.F90` calls the production characteristic slope, positivity guard
and face reconstruction routines in their normal order, compiled with GNU
Fortran 13.3.0, `-O2 -fcheck=all -Wall -Wextra`. For all 144 left/right states,
CGS inputs reproduce the native pressure, normal velocity, density and two
acoustic slopes exactly. Repeating after SI conversion gives maximum errors
of 1.45519e-11 Pa, 5.24025e-14 m/s, 5.55112e-17 kg/m3 and 9.98507e-15
kg/(m2 s), respectively. Both unit systems satisfy the recorded roundoff
envelope of 128 epsilon times max(1, abs(native SI quantity)); the positivity
guard changes no slope. There are no replay failures or skips.

This rules out a measurable reconstruction-operator discrepancy for these
saved native stencils. It does not establish full Riemann-flux equivalence,
all-state limiter equivalence, or spatial convergence. The differing actual
candidate face states require tracing the upstream stencil/sound inputs,
including when coarse/fine preparation or prior updates first diverge.
C2 spatial and long-trajectory qualification remain open.

The candidate primitive-stencil diagnostic completed on 2026-10-26 in
`build/pelef-momentum-stencil-20261026/`: selected serial/MPI build (617.03 s),
eight direct application checker runs (5.40 s), and dt50 trajectory (56.12 s)
passed. Ordinary trajectory outputs remain byte-identical. The serial binary
SHA256 is `935ecfe275eddbcdcd0caa9cd85541129e6b56b1bb3fef88568242896dc19362`;
exact source hashes and commands are in `record.json`. The checker verifies
four offsets per face, donor sound correspondence, selected-cell pressure,
mass-fraction normalization, and unchanged enabled/disabled trajectories.
The new `.mol_stencil.csv` records prepared primitive inputs only when the
existing opt-in MOL trace is enabled. This is not an MPI trace implementation.

`compare_stencil.py` matches 288 native/candidate donor rows at the two y faces
of fine (2,5,8). Initial inputs agree to roundoff. At the first Euler predictor,
pressure differs by at most 6.64761e-7 Pa and O2 mass fraction by 4.42185e-12,
while density and velocity differences remain roundoff-sized. Sound differs by
3.64571e-9 m/s. `replay_candidate.py` reconstructs all 144 actual candidate
face pressures, normal velocities and densities exactly from these saved
inputs; the positivity guard is inactive. Together with the native-input
replay, this identifies differing inputs rather than differing reconstruction
operators for the sampled states, not for all possible states or geometries.

`predictor_origin.py` examines the first 25-ns fine Euler predictor at the
full-volume target. Native total face divergence and the candidate reactor
observer's saved `spatial_old` agree for species to 1.30207e-11 kg/(m3 s).
The inferred nonspatial predictor-source difference is +1.70674e-5 kg/(m3 s)
for O2 and -8.45925e-6 for HO2. This remainder is computed from state change
minus spatial rate, not a direct capture of the initial reaction source.
The MOL stage observer runs before transport merge, so its supplied face flux
must not be used as the total species spatial rate; the saved reactor rate
includes diffusion. No external source is present at this target.

The first pressure/composition discrepancy therefore motivates a direct
initial-reaction-source audit. Its influence on the later transverse momentum
maximum and spatial refinement behavior is not established. No numerical
policy or tolerances changed; C2 remains open. Python syntax, mechanism
provenance, project contract and diff whitespace checks passed for this work.

The 2026-10-26 initial-source audit is in
`build/pelef-initial-source-20261026/`. A standalone replay linked against the
current ordinary application libraries reproduces the inferred first fine
initial-reaction source (including O2 difference +1.70671e-5 kg/(m3 s)). The
native initial plotfile contains the direct `rho_omega_*` authority; a pinned
PelePhysics single-cell replay reproduces these rates to 8.89e-16 kg/(m3 s).
Both integrate 25 ns. Changing native analytic to finite-difference dense
Jacobian changes the native source by at most 1.008e-8, below this discrepancy.

PeleF's configured mass-fraction atol 1e-14 becomes density atol
9.6495871972e-16 kg/m3 at this cell; native 1e-18 g/cm3 is 1e-15 kg/m3.
Matching that tolerance alone retains a maximum source difference 1.70629e-5.
Diagnostic removal of constraints and use of the supplied initial temperature
also leave this difference unchanged. These were isolated replay objects;
production sources and adopted tolerances were not changed. One initial
variant compile failed because the preprocessed compile flags were applied
to raw source; rebuilding with the actual SUNDIALS preprocessor definitions
succeeded. Commands and logs are retained for each variant.

Direct evaluation at the identical initial state finds a maximum species
mass-production-rate difference 8.62632e-6 kg/(m3 s), before time integration.
The existing explicit native-activation bundle was provenance-checked and
regenerated byte-identically, then compiled into the isolated replay. This
reduces the instantaneous RHS difference to 1.05606e-10; matching density
atol as well reduces the initial-source difference to 2.33104e-9 kg/(m3 s).
The bundle changes only the documented native calorie activation conversion;
canonical Cantera chemistry remains unchanged. Tightening only the candidate
tolerance does not demonstrate reference convergence because the native
integration remains fixed at its original tolerance.

Thus the canonical/native activation convention accounts for the dominant
first-predictor chemistry discrepancy at this target. This is not a claim
that it explains the later transverse momentum maximum or failed spatial
convergence. Next compare an explicitly named native-compatible mechanism
trajectory against the canonical baseline under identical mesh/time settings.
Build/link inputs, bundle/generated hashes, replay outputs and variant
comparisons are retained in the audit directory. No full trajectory or MPI
qualification was rerun in this single-cell audit; C2 remains open.

The activation-only ordinary trajectory completed in
`build/pelef-native-activation-trajectory-20261026/` (notification exit 0).
A clean selected serial/MPI build and eight diagnostic checker invocations
passed, followed by the 55.45-s n8/dt50 trajectory through 300 ns. All seven
original reference field gates and the inventory/budget check passed.
Binary SHA256 is `3c3b10a5887a5ece80df13c8c856649747259084e907646991bb5fc9a3e8dd7c`.
The input differs from the canonical diagnostic run only in output path;
chemistry comes from the separately generated explicit native-activation
bundle, with identical Simple transport and unchanged mass-fraction atol.
Preparation and build receipts retain bundle/source/command identities.

At 300 ns the y-momentum Linf changes from 1.49585720477e-8 to
1.49581225184e-8 kg/(m2 s), ratio 0.99996995, at the same fine (2,5,8).
The x/z ratios are 1.00461316 and 0.99555902. These errors are measured over
544 composite fluid cells, not over all 1024 rows in the native level exports.
The first predictor cell-pressure discrepancy falls from +6.64746e-7 to
-2.31375e-9 Pa, but the final sampled discrepancy remains +3.22078e-5 Pa.
Thus activation conversion explains the early chemistry difference but is
not the dominant cause of this later y-momentum maximum. No compatibility
mechanism or tolerance has been adopted as the standard model.

The subsequent isolated native all-axis capture completed exit 0 in
`build/pelec-momentum-allface-20261026/`. It preserves all seven reference
tables and the previous flux/y-face/stencil traces byte-for-byte. Its 216
face-state records and 864 donor records match the actual candidate's six
faces at fine (2,5,8). With the activation-only candidate, initial valid and
halo inputs agree to roundoff. At the first predictor, z-direction halo
pressure differs by as much as 6.87414e-6 Pa, versus 2.35741e-9 Pa in valid
z donors; the halo discrepancy is unchanged on the next SDC spatial call.
This localizes an early input discrepancy to the coarse/fine halo preparation
path, but does not distinguish differing coarse donors from interpolation
operators. Next audit the pre-synchronization coarse inputs/endpoints feeding
these halo cells. These n8 diagnostics do not close C2 spatial convergence.

The target coarse-donor capture in `build/pelef-target-donors-20261026/`
completed with all seven native outputs unchanged and all candidate ordinary
and existing diagnostic CSVs byte-identical. `analyze.py` replays 144 candidate
linear interpolations and center time mixtures; 36 native donor groups are
retained. The isolated candidate build/run took 412 s end-to-end with the
native capture. No production numerical source changed.

A follow-up timestamp capture in `build/pelec-target-donor-time-20261026/`
identifies native interpolation event 2 as level 1, stage 1, requested time
25 ns and interval 25 ns. Seven outputs and the entire native donor trace
remain identical. Its first attempt failed before compilation because the
historical Advance object path differed from the current link; the retry
used the source corresponding to the actual linked object and exited 0.
The first target candidate halo uses ancestor fraction 0.5. At fine (2,5,9),
both native and candidate shared limiter scales are zero in all three axes:
the spatial interpolation returns the coarse center (1,3,5).

`eos_replay.py` and `first_halo_origin.py` replay the actual prepared states
with the current production EOS. Candidate-minus-native coarse center total
energy is -2.81483517e-6 J/m3. Direct total-energy EOS recovery gives a donor
input pressure contribution -9.85383e-7 Pa. The candidate halo path first
recovers temperature from held internal energy, then uses that guess for
hydro's total-energy recovery. Replaying these two calls adds -5.88865e-6 Pa,
reproducing the actual candidate halo pressure. The sum explains the observed
-6.87414e-6 Pa within floating-point EOS residual (recorded explicitly).
The diagnostic assertion uses 128 ulps at the pressure magnitude; this is
not a changed trajectory/physical acceptance gate.

This identifies a local temperature-recovery path discrepancy, in addition
to smaller coarse donor differences, rather than a spatial limiter defect
at this event. Next test a bounded alignment of halo temperature recovery
with the default nonrotating hydro energy authority. Do not infer that this
single-event decomposition explains the entire later momentum Linf or spatial
convergence; those require the same-condition trajectory comparison. C2 is open.

The isolated coarse-halo energy alignment pilot completed in
`build/pelef-halo-energy-alignment-20261026/` with notification exit 0.
Its n8/dt50/300-ns trajectory passes seven original field gates and budget
checks. The first z-halo pressure improves as predicted by the local replay,
but final y-momentum Linf increases from 1.49581225184e-8 to
1.50147472233e-8 kg/(m2 s), ratio 1.00378555. The x/z ratios are
1.00018414/0.98854790. Therefore this broad common-halo modification is not
adopted as a momentum-error fix. The production source remains unchanged;
no MPI or general transport qualification is claimed for the pilot. Its
local pressure improvement alone does not establish trajectory improvement.

The next bounded spatial evidence step uses the current canonical ordinary
application at n16/dt100/100 ns against the saved pinned pilot16 reference.
This fills the intermediate grid for the recorded n8/n32 one-step comparisons;
archived binary identities must remain explicit until a same-source campaign
is rerun. It is not by itself mesh convergence or a long-trajectory result.

The current canonical n16/dt100/100-ns ordinary application completed in
125.30 s in `build/pelef-mesh16-current-20261026/`, using binary
`935ecfe275eddbcdcd0caa9cd85541129e6b56b1bb3fef88568242896dc19362`.
Both initial/final original field gates and independent inventory/budget
checks passed. y-momentum Linf is 6.80369e-9 kg/(m2 s), volume-weighted L2
2.63024e-10; its maximum is fine (1,14,16), adjacent to the upper-z
coarse/fine boundary. Saved n8/n32 comparisons show increasing momentum
reference differences under refinement, but their binary identities differ.
A same-binary rerun is therefore required before interpreting this pattern.

The same-binary n8/100-ns calculation completed in 18.40 s and both field
gates passed. Its first budget check mistakenly retained the archived
300-ns H2-consumption reference while checking the 100-ns output. The driver
stopped before starting n32. The existing saved 100-ns native initial/final
exports were then verified parsed-identical to the first two exports in the
300-ns reference, and its unmodified 100-ns budget contract passed. No
calculation was repeated and no tolerance changed; original failed check
receipts are retained under `build/pelef-mesh-samebinary-20261026/n8/`.
The n32 rerun proceeds separately, saving a 100-ns checkpoint for later
continuation. This remains fixed-reference error evidence, not physical
mesh convergence or a full C2 completion claim.

The identical-binary n8/n16/n32 comparison at dt100/100 ns is now complete.
The n32 rerun finished in 951.70 s, passed initial/final field gates and
independent budgets, and saved a 94,686,586-byte checkpoint with SHA256
`cb89058862d63f0dc62089bc1c1d82c222487ffe4ac54ffd4abe0a647c6e4f45`.
All three grids use binary `935ecfe275eddbcdcd0caa9cd85541129e6b56b1bb3fef88568242896dc19362`.
`build/pelef-mesh-samebinary-20261026/mesh-record-comparison.json` records
volume-weighted L1/L2, maxima and locations. y-momentum Linf is
1.43457e-9, 6.80369e-9 and 1.86897e-7 kg/(m2 s); corresponding L2 values
are 1.33147e-10, 2.63024e-10 and 1.41912e-9. These are increasing
solver-to-reference differences, despite passing the original gates, not
physical discretization errors or evidence of spatial convergence.

To separate step-size dependence, `build/pelef-time32-half-20261026/` runs
the unchanged pinned native binary and current canonical application at
n32/dt50 through the same 100 ns endpoint. Native integration finished
normally in 76.02 s. Initial gas, geometry and temperature are identical to
the dt100 export; initial reaction histories legitimately differ because
the initialization reaction interval changes. The original driver asserted
whole-table identity and stopped after native exports, before the candidate
run. Its failure receipt is retained. The recovery supervisor validates
non-history identity, records all history differences, reuses the completed
native calculation, and starts the candidate once. No model or trajectory
acceptance tolerance was changed. Candidate gates and temporal comparison
remain pending until its terminal receipt is verified.


### 2026-10-26: n32 half-step comparison fails original field gates

The ordinary same-binary n32 run at dt=50 ns reached 100 ns in two coarse
steps (PeleF 930.433 s; pinned PeleC 76.022 s; both exit 0). The initial
field comparison and final inventory gates pass, but the 50-ns field output
fails density/xmom/ymom/total-energy/internal-energy gates, and the 100-ns
output fails xmom/ymom/total-energy/internal-energy gates. No limits changed.
The checkpoint hash was independently verified. Evidence and exact binary/input
identities are in `build/pelef-time32-half-20261026/record.json`.

At 100 ns, y-momentum Linf is 1.1593196079e-6 kg/(m^2 s), versus
1.8689706050e-7 for dt=100 ns (ratio 6.20298). Its maximum is fine
cell (2,1,27), with native volume fraction 1. The x and z momentum mismatch
ratios are 2.99090 and 3.35397. Temperature Linf improves by a factor
0.63087. These are solver differences, not physical discretization errors.
Both solvers' own y-momentum timestep differences are approximately 0.003531;
the much smaller mismatch is not a temporal-order estimate. See
`time-comparison.json` and `probe-trajectories.json` in the same folder.
The comparison script now accepts a completed failed-gate run for diagnostics,
retains its failed gate status, and continues to report qualification=false.

The initial native reaction-history columns vary with the requested timestep;
initial gas/geometry/temperature identity was verified separately. The original
pre-candidate driver assertion failure remains preserved, and the native run
was not repeated. A six-probe, same-binary diagnostic capture is running under
`build/pelef-time32-trace-20261026`, targeting current momentum maxima,
the prior dt=100-ns y maximum and a coarse donor. It must reproduce ordinary
CSV bytes before its stage/face/stencil/reaction/sync traces are interpreted.
No numerical policy has been adopted; C2 remains open.


### 2026-10-26: n32 lower-y halo velocity localized

The six-probe dt=50-ns diagnostic run finished in 896.304 s, exit 0, and all
five ordinary result/budget/step CSV outputs are byte-identical to the previous
run. Trace hashes and source identity are in
`build/pelef-time32-trace-20261026/application/run.json`. At fine (2,1,27),
the six-face hydro divergence closes to 5.83e-11 in momentum-rate units;
actual stage-2 reactor forcing integrates exactly to both saved endpoint
momentum changes. All direct synchronization increments at this cell are zero.
This does not exclude indirect coarse/fine synchronization through halo data.

An isolated pinned-native observer finished in 80.713 s with all three exports
parsed-identical to the original native outputs. See
`build/pelec-time32-trace-20261026/record.json`. Native grown-box face/stencil
captures contain duplicate records: 132 face rows collapse to 72 keys, and
528 stencil rows to 288 keys, only after asserting duplicate values identical.
The sum of native old/corrected total-flux divergence versus candidate actual
forcing reproduces endpoint y-momentum mismatch within 1e-12 at 50/100 ns.
This strengthens the forcing account without establishing interpolation cause.

Initial face fluxes match to 1.46e-11 Pa. At the first Euler predictor, the
lower-y ghost donors (2,-1,27) and (2,0,27) differ by approximately
-1.27207e-4 m/s in u and -1.27144e-4 m/s in v. Same-axis valid donor velocities
remain within 5.42e-14 m/s. The lower-y hydro flux differs by -0.00535744 Pa.
A linked production Riemann replay reproduces all captured candidate y-momentum
face fluxes and native y-normal fluxes within 2e-10 Pa. Replacing only captured
halo velocities with native values reduces the first lower-y face discrepancy
to 1.06634e-6 Pa; replacing only pressure leaves -0.00535825 Pa. Replacing all
halo primitives leaves 2.61760e-7 Pa. See `replay-comparison.json` in the native
trace directory. These are local counterfactual replays, not adopted numerical
changes or proof of a whole-trajectory fix.

A bounded first-50-ns donor capture is now under
`build/pelef-time32-donors-20261026`: candidate interpolation center, neighbors,
time mixtures and result; native interpolation donors with explicit fill time,
level and stage. It must reproduce saved native outputs and candidate first-step
outputs/trace prefixes before separating donor-input and limiter differences.
No production numerical change is adopted. C2 remains open.


### 2026-10-26: shared limiter amplifies a trace-species donor difference

The bounded donor capture completed: native outputs at0/50ns remain identical,
and candidate first-step outputs and four trace prefixes match the saved
canonical run. Candidate runtime was521.415s. Evidence is under
`build/pelef-time32-donors-20261026` (timed and candidate receipts).
Twenty-four candidate interpolation replays pass. Native fill context explicitly
assigns donor event2 to25ns, fine-level predictor; candidate capture5 uses
coarse-time fraction0.5. `origin.py`/`origin.json` reconstruct the actual halo
velocities for both solvers, with native raw-CGS replay exactly matching its
captured velocity. Converting native donors to SI before subtraction changes
velocity by only3.08e-10m/s; this is retained as replay arithmetic sensitivity.

At that event the common z-direction slope factor is controlled by H2O:
PeleF0.4188983387335005 versus nativeCGS0.4245690822381488
(native converted-to-SI0.4245690959673618). The z-neighbor/center H2O densities
are nearly flat, around4.4623e-9kg/m^3; tiny differences in those gradients
change the shared factor, which also scales momentum slopes. Replacing only
the candidate center/six-neighbor H2O values with native values reduces the
u/v halo residuals from approximately1.27e-4m/s to4.21e-10/-3.67e-11m/s.
The interpolation replay therefore identifies trace-species input sensitivity,
not a demonstrated error in the interpolation formula. It does not yet assign
the upstream discrepancy solely to chemistry or justify changing the limiter.

A first50ns n32 activation-basis-only comparison is running in
`build/pelef-time32-activation-20261026`, using the previously provenance-checked
native-activation+Simple executable (SHA3c3b10a5887a5ece80df13c8c856649747259084e907646991bb5fc9a3e8dd7c).
The timestep, tolerances and numerical policies remain the canonical ones.
The completed n8 activation comparison is insufficient to infer the n32 result.
`analyze.py` prepares endpoint and first-predictor comparisons and recomputes the
native50ns H2 reference for unchanged inventory gates after completion.
No production numerical change is adopted; C2 remains open.


### 2026-10-26: activation basis does not resolve n32 limiter sensitivity

The activation-only n32 first50ns run completed in 523.607s, exit0.
Initial fields and correctly recomputed50ns inventory gates pass. Final fields
still fail density/xmom/ymom/total-energy/internal-energy gates. Y-momentum Linf
changes from1.7980646891e-6 to1.7967591189e-6 (ratio0.999273903), while H2O Linf
falls from2.68615e-13 to5.34493e-15kg/m^3. The first predictor lower-y halo v
difference remains-1.27051e-4m/s versus canonical-1.27144e-4. See
`build/pelef-time32-activation-20261026/comparison.json`. Thus activation-basis
alignment alone does not remove the relevant gradient/limiter discrepancy.

Three saved coarse initial states (1,32,13..15) were replayed without CFD under
`build/pelef-time32-donors-20261026/initial_replay`. Initial H2O source is exactly
uniform across these three cells for each variant. Native source is
0.178499550604kg/(m^3s); canonical0.178494203091; native activation0.178499450993;
native activation with matched density atol0.178499550534. This is initial
source evidence only, not a replay of the forced coarse-step endpoint that
sets the interpolation slope.

An isolated first50ns pilot now fixes full-coordinate CVODE density atol to
1e-15kg/m^3, exactly the native1e-18g/cm^3, on the native-activation variant.
Production uses the existing mass-fraction tolerance converted with density.
Only an ignored-build copy of `sundials_forced_reactor_mod.F90` is changed;
the reduced-coordinate routine and production source remain untouched. Build,
commands, hashes and application results are recorded under
`build/pelef-time32-density-atol-20261026`. This is a convention-matching test,
not a changed field gate or an adopted numerical policy. `analyze.py` checks
the same endpoint/halo and inventory metrics after completion. C2 remains open.


### 2026-10-26: matching density atol is not sufficient

The isolated native-activation plus native density-atol run reached50ns in
518.785s, exit0. Initial fields and inventory pass; final reference fields
remain failed. Y-momentum Linf is1.7913107905e-6, ratio0.996243796 to the
canonical result, and first predictor lower-y halo v residual remains
-1.2666605991e-4m/s. Evidence:
`build/pelef-time32-density-atol-20261026/comparison.json`.
Neither activation conversion alone nor its combination with exact density
atol explains the dominant limiter sensitivity; neither variant is adopted.

The next bounded capture, `build/pelef-time32-coarse-source-20261026`, records
coarse cells(1,32,13..15) before/after predictor and reaction updates, including
native transport/total face fluxes and conserved inputs, and candidate actual
forced-reactor inputs. This follows the H2O gradient formation upstream of
interpolation. A prepared isolated replay compares candidate actual forcing
with native forcing inferred from old/corrected total face divergence. Such
inferred forcing must be checked against stage and endpoint evidence before
being treated as native reactor parity. No production source change is made.


### 2026-10-26: transport thermodynamic recovery identified upstream

The coarse-source capture completed in521.261s; ordinary first-step CSVs are
identical. Native0/50ns exports are identical. Twenty-four isolated forced
reactor replays were performed; canonical actual forcing reproduces all six
candidate coarse H2O outputs exactly. At the final coarse iteration, using
native inferred total-face forcing changes the H2O z limiter from0.419062
(native activation with candidate forcing) to0.424544; adding exact native
density atol gives0.424569389, close to native0.424569096 after SI conversion.
The source/forcing evidence is in
`build/pelef-time32-coarse-source-20261026/replay-comparison.json`.

Initial coarse energy-rate discrepancies for z=13,14,15 are dominated by
transport: -512.814978,-664.386994,-780.395196W/m^3. Saved-field production
face-transport replays using held internal-energy EOS reproduce those three
actual discrepancies within1.9e-9W/m^3. The same replay using total energy
minus kinetic energy reduces the native transport-rate residuals to about
9.21e-6W/m^3. Momentum gradients, transport fits and flux formulas are the
same in these paired replays. See `transport-held-comparison.json`,
`transport-total-comparison.json`, and `analysis-provenance.json`.

The periodic y neighbor is coarse(1,1,14), covered by fine cells. Its native
initial total energy minus stored internal energy minus resolved kinetic
energy is0.017592128289J/m^3; at the target coarse(1,32,14) the same residual
is only-8.12e-12J/m^3. Averaging fine conserved states therefore matters to
thermodynamic recovery at this interface. Native `pc_ctoprim` recovers T from
total-minus-kinetic energy, whereas the candidate transport call supplies
held internal energy. At the target's upper-y face, total-EOS replay matches
native energy flux to2.30e-9W/m^2, including thermal and species contributions.
This is an evidenced recovery-policy mismatch upstream of the H2O slope
amplification, not evidence that the shared interpolation formula is wrong.

A single-call transport-total-energy pilot is running to100ns in
`build/pelef-time32-transport-total-20261026`. It changes only an ignored-build
copy of the serial tree transport call; canonical mechanism, original reaction
tolerances and original field/inventory gates are retained. If confirmed,
production implementation needs an explicit ordinary nonrotating MOL path in
serial and MPI, preserving existing held recovery for general transport,
rotating frames and experimental halo policies. MPI transport lives in
`mpi_eb_tree_transport_3d_mod.F90`; merely fixing the serial call is insufficient.
The pilot is not adoption, full matrix evidence, or C2 completion.


### 2026-10-26: transport-total pilot passes; scoped integration in progress

The isolated canonical-mechanism transport-total pilot reached100ns in891.278s,
exit0. All original field gates at0/50/100ns and final inventory gates pass.
At50ns, momentum Linf(x,y,z) is(6.48686e-8,1.20778e-8,1.46195e-8); at100ns it
is(4.11096e-8,1.40471e-8,2.64764e-8)kg/(m^2s). The canonical100ns y mismatch
was1.15932e-6. Evidence:
`build/pelef-time32-transport-total-20261026/application/run.json`.
This verifies the recovery-policy correction for the bounded serial case;
it does not establish mesh convergence, long trajectories, MPI parity or C2.

Production integration adds the optional `primitive_from_total_energy` to
serial/MPI tree transport-halo kernels, defaultfalse. Only ordinary MOL callers
without a rotating frame or separate-momentum policy request it. Normal-face
transport then omits the held-energy override; halo interpolation and existing
wall-source recovery retain their inputs. Explicit requests with a rotating
frame or separate-momentum policy are rejected. MPI includes the new flag in
collective control agreement before mutation.

A serial regression supplies a deliberately inconsistent held thermal field:
total-energy recovery must match an independent no-held reference, while the
default-held result must differ. MPI extends the ownership/depth comparisons
with the new policy and checks rank-disagreement rollback. Build and43 focused
registered tests are queued before rebuilding both selected applications and
running the8-case diagnostic checker. Results are pending, not passed.

Baseline copies and `task-only.patch` under
`build/pelef-transport-total-integration-20261026` isolate this change from the
pre-existing dirty checkout. Provenance preflight, Python syntax,
`git diff --check`, and project contract pass before the build. No compiler,
MPI, backend configuration or generated mechanism source was changed.


### 2026-10-26: first integration test result, 42/43 pass

The focused unit/application-kernel build succeeded in244.757s. The43-test
focused CTest selection ran in3465.634s:42 passed and1 failed.
`unit_eb_tree_mol_spatial_halo_3d` fails its existing exact-equality assertion
`retained energy output preserves tree gas rate` in the hydro-only observer
fixture. The new serial transport recovery regression, MPI total-energy
ownership/depth comparisons, rank-disagreement rollback, and rotating/non-MOL
controls passed. No failed test is counted as passed and no tolerance changed.

The saved pre-change serial tree module is being compiled into an isolated
copy of the failing test to establish whether this failure predates the
transport change. Its evidence will be under `baseline/record.json` in
`build/pelef-transport-total-integration-20261026`. Selected serial/MPI application
build and diagnostic checks continue independently under `selected-record.json`.
The main `record.json` remains failed and its original logs are preserved.
This is not a clean full matrix or completed release qualification.


The isolated pre-transport-change serial tree reproduces the same hydro-only
bitwise assertion failure (exit1,2.724s). The source hash, exact compile/link
commands and failure log are retained under
`build/pelef-transport-total-integration-20261026/baseline`. Therefore that
failure predates this transport change; it remains an unresolved failed check,
not a passing test. The selected serial/MPI application build continues, and a
fresh nonuniform MPI2 checkpoint to MPI3 restart comparison against serial
through200ns is prepared under the same integration directory. Its scope is
restart equivalence, not fixed-reference200ns qualification.


### 2026-10-26: selected applications and restart pass

The selected serial/MPI build completed in710.110s. All8 diagnostic checker
cases pass in6.666s, including ordinary-output neutrality and refusals. The
new nonuniform MPI2 save100ns / MPI3 resume200ns test passes: prefix and final
CSV files are byte-identical to uninterrupted serial, and the saved checkpoint
hash is unchanged. Actual final bytes and checkpoint digest were independently
rechecked. Records: `selected-record.json`, `trace-checks/`, `restart-check/run.json`
under `build/pelef-transport-total-integration-20261026`.

Production n32 dt50ns to100ns is running under
`build/pelef-time32-transport-production-20261026`, with original gates and
comparison to the successful pilot's ordinary bytes. A same-current-binary
n8/n16/n32 dt100ns to100ns campaign is running under
`build/pelef-mesh-transport-total-20261026`. Its n8 native reference uses the
correct100ns smoke contract. Neither pending calculation is reported passed.
The43-test record remains42passes plus1 pre-existing hydro bitwise failure.


The same-binary mesh campaign completed n8 and n16 (18.583s and127.153s),
with both initial/final field gates and inventory gates passing. Its bookkeeping
then attempted to hash an n16 checkpoint that the unchanged input never requested.
The original failed driver receipt is preserved as `initial-driver-record.json`.
`resume_n32.py` verifies existing receipts and the binary hash, records the n16
checkpoint as not requested, and runs only the previously unstarted n32 case.
No successful n8/n16 calculation was repeated. This is a driver bookkeeping
failure, not a failed simulation or a checkpoint-format failure.


### 2026-10-26: production n32 dt50 reproduces successful pilot

The rebuilt production serial application (SHA
`9fa77a2f476637cca9e1f665fcf9bee0a7bd90a6706c02743fd9033261a8f33f`)
reached100ns in971.707s, exit0. All original0/50/100ns field gates and final
inventory gates pass. Five ordinary result/budget/step CSVs are byte-identical
to the successful isolated transport-total pilot. Actual file digests and
100ns output clock were independently checked, and the94MB checkpoint digest
was recorded in
`build/pelef-time32-transport-production-20261026/verified-result.json`.
Thus the scoped production serial change reproduces the pilot's improvement;
this is bounded reference agreement, not a convergence or physical-validation
claim. MPI restart evidence remains separately recorded. The n32 dt100 mesh
campaign continuation is pending, and the pre-existing hydro bitwise unit
failure remains unresolved.


### 2026-10-26: same-binary mesh audit passes gates, not all maximum norms

All n8/n16/n32 dt100ns cases reached100ns and passed initial/final field gates
and inventory gates. n32 runtime978.198s. In
`build/pelef-mesh-transport-total-20261026/mesh-record-comparison.json`,
y-momentum fixed-reference Linf is1.20252e-9,2.02325e-9,2.90961e-8kg/(m^2s).
The corresponding pre-change n32 value was1.86897e-7. Fixed-reference maximum
errors therefore improve but do not monotonically decrease with refinement.
Atn32 the new dt50ns y mismatch is1.40471e-8, about0.483 of dt100ns.

A separate conservative volume restriction compares each solver to itself
across adjacent grids; matching fluid volumes is explicitly checked. See
`mesh-self-comparison.json`. From8/16 to16/32, candidate y-momentum Linf/L2
ratios are0.305273/0.232783, closely matching native0.305273/0.232783.
All six inspected fields' volume-weightedL1/L2 differences decrease. However,
z-momentum Linf grows by1.909647 in both solvers. These inter-grid differences
are not exact physical errors, and the fixed100ns timestep plus two grid pairs
do not establish asymptotic spatial order or close C2.

A same-binary n32 dt25ns to100ns comparison against a fresh fixed-native run
is now under `build/pelef-time32-quarter-20261026`. It retains original field
and inventory limits and separates initial gas identity from timestep-dependent
reaction history. Five output times and a three-timestep self-comparison of
conserved fields/species and temperature are prepared for verification.


### 2026-10-26: third timestep and evolution-only mesh comparison

The n32 dt25ns run completed through100ns: native78.866s and candidate976.237s.
All five original field comparisons and the inventory gate pass. Receipts and
three-timestep comparisons are in `build/pelef-time32-quarter-20261026/record.json`
and `time-comparison.json`. The same production binary is used for100/50/25ns.
All15 fields have decreasing volume-weighted temporal L1/L2 differences.
For y momentum, the ratio of50-vs25 to100-vs50 differences is0.422803/0.351357
in Linf/L2. For z momentum it is1.138084/0.494817: its maximum is nonmonotone.
Pinned PeleC reproduces these trends. These are temporal self-differences,
not exact errors or an asymptotic-order qualification.

The mesh analysis now also compares final-minus-initial evolution increments
before conservative restriction, separating initial-profile sampling from
subsequent evolution. `build/pelef-mesh-transport-total-20261026/self_evolution.py`
and `mesh-evolution-comparison.json` cover all15 fields. From8/16 to16/32,
y-momentum evolution Linf/L2 ratios are0.948318/0.356908, closely matching
native values. Thus the endpoint Linf ratio0.305273 above overstates the
reduction attributable to evolution alone. All15 evolution L1/L2 differences
decrease; z-momentum Linf still grows by1.909647 in both solvers. This does
not establish all-field spatial convergence.

A n32 dt25ns extension uses the accepted100ns checkpoint to reach200ns in
`build/pelef-time32-extend200-20261026`. The fresh fixed-native run has finished
and its0..100ns prefix matches the accepted reference byte-for-byte. Candidate
restart, five100..200ns field gates, cumulative inventory, and preservation
of the source checkpoint remain to be verified. The archived initial CSV is
copied with provenance solely for the cumulative budget; it is not a new
candidate output. C2 remains open, including the pre-existing hydro unit failure.


### 2026-10-26: 200ns extension fails field gates despite budget closure

`build/pelef-time32-extend200-20261026/record.json` records normal completion
(native151.241s, candidate868.747s), but scientific field gates fail. The five
100/125/150/175/200ns comparisons pass only through150ns. At175ns both energy
fields exceed their original scaled limit; at200ns density, z momentum, total
and internal energy fail. The200ns z-momentum maximum is3.27660658e-6kg/(m^2s)
(limit1e-6), and total-energy scaled difference is6.15009016e-8 (limit1e-8).
No thresholds were changed. Cumulative inventory passes. A launcher exit0
therefore does not constitute reference agreement.

Independent `localize.py`/`localization.json` verification confirms source
checkpoint hash unchanged, restored100ns CSV byte-identical, final checkpoint
hash matching its receipt, and all five native prefix exports byte-identical
to the accepted25ns run. It maps the ten largest absolute differences for
rho, all momentum components and total energy at each saved time. All250
selected entries are fine-level k=32 with volume fraction1. The200ns z and
energy maximum is fine(26,14,32); y momentum peaks at(25,14,32). These are
upper coarse/fine interface cells, not cut cells. This location evidence does
not yet attribute the remaining mismatch to interpolation or synchronization.

`build/pelef-time32-extend-trace-20261026` now repeats only the100..200ns
checkpoint segment with six relocated probes (three fine maxima, another
fine maximum, and two parent cells across the upper interface). Ordinary
outputs must match the uninstrumented extension before interpretation.
`build/pelec-time32-extend-trace-20261026` separately builds ignored diagnostic
copies of Diffusion.cpp/MOL.cpp to capture normal-face reconstruction and
all-direction transport/total fluxes; all nine exports must match the unchanged
fixed native reference. Pending analysis compares actual total forcing,
directional flux divergence and reconstructed/stencil pressure and velocity.


The upper-interface native diagnostic completed in195.423s. Independent
`build/pelec-time32-extend-trace-20261026/native_audit.py` verifies all nine
exports byte-identical and the fixed reference executable unchanged. Recorded
old/corrected total face divergence integrates y/z momentum changes at the
three fine probes over100..200ns with maximum residual9.01e-16kg/(m^2s).
Overlapping grown-box captures have identical duplicates:528 face rows reduce
to288 unique entries and2112 stencil rows to1152. No conflicting duplicate
was accepted. Candidate trace identity and cross-solver diagnosis remain pending.
A production Riemann replay helper is rebuilt from current selected libraries
for captured-halo counterfactuals; no counterfactual result is claimed yet.


### 2026-10-26: upper-interface stage balance isolates transport contribution

The candidate trace completed in868.215s. Seven ordinary output files are
independently byte-identical and the source100ns checkpoint is unchanged.
`build/pelec-time32-extend-trace-20261026/comparison.json` closes the y/z
momentum error growth at three fine probes using actual candidate forcing
minus native old/corrected face divergence, with residuals around1e-15.
At fine(26,14,32), all144 recorded direct z-momentum sync increments and
all16 reaction rates are zero. This excludes a direct local synchronization
increment, not an indirect effect through coarse halo donors.

`transport-residual.json` separates the100..200ns z-momentum error growth
(-3.27211458e-6kg/(m^2s)) into hydro(+5.53773089e-7) and the non-hydro residual
(-3.82588767e-6). The first large change occurs over137.5..150ns: total
-1.53769961e-7 versus hydro+1.00096811e-8 and residual-1.63779642e-7.
Thus the pressure/Riemann contribution alone neither explains the sign nor
the dominant magnitude. Direct directional transport capture is needed.

`stencil-comparison.json` and `replay-comparison.json` validate production
Riemann replay against the recorded face fluxes. Replacing upper-halo pressure
can nearly remove some early upper-z-face differences, but later interior
pressure/velocity differences persist, and this counterfactual does not close
the dominant transport residual. No pressure-only correction is adopted.

The bounded next capture is `build/pelef-upper-transport-20261026` (saved100ns
checkpoint to150ns only) and `build/pelec-upper-transport-20261026` (fixed native
prefix through150ns). Ignored diagnostic copies record six transport face
fluxes and the local24:28,12:16,30:35 donor neighborhood. Candidate outputs
must match accepted steps4..6 and native exports must match steps0..6.
Production EOS and face-transport replay helpers are prepared to separate
velocity gradients, temperature-dependent coefficients and directional fluxes;
replay residuals must be checked before interpreting those counterfactuals.


The native upper-transport run completed in133.661s; all seven exports match
fixed-reference hashes. Its5400 local primitive rows cover36 fine stages.
`build/pelec-upper-transport-20261026/replay_native.py` feeds these native
primitives and fully-connected central gradients into the current production
PeleF face-transport kernel. Across216 faces, maximum momentum flux residuals
are8.11e-15Pa or less; energy flux residual is1.06e-9W/m^2. This same-input
local agreement supports inspecting stage inputs/gradients before changing
the transport formula. It is not a coupled-trajectory qualification. Candidate
transport capture and same-input replay remain pending.


### 2026-10-26: transport replay identifies upper-halo tangential gradients

The isolated candidate transport capture completed in419.742s; steps4..6
are byte-identical and the original100ns checkpoint is unchanged. Both local
transport replays pass the source-reproduction check: candidate momentum flux
residual below4.61e-15Pa, native below8.11e-15Pa. Energy flux residuals are
9.10e-13 and1.06e-9W/m^2, respectively. See
`build/pelef-upper-transport-20261026/comparison.json`.

At fine(26,14,32), first-predictor native event35/candidate transport event11
uses the150ns ancestor endpoint. Upper halo(26,14,33) has candidate-minus-native
(u,v,w)=(-4.36555e-4,+1.18014e-3,+7.43e-9)m/s. Its adjacent x-child(25,14,33)
has(+4.37637e-4,-1.18306e-3,-2.52e-8)m/s. The tangential velocity split changes
transverse gradients at the upper-z face and yields z-momentum transport-flux
difference+0.00397061470Pa. Using native transverse gradients in the replay
reduces that difference to-5.15e-7Pa; native temperature alone leaves it near
+0.00397066418Pa. These are captured-input counterfactuals, not trajectory fixes.

`balance.json` independently integrates all six transport faces over137.5..150ns.
The x/y/z directional contributions to z-momentum error growth are
-1.62503e-9/-4.70108e-9/-1.57454e-7kg/(m^2s). Their sum matches the previously
inferred non-hydro residual within4.73e-16. Thus the upper-face transport
contribution is now directly observed, rather than inferred only by subtraction.

`build/pelef-upper-donors-20261026` captures candidate coarse center/six
neighbors and timed native FillPatcher donors for coarse(13,7,17), targeting
fine(25/26,14,33). The native request time and candidate interpolation replay
must be established before assigning the shared limiter's controlling component.
No production numerical correction or tolerance change is made at this stage.


The timed native donor run completed in131.150s; all seven exports match.
Native FillPatcher events35/36 request150ns (fine dt12.5ns, stages1/2).
`build/pelef-upper-donors-20261026/native_origin.py` replays raw-CGS donors at
event35 and reproduces actual transport-halo velocities at(25/26,14,33) within
2.23e-16m/s. The native shared x-limiter factor is0.495740819761 and is
controlled by z momentum, not a trace species. y/z factors are zero. The
candidate controlling component remains to be verified from its capture;
the earlier trace-species mechanism cannot be assumed for this location.


### 2026-10-26: shared limiter amplifies a synchronized coarse momentum difference

The candidate donor capture completed in426.197s and matches accepted steps4..6
and ordinary stage traces. `build/pelef-upper-donors-20261026/analyze.py`
reproduces48 interpolation captures. `origin.py` matches candidate capture42
to the independently recorded event11 transport state at(26,14,33), with
ancestor fraction1, and the timed native150ns donor. Both use a common
x-limiter controlled by z momentum: candidate0.491724909679 versus native
0.495740819761. Replacing only all seven donor z-momentum values by native
values reduces the halo velocity residual to3.83e-10m/s or less. Replacing
H2O or all species leaves the velocity error unchanged. Native SI and raw-CGS
replays agree here; this is not the previous trace-species limiting mechanism.

The x-neighbor/center/x-neighbor z-momentum donor errors at150ns are
-4.78953e-9/-9.12281e-10/-1.82945e-10kg/(m^2s). The left donor is coarse
(12,7,17), just above the child. In `donor-timeline.json`, its error changes
from+4.20031e-10 in the125ns pre-synchronization endpoint to-6.52812e-9 in
the next interval's125ns starting state. This motivates checking the accepted
coarse/fine face transfer, rather than altering the limiter to hide its input
difference. The donor substitution is a local sensitivity result, not a fix.

`build/pelef-upper-sync-20261026` runs only the accepted100..125ns checkpoint
segment with one coarse and four child-face probes. Its matching native
diagnostic is `build/pelec-upper-sync-20261026`; six0..125ns reference exports
must remain identical. The next comparison reconstructs the synchronization
increment from integrated coarse and four fine face transfers, then compares
it with the observed donor jump. No production numerical changes were made.


The upper-sync native diagnostic completed in111.254s with all six exports
identical. `build/pelec-upper-sync-20261026/native-sync.json` reconstructs
coarse(12,7,17)'s125ns momentum correction from the accepted lower-z coarse
transfer and four upper-z child transfers: -1.98192524323e-5kg/(m^2s), matching
the timed donor change within2.78e-16. Candidate synchronization comparison
remains pending; this establishes the native transfer convention first.


### 2026-10-26: synchronization closes; incoming fine transport differs

The canonical upper-sync candidate completed in213.880s. Independently checked
steps4/5, trace hashes and source checkpoint agree with the accepted extension.
`build/pelec-upper-sync-20261026/comparison.json` reconstructs its125ns coarse
correction as-1.98262005799e-5kg/(m^2s), exactly matching both the reported
increment and donor change. Native correction is-1.98192524323e-5; their
-6.94814761e-9 difference agrees with the observed donor error jump within
2.78e-16. Therefore the captured synchronization applies each solver's incoming
fluxes consistently; a synchronization formula defect is not demonstrated.

The coarse integrated-flux difference is7.67e-14, versus fine face differences
[-5.41635e-12,+8.39784e-13,-3.95917e-12,+1.57394e-13] in SI time-integrated
momentum flux density. `flux-split.json` shows the two dominant negative fine
differences are mainly transport residuals (-4.63163e-12 and-4.19796e-12).
This completes an observed amplification path from fine transport transfer,
through coarse synchronization and z-momentum-controlled interpolation, to
larger tangential velocity gradients and subsequent viscous stress differences.
It does not yet identify the original source of the small incoming discrepancy.

A bounded hypothesis test now runs in `build/pelef-upper-activation-20261026`:
clean isolated serial build from current source, previously source-verified
native activation-energy basis (ratio0.9999943013061792), unchanged Simple
transport and chemistry tolerances, n32/dt25ns fresh initial state to150ns.
The activation-only transformation and identical Simple transport data were
reverified. All seven original fixed-reference field gates, final inventory,
coarse donor errors and upper-interface maxima will be compared with current
canonical results. This revisits the known conversion mismatch after the
transport-EOS correction; the earlier pre-correction activation pilot did not
resolve its dominant error. The current test is diagnostic, not production
mechanism adoption or C2 closure. Generated data and the fixed native solver
remain unchanged.


### 2026-10-26: activation-only sensitivity does not explain the target error

The clean diagnostic build completed in422.818s and its0..150ns run in1398.568s.
All seven original field gates and inventory pass, as they do for the canonical
run over this interval. `build/pelef-upper-activation-20261026/comparison.json`
compares actual saved fields: at150ns y-momentum Linf rises from1.01148448e-7
to1.29873615e-7 (ratio1.283990), while z Linf falls from4.47214338e-7 to
2.85650473e-7 and changes location. The target fine(26,14,32) z error worsens
from-1.57310e-7 to-2.19920e-7; relevant coarse donor errors also worsen.
Therefore this change is not adopted as a resolution of the localized mechanism.
Binary/checkpoint hashes and inventory receipt were independently checked.

`build/pelef-upper-density-atol-20261026` now isolates the other known solver
convention: canonical mechanism, unchanged rtol1e-10, but full-coordinate
absolute density tolerance fixed at1e-15SI, matching pinned ode.atol1e-18CGS.
The ordinary caller otherwise supplies input atol1e-14 times local density
(`src/reactive/reactive_3d_mod.F90`). This is an ignored-source diagnostic
replacement only, with fresh initial state and the same150ns field/budget gates.

In parallel, `build/pelef-initial25-trace-20261026` and
`build/pelec-initial25-trace-20261026` capture only the first25ns. The original
z-momentum maximum is already9.25980e-8 at fine(1,8,32) after that first step,
so identifying this earlier seed avoids assuming the later amplification is
the original cause. Existing0/25ns exports provide the identity baseline.


### 2026-10-26: first-step seed precedes the later z-momentum limiter

The initial25ns candidate capture completed in331.895s and independently
matches the canonical steps0/1 and trace hashes. Native capture completed
in34.544s with both exports identical. At fine(1,8,32), initial stencils agree
at roundoff, whereas the first predictor's upper halo differs in tangential
velocity by2.63517e-4m/s before appreciable interior velocity differences.
The two fine substep z-error increments are3.13893e-8 and6.12087e-8kg/(m^2s).
Their hydro contributions are negative and much smaller; the non-hydro
residual dominates. Actual forcing closes the final9.25980e-8 difference.
See `build/pelec-initial25-trace-20261026/initial-balance.json` and stage reports.

Timed native donor capture in `build/pelef-initial25-donors-20261026` completed
in31.963s with both exports identical. Raw-CGS replay at requested12.5ns
(event2) reproduces the actual fine(1,8,33) halo velocity exactly. Here the
shared y-limiter is0.315913188863 and is controlled by H2; x/z factors are zero.
This earlier seed must be distinguished from the150ns z-momentum-controlled
x-limiter. Candidate donor replay and the independent density-atol sensitivity
remain pending; neither chemical causality nor a numerical fix is claimed.


### 2026-10-26: H2 seed separates reaction offset from forcing variation

The initial donor candidate completed in322.265s, matches steps0/1 and stage
traces, and all12 captured interpolations replay. At requested12.5ns the
shared y-limiter is0.316618951344 versus raw-native0.315913188863. Substituting
only H2 donors reduces the tangential velocity discrepancy from2.63517e-4
to1.08313e-7m/s. The remaining difference includes SI-conversion cancellation:
native donors replayed after SI conversion yield factor0.315913478952,
whereas raw CGS replay matches the actual native halo exactly. See
`build/pelef-initial25-donors-20261026/origin.json` and `native-origin.json`.

The density-atol-only coupled run completed in1440.204s and passes all seven
150ns field times and inventory. Its y-momentum Linf ratio to canonical is
0.998649 (only0.135% reduction); z ratio0.745172. It does not resolve the
principal y discrepancy and is not adopted as a demonstrated fix.

`build/pelef-initial25-source-20261026` captures native coarse(1,3/4/5,17)
old/corrected transport and total face fluxes. Saved candidate coarse-center
forcing reproduces the actual candidate H2 exactly in the production reactor
replay. At the final coarse endpoint, replacing all forcing components by
native-derived values changes the three canonical H2 results by approximately
[-7.11e-17,-6.68e-17,+1.349e-15]kg/m^3 relative to their actual values. This
nonuniform forcing effect matters to the nearly-flat H2 slope, independently
of the larger almost-uniform chemical offset.

`triplet.json` verifies that the native reactor using the inferred forcing
reproduces actual native endpoint H2 within1.74e-18kg/m^3. Candidate production
reactor with native activation basis and native density atol reproduces all
three within7.81e-18 when supplied the same native forcing. Canonical chemistry
with the same forcing retains an approximately uniform2.334e-14 offset.
These are local controlled replays, not coupled-CFD qualification.

`build/pelef-initial25-matched-chemistry-20261026` now tests both conventions
together for only0..25ns with original gates and budget. In parallel, the
candidate source capture extends the recorded coarse probes to all three
cells, reusing the already captured native data. This will distinguish the
actual nonreacting forcing variation from the reactor-only offset. No baseline
mechanism, production tolerance or acceptance limit has been changed.


### Initial 25 ns three-donor forcing and combined chemistry decision

The canonical diagnostic capture in `build/pelef-initial25-source-20261026/application`
completed in 316.187 s, exit 0, with step0/step1 byte-identical to the quarter-step
reference candidate. `replay_all.py` completed 30 isolated forced-reactor replays;
all six canonical actual-forcing H2 endpoints exactly reproduce the application.
At coarse y=3/4/5, replacing the full forcing vector with native face-derived
forcing changes corrected H2 density by -7.11237e-17/-6.67869e-17/+1.34875e-15
kg/m3. This is a full forcing-vector counterfactual, not species-only attribution.

`source_split.py` records old/corrected average source decomposition and directional
hydro divergences in `source-split.json`, with input hashes. At coarse(1,5,17),
corrected-average H2 forcing difference is -5.38282559e-8 kg/m3/s: hydro
-5.29494208e-8 and inferred transport -8.78835033e-10. At the first predictor,
z-direction hydro divergence difference is -5.90315446e-8; x/y differences
are +7.37010e-11/-1.01533e-10. The next diagnostic captures the native
reconstructed faces and primitive stencil there. Cancellation in directional
flux divergence requires component-scaled closure checks; absolute energy
closure residual alone is not a momentum/source mismatch.

The separate `build/pelef-initial25-matched-chemistry-20261026` native-activation
plus native constant-density-atol sensitivity passes both original field times
and budget gates. At25ns, y-momentum Linf changes from1.684333998e-8 to
1.143346862e-8; z Linf from9.259795817e-8 to9.199117186e-8. Thus chemistry
alignment does not resolve the initial dominant z maximum. These diagnostics
do not adopt a mechanism, limiter, tolerance or synchronization change and
do not close C2.

The native coarse-face capture `build/pelec-initial25-coarse-face-20261026`
completed in30.182s; both exports have identical SHA256 to the pinned reference.
The first-predictor z-face comparison finds pressure differences about6.986e-7Pa
on the coarse-cell side, versus velocity differences below2.5e-13m/s and density
differences below7e-17kg/m3. `face-comparison.json` retains the comparison.
Pressure substitution in a same-production face replay is the next causal check;
this observation alone does not establish a reconstruction defect.


### Coarse pressure substitution and conservative-state origin

`build/pelec-initial25-coarse-face-20261026/replay.py` replays the captured
four-cell primitive stencils through the production Riemann kernel, verifying
candidate all-component fluxes and native H2 fluxes before substitutions.
At first-predictor coarse(1,5,17), z-direction H2 divergence discrepancy changes
from -5.90311757e-8 to -6.23064100e-13 kg/m3/s when only stencil pressures
are replaced by native values. Velocity, density/sound, and species substitutions
leave this discrepancy essentially unchanged. Corrected-stage pressure substitution
reduces x discrepancy from -3.58308903e-8 to -5.68327607e-11, and z from
-6.99235723e-8 to -2.42488387e-9. These are local counterfactuals, not trajectory
qualification or proof of a reconstruction implementation defect.

`build/pelef-initial25-source-20261026/pressure_origin.py` compares EOS
pressure from captured conservative states. With tighter inversion used only
as a diagnostic, the predictor pressure difference is6.98724762e-7Pa.
Replacing all nine species densities leaves -1.45519e-11Pa; replacing density,
momentum, or total energy separately leaves the original difference. Thus
the species-state discrepancy explains the pressure difference at this cell.
The pinned default inversion stops before the final update at1e-6K; changing
the initial temperature guess can mask this small difference. The old stored
temperature reproduces the candidate actual pressure exactly; a shared converged
inverse confirms the pressure difference persists from the conservative states.
No tighter inversion or chemistry change is adopted.

The native-activation plus constant-density-atol trajectory still has z-Linf
ratio0.993447 at25ns. An output-preserving trace of that same decision case
is running in `build/pelef-initial25-matched-source-20261026`; its purpose is
to localize the remaining discrepancy before choosing a correction.


### Matched chemistry trace retains the neighboring-cell error

`build/pelef-initial25-matched-source-20261026` completed in306.365s, exit0.
Both ordinary outputs are byte-identical to the combined chemistry decision
case. `source_split.py` and `compare_stencils.py` retain hashed-input analyses.
At coarse(1,5,17), first-predictor z H2 hydro discrepancy falls to
-9.92172e-13 kg/m3/s, and all captured z-stencil pressures agree within
2.91e-11Pa. At corrected stage2, however, z hydro discrepancy is
-7.39025e-8 and x is-3.69221e-8. The adjacent coarse(1,5,16) pressure
differs by-4.37168637e-7Pa, while coarse(1,5,17) itself differs by
only3.7835e-10Pa. Neighbor reaction/source capture is therefore the next
localization step, not an EOS tolerance adjustment.

At the first fine predictor, upper halo(1,8,33) transverse velocity difference
remains2.61896993e-4m/s versus canonical2.63517178e-4. This directly
confirms that removal of the coarse first-predictor chemistry discrepancy
does not remove the later interpolation/transport discrepancy. New native
and matched-candidate captures at coarse(1,5,16),(1,5,17),(2,5,17) are
running in `build/pelef-initial25-matched-neighbor-20261026`.


### Direct native reactor input supersedes first-iteration face inference

Matched-neighbor capture completed in308.312s, with both ordinary outputs
byte-identical to the combined decision case. Native source and subsequent
native reactor captures likewise reproduce the two pinned exports.
`build/pelec-initial25-neighbor-reactor-20261026/compare.py` uses actual
`native_reactor.csv` nonreacting source, initial species/temperature/internal
energy, and before/after states. At coarse(1,5,16), actual first-reaction H2
difference is-1.96813e-14kg/m3; matched replay with actual candidate input
reproduces the candidate within3.47e-18, while actual native input reproduces
native within7.81e-18. Across three cells and two reactions the native-input
replay residual is at most1.31e-17. The observed discrepancy is attributable
to input forcing at these captured cells, not the matched reaction solver.

Correction to prior inferred-input interpretation: first-iteration native
energy forcing is reconstructed from intermediate states after computeTemp
in `Advance.cpp:145-165`, not simply the average of face divergences.
At coarse(1,5,16), inferred energy forcing differs from the directly recorded
input by-5.03646e6 W/m3. Therefore first-iteration face-inferred reactor
replays are not valid native-input comparisons. Second-iteration sources use
the direct average; local measured face substitutions remain independent
evidence. `comparison.json` in the native reactor directory is the authority
for the direct-input comparison.

At coarse(1,5,16), first-predictor H2 hydro divergence discrepancy is
-1.53513806e-6 kg/m3/s in x, with y/z near2e-12. Native face/stencil capture
is running in `build/pelec-initial25-neighbor-face-20261026` to distinguish
neighbor primitive inputs from reconstruction. No production change adopted.

Native neighbor-face capture completed with identical pinned exports.
Its `replay.py` reproduces actual candidate all-component and native H2
face fluxes. At coarse(1,5,16), first-predictor x H2 divergence difference
-1.53513756e-6 kg/m3/s falls to-4.12115e-13 with pressure-only substitution.
The x stencil identifies coarse(2,5,16): pressure difference-1.81826181e-5Pa,
density difference1.39e-17kg/m3, maximum species-fraction difference5.55e-17.
Unlike the earlier canonical chemistry difference, this matched-case pressure
difference is not accompanied by a resolved species discrepancy. Conservative
state/EOS input capture is running in `build/pelef-initial25-eos-cell-20261026`.


### EOS-cell conserved state agrees; actual temperature guess remains

`build/pelef-initial25-eos-cell-20261026` native and matched-candidate
captures both reproduce their ordinary reference outputs. `compare_eos.py`
records the conservative-state comparison and EOS replays. At first predictor
coarse(2,5,16), candidate-minus-native energy is2.91038e-11J/m3, density
1.38778e-17kg/m3, species differences at most6.93889e-18kg/m3. Using either
the old or native stored temperature guess, both conservative states recover
pressure within1.02e-10Pa of actual native. The actual candidate stage pressure
remains-1.81826181e-5Pa different. Thus neither conservative-state mismatch
nor these two tested initial guesses explains the stage pressure difference.

An isolated application observer adds `stored_temperature` to the stage CSV
in `build/pelef-initial25-temperature-trace-20261026`. It links unchanged
matched-chemistry numerical libraries and the same isolated reactor override.
Source transformation, commands and binary hashes are recorded. Its ordinary
outputs must match before interpreting the newly captured guess. No production
EOS policy or tolerance was changed.


### Actual stored-temperature replay and isolated admission fix

Temperature observer completed in304.098s; ordinary outputs remain byte-identical
to matched chemistry. At coarse(2,5,16), stored T is2047.6163901598154K both
at stage0 and predictor stage1. Replaying candidate state with this actual
guess reproduces the observed pressure exactly, including the-1.81826181e-5Pa
difference. Even native conservative state with the candidate guess reproduces
that pressure discrepancy to1.46e-11Pa. Old/native guesses instead agree with
native pressure near roundoff. The native default1e-6K inversion stop retains
this sub-tolerance guess difference; no tolerance change is warranted from
this observation alone.

Serial `reset_held_eb_tree_family_3d` and its MPI counterpart currently
write the recovered temperature even in reset_level=-2 admission, despite
retaining valid conserved/held state. An isolated serial source copy in
`build/pelef-initial25-preserve-seed-20261026` preserves temperature when this
admission succeeds; failed admission still follows existing reset/repair.
The same matched-chemistry25ns original-gate decision test is being prepared.
This is not yet a production change or qualification result.

The preserve-seed isolated build completed and the25ns original-gate pilot
is running. A companion `build/pelef-initial25-preserve-seed-canonical-20261026`
links the same tree object against canonical chemistry libraries and the
original density-atol convention, to separate the temperature-history change
from chemistry sensitivity. Each runner includes provenance, both field times,
and cumulative budget checks. Calendar rollover does not create a new goal.


### 2026-10-27 temperature-seed admission implementation

Both isolated25ns pilots pass provenance, original initial/final field gates
and cumulative inventory gates. With canonical chemistry unchanged, z-momentum
Linf changes9.25979582e-8 to1.53152724e-8 (ratio0.165395); y changes
1.68433400e-8 to1.69834421e-8 (ratio1.008318). With matched chemistry,
z Linf changes9.19911719e-8 to4.32738067e-9 (ratio0.0470413), and y
1.14334686e-8 to6.07220828e-9 (ratio0.531091). Canonical chemistry is
retained; the experimental chemistry settings are not adopted. No global
all-component improvement, long-trajectory agreement or convergence is claimed.

The scoped implementation now preserves stored temperature on successful
reset_level=-2 admission in serial and MPI, retaining existing repair on
failure and reset behavior at advancing levels. A three-rank regression
checks initial temperature seed preservation through serial/MPI source-stage
callbacks, using the existing independent held-energy test fixture. Builds
and focused tests are running. Pre-change file snapshots and the task-only
patch are in `build/pelef-preserve-seed-integration-20261027`. Restart and
production trajectory verification remain pending.

Initial integration test pass completed compilation but reported1/3 tests
passed: existing `unit_mpi_eb_pinned_eos_source_3d_3` failed its accepted-forcing
ledger assertion, and the new test's first unsupported `seed` source-kind
combination failed an independent UEINT-source assertion before the intended
temperature check. The seed check now extends the existing
`eos_source corrector_material reactor` case via a fourth argument; no
tolerance or existing assertion changed. Corrected regression and related
reactor tests are running. A separate baseline build of the prechange serial
and MPI modules tests both the existing ledger failure and the new seed
assertion. These failures remain unresolved until those comparisons finish.

The expanded second run passed1/5; both existing reactor-corrector tests and
the initially extended seed case encounter a Godunov reactor floating-point
exception before the MOL seed assertion. The baseline audit remains running.
The new regression now exercises admission directly: its stage0 source
callback checks the unchanged stored temperature, deliberately refuses the
source, and verifies full gas/temperature/held-energy/ledger rollback in
serial and three-rank MPI. Existing coupled tests and their failures are
retained, not reclassified as passing. The focused admission regression and
standard held-reset ledger case are rebuilding in unit_driver3.py.

Baseline audit reproduces the existing pinned-EOS accepted-forcing ledger
failure with prechange serial and MPI modules. This specific failure predates
the temperature fix. The first narrow admission test omitted the mandatory
internal_energy_flux argument paired with internal_energy, causing entry
refusal before callback; both serial/MPI calls now supply that argument.
The unit_driver4 rerun is pending. The separate Godunov reactor FPE remains
unclassified until directly rerun against baseline; no claim is made that
all observed failures are historical.

Final focused regression run (`record5.json`, `unit-tests5.log`) passes2/2:
standard three-rank held-reset ledger and new temperature-seed admission plus
source-refusal rollback. The fourth attempt exposed a test collective-order
error: require() performs MPI_Allreduce, so it cannot be called inside
owner-only loops with an idle rank. The final test aggregates locally and
performs one collective check on every rank; production code did not change
during these test-fixture corrections. Baseline final relink is running to
verify that prechange code fails the intended seed assertion, and to classify
the existing Godunov reactor failures independently.

Baseline-final audit now confirms the new regression fails specifically at
`serial admission preserves stored temperature seed` on prechange modules.
Both existing Godunov reactor-corrector FPEs also reproduce on prechange
modules (`baseline/final-record.json`), separately from the pinned-EOS ledger
failure. None of these existing failures was loosened or removed.

Rebuilt canonical selected serial/MPI targets and eight diagnostic cases pass
(`selected-record.json`, `trace-checks/.../checks.json`); these include expected
rejections and observer-on/off output checks, not eight physical validation
cases. Fresh n8 serial200ns versus MPI2-prefix/MPI3-restart comparison is
running, as are43 focused transport/MOL/callback tests. Fresh n32 dt25ns
trajectory to200ns runs in `build/pelef-preserve-seed-trajectory200-20261027`
with unchanged canonical chemistry and all nine original field-time gates.
Earlier baseline fails175/200ns gates; the new outcome is pending.


### Production200ns trajectory after seed preservation

`build/pelef-preserve-seed-trajectory200-20261027` application completes
in1851.133s with exit0 and cumulative inventory gates passed. Checkpoint hash
was independently verified, and the25ns ordinary output is byte-identical
to the isolated canonical preserve-seed pilot. Original field gates pass
0..150ns, fail rho_E/rho_e at175ns, and density/zmom/rho_E/rho_e at200ns.
The scientific gate status is failed despite successful execution.

Saved-field comparison shows z-Linf ratios versus the prechange trajectory
0.1654/0.1447/0.1414 at25/50/75ns, then0.6541/1.0401/0.9823/1.1001/1.1009
at100/125/150/175/200ns. At200ns, z Linf is3.60735385e-6 versus
3.27660658e-6, both atfine(26,14,32); y is8.61447377e-7 versus
7.82466890e-7. The initial seed correction does not resolve the later
interface amplification. No long-trajectory or C2 qualification is claimed.

Fresh n8 MPI2-to-MPI3 restart reaches200ns with prefix/final ordinary CSVs
byte-identical to continuous serial and the source checkpoint unchanged.
The43-case broad suite remains running; recursive non-MOL two-/three-rank
cases have timed out at300s and are not counted as passes or numerical
mismatches. A separate matched-chemistry+preserve-seed200ns sensitivity is
running in `build/pelef-preserve-seed-matched200-20261027`, recording the
late worst-cell and coarse donor stages. Its chemistry settings remain
experimental and separate from production.

Broad43-test run completed in3590.54s:40 passed, one previously reproduced
spatial-halo retained-energy assertion failed, and recursive non-MOL MPI2/3
register cases timed out at300.12/300.10s. Previous same-suite runs completed
those two cases in295.36/288.37s. This supports investigating timing margin,
not classifying the new timeouts as numerical mismatches. They remain
unpassed in this run. `retry_timeouts.py` is prepared to repeat only these
two with unchanged300s limits after the concurrent trajectory finishes.


### Matched chemistry plus preserved seed: completed200ns sensitivity

`build/pelef-preserve-seed-matched200-20261027/record.json` records exit0,
1868.409s application time and passed inventory gates. Input, executable and
checkpoint SHA256 were independently verified. Original field gates pass
0..150ns;175ns fails zmom/rho_E/rho_e;200ns fails density/ymom/zmom/rho_E/rho_e.
At200ns z Linf is4.986983487e-6,1.5220 times the pre-seed canonical trajectory
and1.38245 times the canonical preserve-seed trajectory. This remains an
experimental sensitivity; no chemistry convention or tolerance is adopted.

`analyze_late.py` matches native face events to candidate depth1 probes4/5
at fine(26,14,32)/(25,14,32), using the actual second-corrector forcing.
All100..200ns y/z endpoint increments close within1.977e-15. At the first
cell z-error growth is-4.979609458e-6; hydro contributes+8.428831284e-7
and the total-minus-hydro transport residual is-5.822492586e-6.
This is a discrete balance attribution, not proof of an interpolation bug.
`compare_late_stencil.py` also matches the reconstructed pressure and all-axis
stencils. At187.5ns predictor, halo(26,14,33) has transverse v difference
+0.00882862724m/s. `replay_late.py` validates528 local Riemann cases against
captured candidate fluxes and native normal fluxes. At the largest upper
z-face discrepancy (event44), pressure substitution reduces the hydro
flux residual from-0.00140076Pa to-0.000908432Pa, and all halo primitives
to-0.000767305Pa; interior differences remain. This replay does not isolate
transport's transverse gradients or establish a trajectory correction.

The two timed-out recursive non-MOL tests are rerunning after the trajectory
finished, with their original300s limits. Notification records are under
`build/pelef-preserve-seed-integration-20261027/timeout-retry-notification`.
Further numerical edits await the test result and exact donor-stage audit;
C2 remains open.


The unchanged-limit retry completed: both recursive non-MOL MPI2/3 tests
again timed out (300.11/300.09s). Both reached layout3/step1 with displayed
composite residual3.55e-16; final completion and remaining rollback checks
were not reached. No pass or numerical-regression classification follows.
The notification wrapper exits0 because the diagnostic driver stores the
CTest result without propagating it; `timeout-retry.json` returncode8 and
`status=failed` are authoritative. This retry eliminates simultaneous long
trajectory execution as a sufficient explanation, but not all host-load or
runtime effects. Original300s limits remain unchanged.

`build/pelef-preserve-seed-donors150-20261027` now runs a fresh canonical
preserve-seed trajectory to150ns with isolated donor logging. It captures
center/six neighbors, temporal fraction and ancestor endpoint states for
fine(25/26,14,33), plus stage/reaction clock markers. Its acceptance is
byte identity of all seven ordinary outputs and exact equality of the
four diagnostic CSV streams through step6 against the already completed
canonical trajectory. `analyze.py` is prepared to replay the common limiter
and retain clock bracketing without guessing physical stages from nearest
values. Build/run source identity and commands are in record.json. No
production numerical source was changed for this observation run.


The canonical preserve-seed donor run completed in1408.097s application
time. All seven ordinary outputs and four diagnostic streams match the
previous canonical trajectory; binary/checkpoint/output hashes were
independently rechecked.144 interpolation captures replay, and54 stage
clock markers match the observer CSV sequence exactly. Raw clocks inside
the solver are root-step-relative; `timeline.py` reconstructs absolute
clocks from coarse stage0 boundaries and checks them against the CSV.

For fine(26,14,33) at the137.5ns predictor requesting150ns, the common
x-limiter is0.491332844249, controlled by mz, versus native0.495740819761.
The captured halo velocity differences are[-4.79174778e-4,+1.29535380e-3,
+9.30715769e-9]m/s. Replacing only the seven mz donor entries with native
values gives[4.7191e-10,-8.3597e-11,1.4292e-15]m/s. Replacing species does
not remove the discrepancy. Both interpolation algorithms reproduce their
captured states; this does not establish an incorrect limiter formula.
At coarse(12,7,17), mz error changes from+2.405056946e-10 before125ns
synchronization to-6.934769679e-9 after it, a-7.175275374e-9 jump.
The same causal amplification route persists after the initial seed fix.

A scoped canonical25ns checkpoint continuation to125ns is running in
`build/pelef-preserve-seed-sync125-20261027`, capturing this coarse face
and the four neighboring fine faces. Its gate is identity of all five
25..125ns ordinary outputs and an unchanged source checkpoint. The prepared
analysis compares accepted transfer integrals with cached native faces,
then separates hydro and transport contributions. No new numerical change
or physical qualification is claimed.

The first sync125 launch refused simultaneous `initial_state_file` and
`restart_file` before stepping (0.114s); its log and failed record are
preserved under application/. `driver2.py` removes the initial-state input
and runs under application2/, with notification2/. The source checkpoint
was unchanged. This setup refusal is not a numerical trajectory failure.


The corrected sync125 continuation completed in874.131s. All five ordinary
25..125ns outputs match the canonical preserve-seed trajectory, and the
source checkpoint/hash and four trace hashes were verified. Accepted
coarse/fine flux integrals reproduce the donor synchronization jump:
candidate correction-1.98264277079e-5 versus native-1.98192524323e-5;
difference-7.17527565142e-9 versus observed-7.17527537386e-9 (closure2.78e-16).
This establishes the discrepancy in the incoming fluxes rather than the
synchronization arithmetic at this face.

Fine(23,13,32) and(23,14,32) contribute integrated transport residuals
-4.78969579e-12 and-4.41187314e-12; their hydro differences are
-6.62005368e-13 and+2.81674784e-13. These are face flux-time integrals,
not direct cell momentum errors. The coarse total-flux integral difference
is8.36518e-14. Analysis is under `build/pelef-preserve-seed-sync125-20261027`.

The next bounded observation pairs
`build/pelef-preserve-seed-transport125-20261027` (canonical25..125ns restart)
with `build/pelec-sync-transport125-20261027` (native0..125ns). It captures
transport states around fine(23,13,32), all face fluxes and transverse
neighbors. Gates require five candidate ordinary outputs and six native
exports to reproduce existing results before local gradient substitution.
No production numerical sources or tolerances change.

Native sync-transport125 observation completed in109.076s; all six0..125ns
exports reproduce the fixed reference. Binary, flux and export hashes
were independently verified. `replay_native.py` evaluates the captured
fine(23,13,32) states through the production transport helper: momentum
face-flux residual maxima are6.41e-15/9.60e-15/7.16e-15Pa; energy flux
residual is1.00e-9W/m2. This verifies the native-side local replay, not the
pending candidate comparison or trajectory qualification.


Candidate sync-transport125 observation completed in867.449s; all five
ordinary outputs reproduce the canonical run, with unchanged restart and
verified trace hashes. Local replay reproduces both momentum transport
fluxes within1e-14Pa. At fine(23,13,32) upper-z face, the100..125ns
integrated transport difference is-4.78969586e-12. Native transverse
velocity gradients reduce it to-3.90681412e-14 (99.18% removal); normal
gradients alone leave-4.75061917e-12 and temperature alone-4.78970162e-12.

The largest instantaneous difference occurs at100ns old stage. Saved
coarse10/11/12,7,17 states reproduce the actual halo(22,13,33) density and
velocity exactly for each solver using the x-limited interpolation.
The common x-limiter is controlled by mz: candidate0.003207258419 versus
native raw-CGS0.004432966188. At coarse11,7,17, mz error is already
4.12640024e-10 at25ns, growing to4.95838589e-9 at100ns. Matched-chemistry
sensitivity does not remove that earliest error (4.93189126e-10 at25ns).
Artifacts are upper-summary.json, saved-x-limiter.json and
coarse-mz-history.json under the candidate observation directory.

The next short pair `build/pelef-coarse11-initial25-20261027` and
`build/pelec-coarse11-initial25-20261027` captures the first25ns coarse
update and four fine faces feeding its lower-z synchronization. Gates
require both0/25ns ordinary outputs to match existing references.
Prepared analyze.py reconstructs native pre-sync momentum from full
face divergence and verifies final momentum after reflux, then separates
candidate own-update and synchronization differences. No numerical
correction is inferred solely from the downstream amplification.


Both first25ns coarse11 observations reproduce their0/25ns outputs;
candidate application302.304s and native29.981s. Full-face native momentum
balance closes to1.33e-15. Candidate own coarse update differs by only
1.03295e-13, whereas synchronization differs by4.12538059e-10, explaining
the final4.12640024e-10 error. The four fine-face hydro integral differences
are[-2.70993e-13,+7.03522e-13,-2.72758e-13,+4.35698e-13]; corresponding
transport residuals are[2.91251e-15,-3.26706e-14,8.82646e-14,-1.16729e-13].
Thus the earliest injection is predominantly hydro, distinct from the
later transverse-transport amplification. `comparison.json` and
`flux-split.json` live under `build/pelef-coarse11-initial25-20261027`.

`build/pelec-coarse11-pressure25-20261027` adds native reconstructed-face
and stencil capture at fine(22,13,32), the largest positive hydro contributor.
Existing candidate probe3 provides the matching data; no candidate rerun
is needed. Prepared replay.py and stencil_compare.py require native output
identity and check local flux reproduction before pressure substitutions.


Native pressure25 capture reproduced both exports (29.724s). Local Riemann
replay reproduces captured fluxes. At fine(22,13,32) upper-z first-corrected
stage, the3.74063384e-5Pa hydro flux difference becomes-1.45519e-11Pa when
only halo pressures are replaced. At the second corrected stage3.77620e-5Pa
becomes zero; replacing velocity/species does not remove the discrepancy.
Artifacts: `build/pelec-coarse11-pressure25-20261027/replay-comparison.json`
and stencil-comparison.json. Its notification directory was accidentally
created outside the repository build tree at
`/home/gtt13/projects/pelec-coarse11-pressure25-20261027/notification`;
scientific outputs remain in the repository build directory.

The next short thermo capture reuses the already built candidate transport
observer for a fresh25ns run in `build/pelef-coarse11-thermo25-20261027`.
A native pair in `build/pelec-coarse11-thermo25-20261027` directly records
conserved S (including stored T and UEINT) before primitive conversion.
This avoids inferring conserved energy from pressure. Prepared compare_eos.py
will test original/converged inverse EOS and exchanged temperature guesses;
converged inversion is diagnostic only, with no tolerance change adopted.


Thermo25 observations reproduce both ordinary outputs (candidate304.129s,
native29.904s). Direct conservative-state EOS replay reproduces candidate
halo pressure at fine(22,13,33) exactly. At first predictor/corrector,
its native pressure difference3.74063384e-5Pa falls to2.93046469e-7Pa when
using the native stored temperature as the inversion guess, or when the
inversion is converged diagnostically. Native captured state reproduces
native pressure within9e-11Pa with either guess. Candidate stored guess
is1976.462840736262K, native1976.462816961302K; candidate converged temperature
is1976.462840012323K. Thus most of the pressure discrepancy is the
inverse-EOS stopping path, while a smaller state/composition difference
remains. See `build/pelef-coarse11-thermo25-20261027/eos-comparison.json`.

Source inspection finds that coarse-fine halo filling validates the
interpolated state through EOS, then publishes the recovered temperature
for another primitive recovery. An isolated pilot at
`build/pelef-halo-seed-pilot25-20261027` retains the interpolated guess after
successful validation for coarse-fine, non-inlet, non-paired-pressure halos.
It uses canonical chemistry, unchanged tolerances and the original0/25ns
field/inventory gates. Production source is untouched; a successful pilot
would still require an explicitly scoped serial/MPI policy, rollback tests,
restart checks and longer-trajectory qualification. Generic halo default
semantics must not be silently changed during integration.


The isolated halo-seed25 pilot passes both original field gates and the
inventory gate. Target coarse11 mz error falls from4.12640024e-10 to
-3.78149268e-11 (90.8% reduction in magnitude); the first predictor halo
pressure difference falls from3.74063384e-5 to2.93031917e-7Pa. Global25ns
momentum Linf is essentially unchanged, so the short pilot is not a
trajectory solution. A fresh200ns pilot is running with all nine gates.

A scoped production implementation adds default-false
`preserve_temperature_seed` to serial/MPI halo filling, preserving successful
validation and rollback. Ordinary nonrotating MOL hydro and total-energy
MOL transport request it; experimental momentum/pressure policies retain
existing behavior. MPI checks rank agreement for the policy. Source
snapshots and task-only.patch are under
`build/pelef-halo-seed-integration-20261027`. The MPI halo regression adds
serial/MPI donor-seed checks, default recovery behavior, inconsistent-policy
refusal and invalid-state rollback on both owner layouts. Tests/builds are
pending; this is implemented but not yet verified. The first unit launcher
failed with undefined Python metadata before building; corrected launcher
runs under unit-notification2. Selected application diagnostics and n8
serial/MPI2-to-MPI3 restart checks run independently. integrated25.py is
prepared to compare the scoped production path to the isolated pilot.


The first compiled seed-policy regression passed the existing serial halo
test but failed both MPI test programs at their new serial narrow-halo
seed assertion. This exposed the unconnected width1 branch; the wide MOL
pilot branch was already covered. The explicitly requested policy now
also preserves validated nonphysical coarse donors in width1, with the
default unchanged. `unit-record.json` retains the failing evidence;
unit_driver2.py/unit-notification3 rerun the same three tests. Selected
application checks started before this narrow-branch addition, so a final
selected_check2.py rebuild/check/restart pass is prepared after those
checks finish. No earlier binary is relabeled as the final source.


After connecting the explicitly requested policy in the narrow serial
branch, unit-record2.json reports3/3 passed: existing serial halo0.25s,
MPI2 halo8.28s and MPI3 halo7.85s. Both MPI programs exercise two owner
layouts, including all patches on rank0, seed retention/default recovery,
rank-policy disagreement and invalid-state rollback. The prior failing
assertion is retained as evidence that the new check detects an unconnected
seed path. Fourteen relevant serial/MPI MOL, transport and reflux cases
are now rebuilding/running in mol-notification; the previously reproduced
serial spatial-halo retained-energy failure remains separately classified.


The first scoped application build and eight diagnostic cases passed.
The n8 serial200ns versus MPI2-prefix/MPI3-restart check passes with
byte-identical prefix/final CSVs and unchanged checkpoint (66.577s combined
checker time). This build started before the later width1 policy addition;
selected_check2.py now rebuilds and verifies the final source separately.
No n32/200ns reference qualification follows from the n8 restart result.


The related14-test MOL/transport/reflux suite completed in371.32s:
13 passed; unit_eb_tree_mol_spatial_halo_3d again failed at the retained
energy output preserves tree gas rate assertion. This exact assertion
was independently reproduced before the halo-seed change in
`build/pelef-transport-total-integration-20261026/baseline/test.log`.
The suite is recorded as failed, not14/14. Inspection shows its two calls
also differ in held-state/coarse-halo preparation, so resolving that API
invariance remains separate work; neither its assertion nor tolerances
were loosened here. Final-source selected/restart and pilot200 remain pending.


Final-source selected-record2.json passes the rebuilt serial/MPI application
checks and eight diagnostic cases. restart-check2/run.json confirms exact
n8 serial200ns versus MPI2-prefix/MPI3-resume CSVs and unchanged checkpoint;
executable/checkpoint hashes were independently verified. integrated25.py
is now running the n32 production binary against the isolated pilot's two
ordinary outputs. This equivalence check and the200ns pilot's original
field gates remain outstanding; no C2 qualification is claimed.


The final production n32 integrated25 check completes in305.423s. Both
0/25ns ordinary outputs and all four stage/sync/reaction/stencil CSVs are
byte-identical to the isolated halo-seed pilot; the production executable
hash is380ba789e9fa02825d2864f4585f0abb9229f6f5d9a8cb0caa1ebc258d1324a6.
Evidence: integrated25-equivalence.json under the integration directory.
This establishes scoped implementation equivalence for25ns, not equivalence
for the still-running200ns pilot or complete fixed-PeleC qualification.


### Halo-seed200ns pilot: momentum gates improve, energy remains outside

The fresh isolated200ns pilot completes in1836.396s with inventory gates
passed and checkpoint hash independently verified. All original field
gates pass0..175ns. At200ns only rho_E/rho_e fail; density and all three
momentum gates pass. Versus the current canonical admission-seed baseline,
z Linf falls3.607353855e-6 ->8.055036132e-7 (ratio0.223295), while energy
Linf falls0.01613097746 ->0.00285340744 (ratio0.176890). The remaining
scaled energy error1.197653727e-8 exceeds the unchanged1e-8 limit. The
z-maximum moves to fine(9,10,32); energy maximum remains fine(26,14,32).
This is improved reference agreement, not C2 completion or physical validation.

`build/pelef-halo-seed-production200-20261027` now continues the verified
production25ns checkpoint to200ns. That checkpoint is byte-identical to
the isolated25ns pilot. Its step0 CSV is explicitly copied from the same
verified production prefix, with path/hash recorded, for the full inventory
comparison; it is not reported as a fresh0ns run. The continuation checks
eight ordinary outputs against the isolated200ns pilot, final budget CSV,
all nine field tables and the original inventory checker. It records the
new momentum maximum and late energy maximum plus nearby coarse donors.
Prepared analyze_late.py includes energy forcing closure as a diagnostic;
a residual is retained rather than automatically attributed to face flux.


### Production halo-seed continuation and remaining energy balance (2026-10-27)

`build/pelef-halo-seed-production200-20261027/verified.json` independently
rechecks the production executable and source checkpoint hashes, all eight
25..200ns CSVs against the isolated pilot, and the final budget CSV.
The continuation completed in1510.576s. Original field gates pass0..175ns;
only rho_E/rho_e fail at200ns (scaled1.1976537e-8 versus unchanged1e-8).
The cumulative inventory checker passes. This confirms the production
integration, not C2 convergence or physical validation.

`late-flux-comparison.json` matches actual second-corrector forcing with
cached native face divergence. At fine(26,14,32),100..200ns total-energy
error growth is-0.002848078205716 J/m3; integrated forcing difference is
-0.002848078206804 J/m3 (residual1.09e-12 J/m3). Integrated hydro differences
are +0.000413340192362,+0.000076773342118,+0.000351061781811 J/m3 in x/y/z;
the total-minus-hydro transport difference is-0.003689253523095 J/m3.
This identifies a spatial transport contribution, not an exact transport
operator defect. No tolerance, chemistry convention or limiter was changed.

`build/pelef-halo-seed-transport150-20261027` captures the actual transport
halo states and six face fluxes around that cell from the verified25ns
checkpoint through150ns. Its ordinary outputs must reproduce production
before replay results are interpreted. Cached native0..150ns transport
captures provide the matched reference. The prepared replay substitutes
temperature, composition, density, velocity and gradients separately;
source replay residuals remain explicit rather than assumed zero.


The transport150 capture completed in1067.180s and all six25..150ns ordinary
CSVs independently match production (`verified.json`). Local transport replay
reproduces actual energy fluxes with maximum residual1.82e-12 W/m2 for the
candidate and1.06e-9 W/m2 for native. `summary.json` integrates old/corrected
face differences over25..150ns, using6.25ns weights and0.01/64m cell width.
At fine(26,14,32), the transport energy-error increment is-1.6532749e-4 J/m3:
z contributes-1.6626832e-4, with small opposing x/y contributions. Mode splits
are-1.0206336e-4 conduction,-6.3290236e-5 species/barodiffusion and+2.61027e-8
viscous work. Substituting native face temperatures reduces the total to
-6.3239288e-5; substituting species reduces it to-1.0150206e-4. Gradient
substitution barely changes it. These are local replay sensitivities, not
independent additive causal corrections or a validated implementation fix.
The worst instantaneous upper-z face difference is3.920071 W/m2 at native
event35 (candidate event29, predictor requesting150ns); its upper halo
T differs by-0.001183112K. A combined donor/transport capture is running in
`build/pelef-halo-seed-donors150-20261027` to identify the upstream interpolation
inputs. It retains canonical chemistry and the unchanged acceptance gates.


The combined donor capture completed in1108.799s; all six ordinary outputs
and both transport traces independently match the preceding production and
transport-only runs (`pelef-halo-seed-donors150-20261027/verified.json`). All120
captured interpolations replay successfully, with observer-verified stage
clocks including the25ns restart offset. Native event35 matches the predictor
requesting150ns. At halo(26,14,33), the x shared limiter is0.4965109695233534
versus native0.4957408197605918, controlled by mz. Replacing only the seven
coarse mz donor values reduces the actual temperature discrepancy from
-0.001183111687K to+5.12998e-8K, and H2 mass-fraction discrepancy from
-4.72781e-8 to+5.56818e-12. Replacing donor temperature, energy or species
alone does not remove it (`origin.json`, `thermo-origin.json`). This supports
amplification of small input differences; it does not demonstrate an incorrect
interpolation formula. The coarse(12,7,17) mz error jumps from1.39069005e-10
before125ns synchronization to1.01775530e-9 afterward. The next capture,
`pelef-halo-seed-transport125-20261027`, records the four contributing fine
faces and transport states for a same-stage comparison to cached native data.
C2 remains open; no additional production numerical change was made.


The halo-seed transport125 run completed and its five ordinary CSVs and trace
hashes were independently verified. `sync-comparison.json` reproduces the
125ns coarse mz discrepancy jump8.7868629e-10 from incoming transfers
8.7868601e-10 (residual2.78e-16). Fine(23,13,32) and(23,14,32) dominate,
with integrated transport differences5.93037e-13 and6.65726e-13 Pa s;
hydro is small. Local transport replay at(23,13,32) reproduces both momentum
fluxes within1e-14Pa. Replacing transverse gradients changes its integrated
upper-z difference from5.93037e-13 to-6.40725e-15 Pa s (98.9% reduction).
Normal-gradient or temperature substitution does not remove it. This closes
the remaining amplification route through incoming transport, not a defect
in reflux arithmetic. Results are in`pelef-halo-seed-transport125-20261027`.

A short isolated25ns decision test in`pelef-halo-seed-matched25-20261027`
combines current production objects (including both temperature-seed fixes)
with the previously verified native activation-basis mechanism archive and
native density-atol reactor object. It records their exact hashes and links
only those two sensitivity objects; production chemistry remains unchanged.
The decision compares the upstream coarse mz errors and original field and
inventory gates before considering any longer chemistry sensitivity run.


The current-halo-seed matched-chemistry25ns decision test completed, passed
both original field tables and inventory gates, with binary/input-object
hashes independently verified. It does not justify a long sensitivity run:
coarse(11,7,17) mz error changes from-3.78149e-11 to+4.29090e-11 (absolute
ratio1.13471), and coarse(12,7,17) from+2.37987e-11 to-2.68126e-11
(ratio1.12664). Coarse13 improves, and fine(22,13,32) improves27.5%, so
improvement is not uniform. Production chemistry remains unchanged.
`cell-comparison.json` and both canonical/matched `pressure-replay.json`
retain the decision evidence. Riemann replays reproduce actual fluxes;
replacing halo pressure removes corrected upper-face discrepancies, while
interior-face pressure differences remain. The next short canonical run,
`pelef-halo-seed-thermo25-20261027`, reuses the verified transport observer
binary to capture halo conserved states and stored temperatures against
cached native conserved-state data. It must reproduce both ordinary CSVs
before EOS-path interpretation. No new numerical fix was adopted here.


The fresh canonical halo-seed thermo25 capture passes both original field
and inventory gates and independently reproduces both production CSVs.
`eos-comparison.json` shows that the first predictor/corrector halo pressure
residual2.93032e-7Pa persists with native temperature guesses and converged
inverse EOS. Native conserved-state replay agrees within5.83e-11Pa. Thus
this remaining residual comes from conserved-state differences, not the
previously corrected seed publication. `state-decomposition.json` retains
nonadditive EOS input substitutions: replacing energy alone gives+1.46760e-4Pa,
and density/species together gives-1.46739e-4Pa. Their near cancellation
requires tracing the complete interpolated state. Density-only substitution
was rejected by species consistency and is excluded, not interpreted.
The paired first25ns donor capture in`pelef-halo-seed-donors25-20261027`
and`pelec-coarse11-donors25-20261027` records the actual coarse11,7,17
interpolation stencil and timed halo22,13,33 input. Both ordinary-output
identities are required before assigning the early limiter cause.


The paired initial25ns donor capture completed and reproduces both solvers'
ordinary tables. Analysis replayed24 interpolations and verified all9 stage
markers against observer clocks. The first predictor halo22,13,33 shared
x-limiter is0.03643392639 versus native0.03642882572, controlled by mz.
Replacing donor mz reduces the pressure discrepancy2.93105e-7Pa to2.19836e-7;
replacing species reduces it to7.32834e-8. Neither alone removes the early
residual (`pelef-halo-seed-donors25-20261027/origin.json`).

The existing thermo capture also exposes an interior-cell seed dependence:
at fine22,13,32 first predictor, actual pressure error2.02141e-7Pa becomes
-3.17232e-9Pa using the native stored temperature. A converged inversion does
not reproduce native's finite-stop path (`interior-eos-comparison.json`).
Source inspection finds reactor endpoint recovery followed by recovery in
`react_sdc_state_3d`, then the EB tree's post-reset recovery. Native
ReactorUtils unflatten recovers temperature, followed by Advance computeTemp.
`pelef-reactor-seed-pilot25-20261027` tests preserving the reactor-exit seed
at the intermediate publication while retaining EOS validation. It changes
only an isolated object, with canonical chemistry; production is unchanged.
Its result must determine whether an explicitly scoped serial/MPI fix is
justified. Generic, rotating and other application behavior is not qualified
by this pilot.


The reactor-exit seed25 pilot passes original field/inventory gates, but both
ordinary CSVs are byte-identical to production; target momentum and pressure
are unchanged. This intermediate-publication change is not adopted.

A separate code-level discrepancy is now identified: pinned
`ReactorUtils.H::fKernelSpec` recovers temperature from
(rhoe_init+t*rhoe_source)/rho at every RHS call, using the temperature
coordinate as the inverse-EOS seed. Current `full_cvode_rhs` instead directly
uses that coordinate for rates and thermodynamics. The isolated
`pelef-full-rhs-energy25-20261027` pilot routes RHS evaluation through the
existing signed-trial forced reactor recovery, keeping canonical chemistry,
endpoint admission and original gates. Its existing approximate full Jacobian
is retained only for the decision test; derivative semantics and signed-trial
coverage require explicit verification before production adoption. This is a
new RHS-policy hypothesis, not evidence that the prior seed fix failed.

The first RHS pilot build stopped before execution because its isolated compile
command omitted the pinned SUNDIALS_INT32_T definition. The definition was
verified in the production Ninja build and restored; notification2 owns the
retry. The failed build is retained and is not counted as a solver result.


### Full-coordinate RHS energy constraint implementation (2026-10-27)

The affine-energy RHS25ns pilot completed and passes both field tables and
inventory gates. It does not improve the localized momentum discrepancy:
coarse11/12 absolute error ratios are1.00694/1.00083 and fine22,13,32 is1.41180.
This is kept separate from the independently identified RHS-policy mismatch.

`src/chemistry/sundials_forced_reactor_mod.F90` now uses the existing
signed-trial energy recovery for each full-coordinate RHS evaluation, exposing
`evaluate_forced_reactor_full_rhs` for direct physical-kernel verification.
Rates and heat capacity use the recovered temperature; the ODE temperature
coordinate supplies its inverse-EOS seed. The reduced-coordinate backend and
endpoint rollback/admission stay unchanged. The old independent-temperature
hybrid Jacobian is removed for this backend; CVODE's dense numerical
Jacobian differentiates the actual constrained RHS. No mechanism or gate
tolerance was changed.

The full and reduced forced-reactor tests pass2/2. The added independent
oracle uses equal-thermodynamics A->B with k(T)=0.01*T, prescribed energy at
1000K and displaced temperature coordinates: the rate remains8 in mass-density
units and the thermal source is zero. Existing signed-intermediate,
nonnegative-endpoint, mass/energy and rollback cases also pass. An isolated
negative control reinstates the prior independent-temperature RHS policy to
check that the new oracle detects it; it fails at the energy-temperature assertion as expected. Exact task edits are isolated in
`pelef-full-rhs-integration-20261027/task-only.patch` against pre-edit backups.
The notification pipeline rebuilds serial/MPI selected applications, runs
observer diagnostics and rank-changing restart, then executes a fresh canonical
n32/dt25ns200ns reference comparison. These are pending, not passed; C2 remains
open regardless of the focused kernel result.


The numerical-Jacobian full-RHS integration pipeline completed: unit2/2,
8 application diagnostics and n8 MPI2->MPI3 restart pass. Restart CSV/checkpoint
identities were independently rechecked. The fresh n32 trajectory completes
with inventory closure but FAILS field gates from25ns onward. At25ns density
scaled error is2.04193e-6 (limit1e-8), x-momentum maximum3.07122e-4
(limit1e-6), and energy scaled error2.85718e-6 (limit1e-8). This is a failed
reference result, not qualification. `field-summary.json` retains all times.

Pinned `ReactorCvodeJacobian.cpp` uses the ODE temperature coordinate in
RTY2JAC even though fKernelSpec recovers temperature from affine energy.
The production full backend now retains this independent-temperature
approximate Jacobian policy: analytic density columns and a temperature
column evaluated with a separate independent-temperature kernel. The actual
RHS continues to use the tested energy constraint. The numerical-Jacobian
source is preserved as`numerical-jacobian-source.F90`; task-only final edits
are in`pelef-full-rhs-nativejac-integration-20261027/task-only.patch`.
A short25ns pipeline rebuilds/retests unit, application and restart behavior
before further trajectory work. Those latest-policy checks are pending.


The independent-temperature approximate-Jacobian implementation passes the
final-source2/2 reactor tests,8 application diagnostics and rank-changing
restart. Source, binary and checkpoint hashes plus prefix/final restart CSV
identities were independently rechecked (`nativejac-integration/verified.json`
under the dated build directory). The original0/25ns field and inventory gates
pass. Local coarse11/12 mz error magnitudes rise3.08/3.29%, while fine22,13,32
falls14.4%; this is not uniformly improved accuracy. The same-binary25ns
checkpoint now continues to200ns in`pelef-full-rhs-nativejac200-20261027`.
Step0 is explicitly copied from its verified fresh prefix for full inventory
accounting. No prior sensitivity output is substituted. The comparison script
uses the original nine gates and pre-RHS halo-seed trajectory as baseline.
Notification metadata is located outside the checkout at
`/home/gtt13/projects/pelef-full-rhs-nativejac200-20261027/notification`;
scientific output remains under the repository build directory.


The nativejac200 continuation completed with inventory closure and verified
binary/input/source-checkpoint/final-checkpoint hashes. All original field
gates pass0..175ns; only rho_E/rho_e fail200ns. Energy scaled error is
1.21649977e-8 versus1e-8, absolute0.00289830810 J/m3 at fine26,14,32.
z-momentum maximum is9.22008218e-7 at fine9,10,32, below1e-6. Relative to
pre-RHS halo-seed, energy grows1.57% and z maximum14.46%, so no trajectory
accuracy improvement is claimed. `baseline-comparison.json` records all nine
times. Current late-stage forcing closes the100..200ns energy error growth
at fine26,14,32 within3.97e-11 J/m3 (`late-flux-comparison.json`).

`pelef-full-rhs-matched200-20261027` now tests the already-generated native
activation-basis archive and native density atol with the CURRENT energy-
constrained RHS and independent-temperature approximate Jacobian. The older
chemistry comparison used the old RHS and is not a substitute for this
combination. This is an isolated sensitivity, not tolerance relaxation or
production mechanism adoption. Original field/inventory gates stay fixed;
its baseline is the current canonical nativejac200 trajectory.


The current-RHS native-chemistry200ns sensitivity completed and its binary,
mechanism archive, input and checkpoint hashes were independently verified.
Inventory closes and all field gates pass through175ns, but rho_E/rho_e
still fail200ns. Energy scaled error1.37673454e-8 (limit1e-8) is13.17% larger
than current canonical; y maximum grows13.21%, z maximum falls1.68%.
No production chemistry/tolerance change is adopted. See
`pelef-full-rhs-matched200-20261027/baseline-comparison.json`.

`pelef-full-rhs-mesh-20261027` refreshes the8/16/32 grid comparison using the
same CURRENT executable, unchanged initial profiles, dt100ns and100ns endpoint.
It checks original per-grid reference fields and inventories, then computes
conservative adjacent-grid restrictions for both candidate and saved native
solutions. Endpoint norms and evolution-increment norms are reported separately
to expose initial-profile sampling. This is not a time-refined spatial-order
claim, and previous-source mesh results do not qualify the current RHS.
The new mesh calculation is pending.


The current-source mesh100ns refresh passes all six initial/final field tables
and all three inventory checks, with same-binary and analysis-input hashes
independently verified. Runtime8/16/32 is22.30/155.29/1171.08s. Conservative
16->32 versus8->16 endpoint differences give z-momentum ratios Linf1.909647,
L1=0.400266,L2=0.461234. Native gives1.909647,0.400266,0.461234 respectively.
Subtracting initial sampling leaves candidate y-Linf ratio0.948318 versus
native0.948318; average errors decrease. Thus weak maximum convergence is
shared with the fixed reference, not cured by the RHS correction.

`pelef-full-rhs-mesh25-20261027` repeats the comparison at dt25ns/100ns.
Only8/16 candidate/native cases are new;32 reuses verified current-source
0/100ns tables from the200ns continuation with explicit copy hashes. No new
32-grid100ns inventory result is claimed from that reuse. Endpoint and
initial-subtracted evolution norms will distinguish time-step contamination
from the existing interface-dominated spatial trend.


The dt25ns mesh comparison completed: new8/16 native/candidate runs pass
all field and inventory gates;32 uses hash-verified same-binary0/100ns tables.
The z-momentum adjacent-grid Linf ratio remains1.728810 (native1.728810),
while L1/L2 ratios fall to0.454989/0.474147. Initial-subtracted y-Linf ratio
is0.948968 in both. Thus the maximum issue persists at one-quarter the step;
this is not spatial maximum-norm qualification. All analysis input hashes
and32 reuse hashes were independently verified.

`pelef-current-native-reactor-replay-20261027` directly replays six captured
native reactor inputs (three coarse cells, two stages), including actual
rhoY,T,rhoe and external sources, with current libraries. Canonical settings
give max species difference4.48e-13 kg/m3 and T difference1.45e-8K. Native
activation basis plus native density-atol reduces these to7.08e-16 kg/m3
and4.55e-12K. These bounded same-input results support upstream forcing and
interpolation amplification; they do not qualify all reacting trajectories.

`pelef-full-rhs-time12-200-20261027` runs matched original native and current
canonical32-grid cases at dt12.5ns through200ns, retaining17 field times and
original inventory gates. The new comparison will test whether the remaining
dt25ns energy discrepancy decreases under temporal refinement. It cannot
retroactively turn the dt25ns failed gate into a pass. No tolerance or model
change accompanies this study.

### 2026-10-27: fixed-step final-clock roundoff and saved-state recovery

The 32-cell, 12.5-ns-step trajectory reached step 16 with clock
`1.99999999999999938e-7 s`, then returned exit 1 against the requested
`2e-7 s` because the application used a strict final-clock comparison.
The failed original run remains recorded in
`build/pelef-full-rhs-time12-200-20261027/record.json`; it is not a successful
uninterrupted application run.

The EB-tree front end now uses one final-time predicate at both loop exit
and completion validation, allowing `8*epsilon*max(abs(final_time),abs(time))`.
It does not change the recorded clock, numerical state, timestep, or field
acceptance limits. The focused regression reproduces the original failure;
with the fix, serial and MPI2 complete 16 steps, MPI3 can restore the completed
serial checkpoint without a further step, and serial/MPI2 still reject a
15-step prefix. All five cases passed. The existing MOL cell-trace checker
also passed. Mechanism provenance (including the selected Simple bundle),
Python syntax, project contract and `git diff --check` passed. This is focused
application evidence, not a full configuration matrix.

Task-only patch, baseline failure, build/test logs, compiler/MPI identities
and binary/source hashes are in `build/pelef-eb-stop-roundoff-20261027/`.
The production checkpoint was restored into `termination-recovery3/` under
the trajectory directory. Its final CSV and cumulative budget CSV are byte
identical to the original step-16 snapshots, and the input checkpoint hash
is unchanged. Recovery attempts 1 and 2 are retained: the first input
incorrectly aliased restart/checkpoint paths; the second completed the app
but its postprocessor requested an absent unsuffixed budget filename.
Neither required rerunning or modifying the saved numerical trajectory.

The recovered run passes the original inventory gates. Of 17 field times,
0 through 112.5 ns pass; z momentum fails from 125 ns, total/internal energy
join it from 162.5 ns, and density also fails at 200 ns. At 200 ns,
z-momentum maximum error is `3.31563798908823e-6 kg m^-2 s^-1` at fine
(12,7,32); total-energy scaled error is `5.2119923486675483e-8` at fine
(12,8,32), against unchanged limits `1e-6` and `1e-8`. These exceed the
25-ns-step errors (`9.22008217564e-7` and `1.21649977175e-8`).
The first z-momentum failure at 125 ns occurs at fine (14,7,32); the maximum
locations over time are in `maxima-timeline.json`. Thus timestep halving
moves/amplifies reference discrepancies at the upper coarse/fine interface;
it does not establish convergence or identify the causal operator by itself.
The next diagnostic must capture these cells and their coarse donors,
rather than assuming the prior 25-ns worst cell remains authoritative.
The original failed-run record is unchanged; `recovery-record.json` records
state/ledger identity, 17 field comparisons and the inventory result.

### 2026-10-27: 12.5-ns interface discrepancy localized to transverse gradients

`build/pelef-time12-interface125-20261027/` captures fine (14,7,32),
(12,7,32), (12,8,32) through 125 ns. Candidate and native runs completed
in 1676.77 and 141.44 s respectively; all 11 candidate snapshots are byte
identical and all 11 native exports parsed-identical to the original
12.5-ns trajectory. The notification returned exit 1 only because the
pressure postprocessor rejected duplicate face records. All 360 pairs were
exact duplicates; the repaired parser retains one only after exact equality
verification. Both original failure and recovered analyses are retained.

At (14,7,32), the 112.5–125-ns z-momentum error increment is
`-1.55917842434e-6 kg m^-2 s^-1`. Integrated final-corrector face forcing
reproduces it with residual `-8.50e-17`. The first fine substep contributes
only `+1.8841e-10`; the second contributes `-1.55936683575e-6`.
Transport contributes `-1.60477623138e-6` across both substeps, including
`-1.54496783763e-6` from z faces. Local transport replays reproduce native
and candidate momentum fluxes within `9.5e-15 Pa`. Substituting native
face temperatures leaves the discrepancy essentially unchanged;
substituting all velocity gradients leaves `9.60e-16`, and substituting
only transverse gradients leaves `8.01e-9`. This attributes the local
transport discrepancy to gradients, not transport coefficients, without
claiming its upstream interpolation cause is proved.

The largest face discrepancy occurs at event 60, upper z face: `0.0774407 Pa`.
Adjacent halo cells (13/14,7/8,33) differ in x velocity by about
`0.0215937 m/s`, despite density differences of about `7.11e-9 kg/m^3`.
Pressure/flux/replay tables and the exact instantaneous/integrated definitions
are in `pressure-comparison.json`, `late-flux-comparison.json`,
`comparison.json`, and `transport-summary.json`. The step-9 checkpoint at
112.5 ns is retained. A one-step donor capture in
`build/pelef-time12-donors125-20261027/` is pending to distinguish temporal
coarse-state interpolation and common slope limiting; C2 remains open.

### 2026-10-27: donor limiter amplification verified on saved trajectory

The one-step `build/pelef-time12-donors125-20261027/` capture completed with
exit 0. Both restored candidate snapshots (steps 9/10) and all 11 native
exports match their ordinary trajectories. All 48 captured interpolations
replay, their marker clocks match the application observer, and reconstructed
candidate-minus-native halo velocity differences agree with the independent
transport-state capture within `4.45e-15 m/s`.

For fine (14,7,33) at native event 60, the y slope's common factor is
`0.22629602778325691` in the candidate versus `0.29958627088803685` in native;
z momentum controls both factors. The coarse center z-momentum discrepancy
is `7.81544e-9 kg m^-2 s^-1`. Replacing only the seven z-momentum donor
values with native values changes halo x-velocity error from
`0.02159367109` to `4.54139e-10 m/s`; replacing all other donor components
leaves `0.02159367063 m/s`. This is verified amplification of existing
state differences by the shared limiter, not evidence that loosening or
replacing the limiter is valid.

Coarse (7,4,17) endpoint z-momentum error increases from `9.86835e-12`
at 100 ns to `8.43658e-9` at 112.5 ns. At 125 ns its post-synchronization
error is `1.99309e-7`, distinct from the pre-synchronization interpolation
input above. `paired-donors.json`, `coarse-endpoint-errors.json`, and
`verified.json` retain inputs, substitutions, clock checks and identities.
A focused ordinary-output-equivalent capture of the earlier 100–112.5-ns
coarse/fine flux transfer is running in
`build/pelef-time12-sync112-20261027/`; no new numerical policy is adopted.

### 2026-10-27: 100–112.5-ns coarse synchronization closure

`build/pelef-time12-sync112-20261027/` completed with exit 0 and reproduces
10 candidate snapshots and 10 native exports. For coarse (7,4,17), the
lower-z synchronization correction differs by `8.42661149e-9`, matching
the observed donor error jump `8.42661219e-9` within `6.94e-16`.
The dominant incoming differences are fine (13,7,32) and (13,8,32),
respectively `5.35369e-12` and `5.28132e-12 Pa s`. Their transport parts
are `5.42120e-12` and `5.34828e-12 Pa s`; hydro contributions are smaller
and opposite in sign.

A local transport replay using the already saved states reproduces those
integrated differences within `1.49e-19 Pa s`. Native velocity-gradient
substitution reduces them below `6.7e-19 Pa s`; temperature substitution
has negligible effect. Transverse-gradient substitution leaves about
`-2.1e-14 Pa s` for each face. These close the preceding amplification
chain through accepted face transfer, without treating synchronization
itself as an implementation defect. At the upper halo, (12,7/8,33) has
x-velocity error `2.07169e-4 m/s` at event 54. The next one-step capture
uses the saved 100-ns checkpoint and donor (6,4,17) to locate that input
error; record directory `build/pelef-time12-donors112-20261027/`.

The sync run's copied `checkpoint_intent` text is stale: the actual input
uses interval 8 and retains step 8 at 100 ns. `verified.json` explicitly
records this correction and the checkpoint hash. Field tolerances and
numerical policies remain unchanged; C2 is still open.

### 2026-10-27: actual identical-input interpolation kernels agree

The donor-(6,4,17) run in `build/pelef-time12-donors112-20261027/`
completed with exit 0; two candidate snapshots and ten native exports match
ordinary output. At event 54, the x common slope factors are
`0.36846354585` versus `0.36963768103`, again controlled by z momentum.
Replacing only z-momentum donor inputs reduces halo (12,7,33) x-velocity
error from `2.07168934e-4` to `6.49078e-10 m/s`.

A direct executable comparison now calls the actual linked PeleF
`interpolate_eb_linear_halo_3d` and the pinned AMReX inline
`mf_cell_cons_lin_interp_llslope`/`mf_cell_cons_lin_interp` functions on
identical captured donor inputs. Across 384 combinations (two captured
windows, both solvers' donor inputs, SI and CGS arithmetic), all 16 output
components are bit identical. This covers interior ratio-2 interpolation
at the captured offsets; it does not qualify physical boundaries or the
full hierarchy. Sources, compile/link commands, binaries, input table and
hashes are in `build/pelef-halo-identical-input-20261027/`. The first native
harness compilation lacked the public wrapper header; the corrected harness
includes `AMReX_MFInterp_C.H` and compiles successfully. No production
interpolation change is warranted by these results.

The saved coarse-field timeline identifies an earlier seed region:
coarse (3,3,17) already has z-momentum error `1.41318e-10` at 12.5 ns,
versus only `7.65e-15` at (6,4,17). The former grows to `-1.02523e-9`
at 25 ns. A single initial-step capture at coarse (3,3,17) and its four
fine faces is running in `build/pelef-time12-origin12-20261027/`, to
identify the initial forcing discrepancy before the later amplification
chain. C2 and the original field/convergence requirements remain open.

### 2026-10-27: first-step seed localized before late-time amplification

`build/pelef-time12-origin12-20261027/` completed successfully and reproduces
both initial/final candidate snapshots and native exports for 12.5 ns.
At coarse (3,3,17), the reconstructed synchronization correction difference
is `1.41324e-10`, close to the final field error `1.41318e-10`; native
pre-synchronization state was not captured here, so this is not presented
as an independently observed native pre/post synchronization closure.
All four incoming fine upper-z differences are predominantly transport.
For fine (5,5,32), the transport contribution is `3.22211e-14 Pa s`;
independent replay gives `3.22210e-14`. Native gradient substitution leaves
`4.52e-21`, while temperature substitution has negligible impact.

The initial transport-state velocities match to about `2e-15 m/s` in the
captured neighborhood. At the first predictor, halo velocities differ by
up to `2.95e-6 m/s`; upper-z transport stress at (5,5,32) already differs
by `3.45046e-6 Pa`. This shifts the next diagnostic to the first coarse
update and its interpolation inputs, rather than another late-time
transfer. `build/pelef-time12-donors12-20261027/` is capturing donor (3,3,17)
and fine (5/6,5/6,33) for that initial step. No tolerance or production
numerical-policy change has been made.

### 2026-10-27: initial limiter seed and chemistry sensitivity discriminator

The first-step donor capture `build/pelef-time12-donors12-20261027/`
completed with exit 0 and ordinary-output identity. At event 2, coarse
(3,3,17) and its neighbors have z momentum around `3e-8` to `6e-8`;
candidate/native differences are only `2e-15` to `7e-14`. These change
the common y slope factor from native `0.84467723554919` to candidate
`0.8446725298248001`. At halo (6,5,33), x/y velocity differences are
`5.0709e-7`/`1.77484e-6 m/s`; substituting only the z-momentum donors
reduces them to `1.01e-12`/`1.38e-14 m/s`.

The root predictor at (3,3,17) has pressure discrepancy `3.1371e-7 Pa`
versus an initial difference `-1.4552e-11 Pa`. Captured donor species
and temperature also differ slightly (center T `3.24e-9 K`). This does
not yet assign the seed to chemistry: spatial arithmetic and chemistry
both remain possible contributors. Saved 25-ns sensitivity outputs with
native activation/native density-atol reduce post-step coarse (3,3,17)
z-momentum error from `2.60e-13` to `-5.07e-16`, although that same
sensitivity failed the later 200-ns energy gate and was not adopted.
A one-step, otherwise identical 12.5-ns sensitivity now runs in
`build/pelef-time12-nativechem12-20261027/` to distinguish its effect on
initial pressure and momentum. It uses the existing identity-checked
isolated executable; production tolerances and mechanism are unchanged.

Adding the initial-step donors to the actual-kernel replay expands it to
576 cases: the maximum output difference is one ULP, and maximum velocity
difference `1.14e-15 m/s`. The original 384-case bit-identical result is
preserved separately as `comparison384.json`; the expanded result does
not support bit identity for every initial-state input, but its rounding
difference is far below the observed halo velocity discrepancy.

### 2026-10-27: native chemistry parameters reduce the initial pressure seed

`build/pelef-time12-nativechem12-20261027/` completed in 240.20 s with
identical initial fields and both original field gates passing. Matching
both the pinned activation basis and native density absolute tolerance
reduces coarse (3,3,17) predictor pressure error from `3.13710e-7` to
`-1.01863e-10 Pa`; the corrected-stage error is `1.60071e-10 Pa`.
The synchronization correction difference falls from `1.41324e-10` to
`-6.22211e-12`. Across the domain, maximum z-momentum error decreases from
`2.86929e-9` to `3.04892e-10`, but total-energy maximum error increases
from `6.05392e-6` to `1.33833e-5 J/m^3`. This is a diagnostic improvement
of the localized seed, not uniform improvement or adoption evidence.
The previously measured 200-ns/25-ns-step sensitivity failure still applies.

Two independent one-step controls, activation-only and tolerance-only,
are running under `build/pelef-time12-chem-factorial-20261027/` to separate
the two effects. They use isolated links of the already generated native
activation archive or the already verified reactor tolerance object;
canonical mechanism sources, generated files and production tolerances
are unchanged. `tools/prepare_lidryer_native_activation.py` documents and
checks the source-derived calorie gas-constant conversion used by the
compatibility archive; it is not a physical calibration.

### 2026-10-27: activation-only versus tolerance-only controlled comparison

Both one-step controls completed successfully with identical initial fields
and original initial/final field gates passing. The factorial record is
`build/pelef-time12-chem-factorial-20261027/comparison.json`.

| Variant | Coarse predictor pressure error (Pa) | Global max z-momentum error | Global max total-energy error (J/m^3) |
| --- | ---: | ---: | ---: |
| Canonical | 3.13710e-7 | 2.86929e-9 | 6.05392e-6 |
| Native activation only | -1.55705e-9 | 3.62216e-10 | 1.20849e-5 |
| Native density atol only | 3.15180e-7 | 2.93276e-9 | 7.47709e-6 |
| Both | -1.01863e-10 | 3.04892e-10 | 1.33833e-5 |

Thus activation conversion, rather than tightening the density tolerance,
is the dominant discriminator for this initial pressure seed. The energy
maximum does not improve, and this single step is not trajectory acceptance.
Direct bundle comparison independently confirms identical species, thermo
and transport data, with only the 23 activation-energy entries differing
by the source-derived ratio `0.9999943013061792`. The pinned native
calorie gas constant is `1.9872155831739962 cal mol^-1 K^-1`, whereas the
canonical importer uses the standard SI/calorie conversion. This is a
reference compatibility distinction, not a proposed change to the canonical
SI mechanism. `verified.json` retains this audit and actual binary hashes.

The activation-only executable, retaining canonical solver tolerances, is
continuing its verified 12.5-ns checkpoint through 200 ns in
`build/pelef-time12-activation200-20261027/`. The run checks all 17 original
field times, original inventory gates, prefix byte identity and checkpoint
immutability. It has not yet passed those gates or been adopted. Existing
combined-parameter 25-ns/200-ns failures remain valid separate evidence.

### 2026-10-27: activation-only 32-cell trajectory passes through 200 ns

`build/pelef-time12-activation200-20261027/` completed its 12.5-to-200-ns
continuation in 2706.17 s. All 17 original fixed-PeleC field times and the
original inventory gates pass. The same-binary step-1 prefix is byte
identical, the restart is unchanged, and all recorded checker input hashes,
executable, input and final checkpoint identities were independently verified.
At 200 ns, z-momentum maximum error is `3.56160409e-8` (limit `1e-6`),
energy scaled error `5.10798271e-10` (limit `1e-8`), and density scaled error
`4.00439906e-10` (limit `1e-8`). Canonical SI errors at that time were
`3.31563799e-6`, `5.21199235e-8`, and above the original density limit.
`baseline-comparison.json` retains every time and field, not only maxima
selected for improvement. The reference solver/data and acceptance limits
were not modified; only the explicit candidate activation conversion differs.

This establishes a bounded native-activation compatibility trajectory,
not adoption of a replacement canonical mechanism, mesh/time convergence,
full MPI qualification, or physical validation. A fresh 25-ns-step trajectory
uses the same executable and tolerances in
`build/pelef-time25-activation200-20261027/`. The same activation archive
is linked into serial/MPI applications for nonuniform MPI2-to-MPI3 restart
and trace checks in `build/pelef-activation-validation-20261027/`.
Both records remain pending until their actual results are verified.

The compatibility serial/MPI checks subsequently completed with exit 0.
The nonuniform 8-cell MPI2 prefix and MPI3 continuation match the serial
prefix/final CSVs byte for byte; the checkpoint and both executable hashes
were independently verified. The existing MOL trace checker also passed.
`build/pelef-activation-validation-20261027/verified.json` retains the scope:
this is focused restart/diagnostic equivalence, not a 32-cell MPI/native
trajectory comparison or a full configuration matrix. The 25-ns serial
trajectory is still pending.

### 2026-10-27: 25-ns trajectory limitation and all-species temporal analysis

The fresh activation-only 25-ns trajectory completed in 2316.37 s. Original
field gates pass at 0–150 ns, but total/internal energy fail at 175 and
200 ns. At 200 ns, scaled energy error is `1.43247274267e-8` against
`1e-8`, with maximum at fine (26,14,32). Z-momentum maximum error
`7.63766219e-7` passes; inventory also passes. Checker-input hashes were
independently verified in `build/pelef-time25-activation200-20261027/verified.json`.
This does not revoke the separately verified 12.5-ns pass, but rules out
claiming unrestricted step-size acceptance for the compatibility setting.

`tools/analyze_eb_temporal_convergence.py` now computes same-time,
volume-weighted Linf/L1/L2 temporal self-differences for density, three
momenta, total/retained internal energy, temperature and all nine species
mass fractions. It requires matching composite keys/volumes, finite fields
and matching native/application clocks; it records every input hash and
only reports observed orders when three equally refined spacings exist.
It does not replace field or inventory gates. A two-cell unequal-volume
norm check, Python syntax, mechanism provenance, project contract and
`git diff --check` passed. The existing 25/12.5-ns data produced all16
observables at50/100/150/200ns; at200ns, H2 L1 self-differences are
`8.77714e-10` candidate and `8.77755e-10` native. Two spacings do not
establish an order.

A third spacing, 6.25 ns through100ns, is running in
`build/pelef-time6-activation100-20261027/`, with the original native solver,
17 field checks and original inventory gates. Its saved temporal spec
compares25/12.5/6.25ns at50 and100ns, where the coarser trajectories already
have verified reference agreement. The later25ns failures remain open;
no acceptance threshold was relaxed.

### 2026-10-27: three-spacing result and native reproducibility control

The 6.25-ns compatibility run completed with exit 0 and passes inventory,
but original field gates first fail at 50 ns (z momentum and total/internal
energy); density also fails from81.25ns. At100ns, z-momentum error is
`1.13160909e-6`, energy scaled error `1.56296465e-8`, density scaled error
`1.45324569e-8`. All17 checker-input hashes and inventory inputs were
independently verified. Thus the12.5ns trajectory pass does not establish
monotone improvement of reference agreement under time refinement.

The three-spacing analyzer completed for all16 observables at50/100ns.
At50ns, volume-weighted L1 observed orders range0.703–1.238 (candidate)
and0.703–1.238 (native); at100ns they range0.535–1.158 and0.534–1.166.
Linf orders are nonpositive for H2/O2/HO2/N2 at50ns and HO2/H2O2 at100ns
in both solvers. The largest candidate/native observed-order differences
across all fields and norms are0.00290 and0.00868 respectively. These
are bounded shared temporal trends, not all-norm convergence qualification.
`temporal-three-spacings.json` records complete metrics, masks, field
conventions and hashes; `verified.json` records the independent audit.

Further blind grid/time sweeps are deferred. A native-only control in
`build/pelec-roundoff-control-20261027/` reruns the identical pinned binary,
input and16-step schedule, alongside a run with each supplied initial
profile temperature moved upward by exactly one binary64 ULP. Other input
bytes are preserved. It will verify the repeat against all17 native exports,
measure the actual initialized perturbation, and compare the perturbed
trajectory with the unchanged original field limits. This assesses native
sensitivity; it does not replace the frozen baseline or authorize any gate
relaxation. C2 remains open.

### 2026-10-28: native roundoff control recovered; pressure arithmetic isolated

Both native controls completed normally on2026-10-27 (repeat137.20s,
perturbed139.86s), but their export/comparison process disappeared after
step0. No live worker remained on resumption. Only the saved plotfile
exports and comparisons were recovered; the original record remains
unchanged. `recovery-record.json` and `verified.json` under
`build/pelec-roundoff-control-20261027/` verify all17 repeat exports are
parsed-identical to baseline. The33 profile temperatures each differ by
exactly one upward binary64 ULP, with all other profile fields unchanged.
The actual initialized maximum temperature difference is `9.09495e-13 K`.
At100ns, perturbed native z-momentum difference is `1.27552e-7` and scaled
energy difference `1.76138e-9`; all17 times remain inside the original gates.
Thus native sensitivity exists, but this control does not explain away the
larger candidate failures or justify relaxed thresholds.

The temporal analyzer now has a known linear-time, unequal-volume fixture
covering SI/CGS conversion, native non-composite-cell exclusion, species
mass fractions, all norms/orders, and rejection of mismatched clocks.
`unit_eb_temporal_analysis` is registered with CTest and passed1/1 after
reconfiguration; direct unittest, Python syntax, project contract and
`git diff --check` also passed. No broad numerical matrix is claimed.

`build/pelef-pressure-arithmetic-20261028/` directly calls current
`mixture_pressure` and pinned `CKPY` on2592 identical captured rho/T/Y
inputs. Current evaluation differs by at most `4.36557e-11 Pa`; a direct
SI sum/product order reduces that to `2.91038e-11 Pa`, while pinned CGS
units/product order matches all2592 cases exactly. This isolates pressure
arithmetic only, not inverse-temperature recovery or coupled flux parity.
The harness's first compile lacked the selected-mechanism module include;
adding that include resolved it, with both attempts retained.

An isolated mixture-pressure object using the matched CGS product order
is running through50ns at6.25ns steps in
`build/pelef-cgs-pressure50-20261028/`. It links the verified activation
compatibility archive with unchanged solver tolerances and checks all9
original field times plus inventory. The original production physics
source is unchanged; this is a causal pilot before any adoption decision.


### 2026-10-28: CGS pressure pilot rejected; inverse-EOS replay audited

The isolated `build/pelef-cgs-pressure50-20261028/` application completed
8 steps /50ns in1040.316s. All9 field comparisons were independently rerun
against the original reference and unchanged limits; steps0–7 pass, while
step8 fails z-momentum and total/held internal energy. Inventory passes.
At50ns, maximum z-momentum error is `1.4009952149e-6 kg m^-2 s^-1` and
maximum total-energy error is `0.0056654770 J m^-3`, respectively1.21801
and1.21769 times the activation-only errors. `verified.json` records
source/binary/input and comparison hashes. Production mixture source bytes
still match the pre-pilot build record. The pressure change is **not adopted**.
The isolated prototype also needs failure-flag auditing before any future
reuse; no robustness qualification is inferred from this trajectory.

`build/pelef-inverse-eos-replay-20261028/` calls actual linked Fortran
`temperature_from_internal_energy` and pinned `GET_T_GIVEN_EY`, followed
by each implementation's pressure function. Inputs are the1848 saved center
and six-neighbor states from the12.5/112.5/125ns donor captures, including
both candidate and native states. Each kernel receives the same SI
internal-energy target `(rhoE - |momentum|^2/(2*rho))/rho`, mass fractions,
and stored temperature seed; the native target is converted to CGS.
The Fortran replay enables pinned extrapolation but not converged inverse.
Maximum differences are `1.8189894035e-12 K` and `8.7311491370e-11 Pa`.
The audit checks all1848 row IDs, finite positive temperatures/pressures,
native return statuses, all recorded input/source/binary hashes, and
independently recomputes the summary from saved raw outputs. This excludes
a larger inverse-EOS discrepancy on these inputs only; it does not exclude
amplification of small differences during coupled evolution.

A new isolated native diagnostic in
`build/pelec-pressure-order-control-20261028/` replaces only the nonrotating
primitive pressure evaluation in `Diffusion.cpp` with the candidate SI
sum/reciprocal/product order. Temperature, composition, sound speed/qaux,
and other native primitives remain native. Original `MOL.cpp` is compiled
alongside the diagnostic object. The frozen native source/binary and
reference exports remain unchanged. The control uses the original32-cell,
6.25ns-step inputs through100ns and will compare all17 outputs with the
original limits. It measures pressure-only arithmetic sensitivity, not a
replacement reference or complete candidate-EOS equivalence. The notification
record owns job identity and continuation. C2 remains open.

Provenance preflight for `mechanisms/h2o2_full.json` passed with upstream
SHA256 `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Diagnostic Python syntax, `git diff --check`, and
`python3 tools/check_project_contract.py --source .` passed. The first
contract-check invocation omitted required `--source` and returned a usage
error; the corrected invocation passed. No numerical production change or
full configuration matrix is claimed for this diagnostic increment.


### 2026-10-28: native pressure-order control reproduces the failure pattern

The pressure-only native control completed16 steps through100ns in126.570s,
with the native terminal marker and final plotfile present. All17 comparisons,
source/binary identities, four unchanged input files and exported CSV hashes
were audited in `build/pelec-pressure-order-control-20261028/verified.json`.
The initial fields are identical. The control passes through43.75ns, then
fails z-momentum and total/held energy from50ns; density also fails from81.25ns.
This is the same failure sequence as the activation-compatible candidate.

At50ns, native-control z-momentum difference is `1.3275718036e-6 kg m^-2 s^-1`
versus candidate `1.1502336451e-6`, both maximal at fine(25,9,32). Energy
absolute differences are `0.0053606793` and `0.0046526351 J m^-3`, both
maximal at fine(26,10,32). The unweighted signed error-vector cosine over
32256 composite cells is0.99894157 for z-momentum and0.99988341 for energy.
Here each vector is that run minus the unchanged pinned native export,
with SI conversion; cosine is dot(a,b)/(norm(a)*norm(b)), not a fitted
causal fraction or volume-weighted norm. At100ns the cosines remain0.99815444
and0.99924787. The first6.25ns differences were much less aligned, so these
late values establish similar amplified patterns, not identical seeds.

Thus pressure arithmetic alone in native can trigger this failure pattern.
This strengthens shared-algorithm sensitivity evidence beyond the earlier
one-ULP initial-temperature control, which stayed inside all gates. It does
not prove that all candidate discrepancies are roundoff or authorize changing
limits or replacing the original reference. No production change is adopted.

The next diagnostic adds read-only interpolation donor capture to original
and pressure-order native binaries at coarse(13,5,17) and(14,4,17), beside
the observed upper-interface maxima. Both runs stop at50ns and must reproduce
all9 respective uninstrumented exports exactly. Their matched context/time/
stage stencils will be replayed with original z-momentum substitution to
quantify shared-factor amplification. Artifacts and notification are under
`build/pelec-pressure-order-donors-20261028/`. C2 remains open.


### 2026-10-28: native pressure perturbation reaches the shared halo limiter

Both donor-capture native runs completed through50ns and independently
reproduced all9 respective uninstrumented exports,18 tables total, in
`build/pelec-pressure-order-donors-20261028/`. At event41, fine-level interval
start43.75ns /dt3.125ns /stage1, coarse(14,4,17)'s y-direction common limiter
changes from0.211480329465 to0.231677594114. In both runs z-momentum alone
controls that factor. At offset(-0.25,+0.25,-0.25), replayed halo velocity
differences are `(0.00293150858,-0.00668319765,-1.00465541e-6) m/s`.
Replacing only the center/six-neighbor z-momentum values with original
native values reduces the x/y differences to `-6.24e-12` / `2.53e-10 m/s`.
The center z-momentum perturbation is `-6.12425617e-8 kg m^-2 s^-1`.
These are matched-stage stencil interventions, not changes adopted in CFD.

A follow-on diagnostic in `build/pelec-pressure-order-flux-20261028/` records
transport and total face fluxes in all three directions at fine(25,9,32),
(26,10,32),(27,7,32), with pressure/time/dt/stage identity. Both original and
pressure-order runs must retain all18 output-table matches. The prepared
analyzer will integrate old/corrector flux differences across each pair of
fine intervals and compare them to momentum-difference growth. This tests
the downstream budget closure, without changing original thresholds.

Actual compiled AMReX/PeleF halo kernels were also replayed on1152 original,
perturbed and z-momentum-substituted stencils/offsets. Maximum component
difference is4 ULP and maximum velocity difference is2.274e-15 m/s; the
independent formula exactly reproduces native velocities on these cases.
`verified.json` records identities and all18 table hashes. This is bounded
interior interpolation equivalence, not full boundary/trajectory qualification.


### 2026-10-28: native pressure-order momentum growth closes from face fluxes

Both directional-flux captures preserve all18 respective native exports.
Source/input/binary/trace/export hashes and native final markers were audited
in `build/pelec-pressure-order-flux-20261028/verified.json`. All72 momentum
increments (8 root steps x3 full-volume fine cells x3 components) close from
recorded total face fluxes, with maximum absolute residual
`2.517896756453e-16 kg m^-2 s^-1`. Each root step integrates two fine
intervals using half old plus half final-corrector divergence. Transport is
phase0; hydro is total phase1 minus phase0. Native extensive fluxes are divided
by native cell volume and converted to SI momentum; no fitted scale is used.

At fine(25,9,32),43.75→50ns z-momentum difference grows by
`1.239160655881e-6 kg m^-2 s^-1`; integrated flux difference is
`1.239160656051e-6`, residual `-1.6971e-16`. The z-direction transport
contribution is `1.214801538834e-6`; z-hydro opposes it at `-2.424921384584e-8`.
The x/y transport contributions are `2.43453e-8` and `2.46734e-8`.
This closes the downstream native arithmetic-perturbation budget. It is not
proof that every candidate seed difference has the same origin.

In particular, at the first6.25ns the candidate z-momentum maximum was
`3.63150e-9`, versus `1.76058e-11` for the pressure-order native control, and
their error vectors were not aligned. The next focused capture therefore
compares the unchanged activation-compatible candidate with original native
at fine(32,8,32),(32,7,32),(31,8,32) through only6.25ns. It records each
stage's pressure/hydro flux and reaction forcing, and requires all4 initial/
final outputs to reproduce prior results. Artifacts/notification:
`build/pelef-time6-first-stage-20261028/`. No production algorithm or gate
is changed; C2 remains open.


### 2026-10-28: first-step candidate discrepancy is also transport-dominated

`build/pelef-time6-first-stage-20261028/` completed the first6.25ns and
reproduced all4 initial/final candidate/native outputs. Binary/trace/export
identities were audited. Stage0/2 forcing closes all9 momentum changes at
fine(32,8,32),(32,7,32),(31,8,32); energy residuals are retained separately
because energy may undergo recovery outside the captured forcing.
At(32,8,32), the two fine intervals contribute `1.21631165e-9` and
`2.41518697e-9 kg m^-2 s^-1` to z-momentum difference. Their non-hydro
contributions are `1.22161690e-9` and `2.43611897e-9`; hydro is small and
opposes the difference. Initial cell pressure agrees exactly at the first
two probes. Later stage pressure differences do not by themselves identify
the transport source, so no pressure-only production fix follows.

The next capture records both solvers' transport primitives/states over a
6x6x6 neighborhood of(32,8,32), plus its six face fluxes. Actual transport
kernel replay and native-temperature/gradient substitution will discriminate
the local source. The original activation-compatible binary is reproduced
through an isolated diagnostic tree object; no production numerical edit.
Artifacts/notification: `build/pelef-time6-first-transport-20261028/`.
C2 remains open with unchanged original gates.


### 2026-10-28: first predictor mismatch comes from transverse velocity gradients

The first6.25ns transport capture preserves all4 prior initial/final exports.
Actual-kernel replay in `build/pelef-time6-first-transport-20261028/` reproduces
native and candidate momentum face fluxes to at most7.11e-15 Pa. Native
energy-flux residual reaches1.25e-9 W m^-2; it is retained separately and not
used for the momentum attribution. Hashes and four output identities were
audited in `verified.json`.

For fine(32,8,32)'s upper z-face, event1 agrees to roundoff. Event2 (first
predictor) has a candidate-minus-native z-momentum stress difference
`-1.17620074516e-4 Pa`. Substituting native temperatures leaves this unchanged;
substituting all native velocity gradients gives `-8.60e-15 Pa`. Substituting
only transverse gradients gives `5.74e-14 Pa`, whereas normal gradients
leave `-1.17620074582e-4 Pa`. Thus the first discrepancy is in transverse
velocity reconstruction, not the transport coefficient's temperature input.
The relevant transverse derivative differences include du/dx=0.07611146 s^-1
and dv/dy=0.00070592 s^-1. Across the captured neighborhood, the largest
velocity difference is at corner halo(33,8,33): du=4.75340002e-5 m/s,
dv=5.88077732e-4 m/s, while initial event1 velocities agree to roundoff.

The next paired donor capture targets coarse(17,4,17), with fine halo
(33/34,7/8,33). Candidate logs retain ancestor old/new states, time fractions,
center/six neighbors and bracketing stage/reaction markers; native retains
matching interpolation context. Artifacts/notification:
`build/pelef-time6-first-donors-20261028/`. The purpose is to distinguish an
input-state/temporal-interpolation discrepancy from shared-limiter amplification,
without tuning the limiter or replacing the fixed native baseline. C2 remains open.


### 2026-10-28: H2O, not momentum, controls the first corner-halo mismatch

`build/pelef-time6-first-donors-20261028/` preserves all4 ordinary outputs.
The48 candidate stencil replays and6 native contexts were checked against
fine-level trace stages. The initial analysis expected coarse-stage CSV
rows despite configuring only fine probes; its marker check failed. The
corrected check compares the6 configured fine stages and independently
requires native context time=start for stage0 and start+dt for stages1/2,
with identical level/stage/dt. No time-nearest assignment is used.
Candidate ancestor time fractions are0,0.5,0.5,0.5,1,1; center-state temporal
mix residual is exactly0 on all captures. This does not independently
validate all native ancestor states, but no candidate time-mix discrepancy
is seen here.

At coarse(17,4,17), H2O alone controls the y-limiter after the predictor:
candidate0.616919955726 vs native0.618494885437. The center H2O density
at fraction0.5 is `4.19899469515e-12 kg m^-3`, differing from native by
`3.35026140816e-16 kg m^-3`. Its neighboring differences are almost uniform,
but small differential differences affect the limiter ratio of very small
slopes. For halo(33,8,33), replacing only H2O in the center/six neighbors
reduces du/dv differences from `4.75334437e-5` / `5.88070848e-4 m/s` to
`-2.998e-15` / `8.882e-16 m/s`. Replacing z-momentum has no such effect.
This distinguishes the initial trace-species-controlled seed from the later
z-momentum-controlled amplification identified in native controls.

The next controlled run retains this same instrumented activation-compatible
binary and changes only input `atol` from1e-14 to1e-15. This was initially
misidentified as native density-atol matching; the caller multiplies input
`atol` by old density, as corrected in the following entry. It will
compare original fixed-reference field gates and the matched donor records.
Temperature-coordinate tolerance and integration trajectories are not claimed
to become identical solely from this density-unit match. No production
threshold or acceptance criterion is relaxed. Artifacts/notification:
`build/pelef-time6-native-atol-donors-20261028/`. C2 remains open.


### 2026-10-28: tolerance interpretation corrected; tighter-input pilot rejected

`src/reactive/reactive_3d_mod.F90:509` passes `atol*old_state(1,i,j,k)` to
the full-coordinate reactor. Therefore the completed directory named
`build/pelef-time6-native-atol-donors-20261028/` is a tenfold-tighter *input*
atol trial, **not** a native density-atol match. At the captured donor, the
original effective density tolerance is `1.038175325166e-15 kg m^-3`; the trial
uses `1.038175325166e-16`. Native1e-18 g cm^-3 corresponds to1e-15 kg m^-3.
The unchanged binary and exact one-option input delta were audited, with the
correction preserved in `verified-correction.json`; original run records remain.

The trial reduces the local corner-halo dv difference from5.88071e-4 to
4.11853e-4 m/s and changes the y-factor to0.617391889805 (native0.618494885437),
but fails original final density/x/y-momentum/energy gates after6.25ns.
Maximum x-momentum difference becomes2.34903e-5 kg m^-2 s^-1 and energy
absolute difference0.01641573 J m^-3. It is not adopted. A local reduction
does not establish global improvement or justify acceptance changes.

The follow-on `build/pelef-time6-fixed-density-atol-20261028/` links the
existing isolated reactor object that supplies constant1e-15 directly to
CVODE, keeping input `atol=1e-14` and all other settings. This is the actual
fixed-density tolerance control, still diagnostic only. The source differs
only at the full-coordinate `FCVodeSStolerances` call; all original gates
remain unchanged. Compare the donor response and global fields before any
further implementation decision. C2 remains open.


### 2026-10-28: fixed effective density-atol passes first-step gates, not halo equality

The actual constant1e-15 kg m^-3 control in
`build/pelef-time6-fixed-density-atol-20261028/` completed6.25ns and passed
both original initial/final field gates. Its input is unchanged except output
paths; source audit confirms the isolated reactor changes only the full
CVODE scalar-atol argument. Binary/object/trace/checker hashes were verified.
The corner y-limiter becomes0.619758431131 (native0.618494885437); du/dv
errors reverse sign to `-3.81355e-5` / `-4.71802e-4 m/s`. Global z-momentum
maximum is2.88970e-9 and energy maximum3.83742e-5 J m^-3. Thus matching
this tolerance does not remove the trace-H2O limiter seed. This is a bounded
sensitivity result, not a promoted setting or C2 closure.

The next discriminator captures actual pinned native reactor inputs and
endpoints at coarse(17,4,17) plus its six neighbors, with original outputs
required unchanged. The existing actual full-coordinate Fortran reactor
will replay identical rhoY/T/energy/external sources/duration using the
activation-compatible mechanism and both density-scaled and constant atol.
This separates reactor error from upstream forcing differences without
another full candidate trajectory. Artifacts/notification:
`build/pelec-time6-donor-reactor-20261028/`. Original reference and gates remain.


### 2026-10-28: identical reactor inputs reproduce the limiter seed

`build/pelec-time6-donor-reactor-20261028/` records both root reactions at
7 coarse donor cells, preserving both original exports. Actual Fortran
full-coordinate replay uses the same native rhoY/T/energy/external sources
and duration for28 cases (14 cell/events x2 atol policies). Native terminal
markers, source/binary/trace/export hashes and replay finiteness were audited.
Maximum H2O endpoint difference is7.81397279e-16 kg m^-3 with density-scaled
atol and8.68113128e-22 with native constant density-atol. The first value
matches the full candidate donor discrepancy's magnitude, so this seed need
not be attributed to a new hydro forcing defect.

For native root reaction event4 at coarse(17,4,17), H2O neighbor differences
are `+2.76190243e-20` and `-1.51001850e-19 kg m^-3` along y. Replacing only
H2O with identical-input Fortran endpoints predicts common factors
0.616951110074 (density-scaled) and0.619755812325 (constant), close to the
full-run factors0.616919955726 and0.619758431131. Native remains0.618494885437.
Thus very small chemically admissible endpoint differences act on much
smaller spatial differences and alter the common slope factor. This is an
ill-conditioned shared-factor diagnostic, not evidence that either endpoint
chemistry is physically wrong, nor permission to ignore original gates.

`tools/analyze_eb_limiter_conditioning.py` now audits explicit paired SI
center/six-neighbor stencils. It reports directional slopes, controlling
components, factors, donor differences, velocity differences and every
single-component substitution. It requires pre-matched contexts, refuses
nonfinite/invalid-density/shape/offset inputs, and emits qualification=false
with input/tool hashes. No component threshold, limiter policy, tolerance
or numerical production code is changed. `conditioning-input.json` and
`conditioning.json` retain the four actual/full-run and same-input H2O
comparisons. Synthetic near-uniform trace-species control and refusal tests
pass (2 unittest methods); Python syntax passes. CTest is registered as
`unit_eb_limiter_conditioning`; its separate result is recorded below.

C2 still requires the original trajectories and convergence gates. The
remaining implementation question is whether native-compatible RHS/Jacobian
arithmetic removes a material part of these same-input reactor differences;
shared-limiter sensitivity alone does not justify tuning or replacing the
fixed reference. Broader C3–C5 qualification is unchanged.

CTest `unit_eb_limiter_conditioning` passed1/1 after reconfiguration.
Mechanism provenance, project contract and diff whitespace checks also passed.
Exact tool/test identities, Python version and commands are retained in
`build/pelec-time6-donor-reactor-20261028/tool-validation.json`. No broader
configuration matrix was rerun. No local diagnostic job remains pending.


### 2026-10-28: Jacobian diagonal ablation and equilibrium-reference arithmetic

Source inspection of pinned CPU dense `cJac` finds no write to the T,T
entry, whereas PeleF computes the full temperature column by centered
independent-temperature RHS differences. A one-entry diagnostic sets only
PeleF's T,T entry to zero in `build/pelef-native-jtt-replay-20261028/`.
Same-input H2O maximum error remains7.81397278e-16 kg m^-3 (density-scaled)
and8.67945107e-22 (constant atol), effectively unchanged from the prior
7.81397279e-16 /8.68113128e-22. This is not adopted and does not establish
full native Jacobian equivalence. The first build omitted the pinned
SUNDIALS_INT32_T macro and failed with integer-kind mismatches; the failed
log is preserved, and the corrected compile/link and replay completed.

`build/pelef-donor-rate-jacobian-20261028/` then compares28 captured pre/post
states at identical fixed T/rhoY. It calls native CKWYR and the same
CKYTCR/DWDOT route used by CPU cJac for the9x9 density block, applying
molecular-weight conversion, versus actual Fortran production and analytic
density-Jacobian functions with activation-compatible kinetics. Maximum
norm-scaled production difference is8.77547e-11; density-Jacobian difference
is1.19717e-10. Here scaling divides the maximum absolute entry difference
by the largest native absolute entry per state, then takes the maximum over
states. This comparison excludes forced energy recovery and the temperature
column; it is not a full RHS/Jacobian certification.

Pinned LiDryer uses8.31446 for equilibrium reference concentration, while
PeleF uses8.31446261815324. An isolated kinetics object changes only that
reference-concentration denominator to8314.46 in SI kmol units, retaining
canonical production source. In `build/pelef-native-equilibrium-donor-20261028/`,
the fixed-state production/density-Jacobian discrepancies fall to3.12257e-14
and1.97870e-13. Constant-atol same-input H2O endpoint difference falls to
1.08149e-22 kg m^-3, but the inferred y-factor remains0.619714595065 versus
native0.618494885437 (prior constant-atol factor0.619755812325). Reducing
absolute chemistry error therefore does not remove the ill-conditioned
spatial-ratio discrepancy. No production equilibrium convention, tolerance,
limiter policy or acceptance threshold is changed.

All three records' source/input/binary hashes were checked in their
`verified.json` files. Diagnostic Python syntax, project contract and
`git diff --check` passed. No broad matrix or physical validation is claimed.
The remaining full temperature-column/Jacobian convention should be compared
before deciding on a compatible numerical correction; C2 remains open.


### 2026-10-28: temperature-column, thermal-row and CVODE-coordinate controls

`build/pelef-temperature-column-20261028/` compares the nine species entries
of the temperature column on the same28 captured states, including their
external species source in finite-difference arithmetic. Native analytic
DWDOT values are converted to SI density-rate/K. At production deltaT
=max(1e-4,1e-6*T), maximum column-norm-scaled differences are5.44098e-10
for native rates differenced numerically,5.01564e-10 for activation-compatible
Fortran, and5.54131e-10 with the native equilibrium-reference constant.
Halving deltaT increases these to1.72651e-9,1.48035e-9,1.46817e-9 respectively.
This distinguishes analytic-versus-centered arithmetic from a uniquely
Fortran temperature-column defect. No deltaT change is adopted. The thermal
row is explicitly outside this column comparison.

`build/pelef-native-thermal-row-replay-20261028/` separately removes external
forcing from the approximate thermal density row, computing its chemical
thermal rate at the trial temperature, and keeps the native-style zero T,T
entry. The same28 reactor cases retain maximum H2O discrepancies7.81397e-16
and8.67985e-22 kg m^-3 for density-scaled/constant atol. This does not materially
change the prior endpoint difference and is not promoted.

`build/pelef-cgs-vector-replay-20261028/` uses internal CVODE density
coordinates in g cm^-3, converting inputs, species derivatives, cross-block
Jacobian factors, endpoint admission and density atol consistently. Physical
thermodynamics/kinetics remain SI. The scalar absolute tolerance on temperature
is also scaled down, as recorded in the scope; this is not a pure temperature-
tolerance invariant experiment. Density physical tolerance remains unchanged.
Both input policies complete all28 cases, with H2O maxima7.81397e-16 and
8.68908e-22 kg m^-3. Coordinate conversion therefore does not resolve this
captured seed either. No production code or checkpoint policy is changed.

Each result has verified source/input/binary hashes. Python syntax, project
contract and diff whitespace checks passed. These are bounded controls,
not a whole-CVODE matrix, trajectory qualification or physical validation.
Together they support the trace-H2O shared-limiter conditioning diagnosis;
the original n32 short-step trajectory/convergence failure remains open.
Further arbitrary tolerance/rounding sweeps are not justified by these data.


### 2026-10-28: larger-grid MPI trajectory qualification started

A missing CPU/MPI evidence item is the n32/12.5ns/200ns application trajectory
that passed in serial. The frozen activation-compatible MPI executable from
`build/pelef-activation-validation-20261027/` is reused only after binary and
mechanism-archive SHA verification; no new current-source rebuild is claimed.
The fresh initial state, physics settings, native reference contract and17
output times are unchanged. Completion requires original native field gates,
byte identity with the recorded serial trajectory, and the original inventory
checker. This comparison cannot close the separate6.25ns/convergence failure.

The first attempt in `build/pelef-mpi2-time12-activation200-20261028/` exited
before simulation with `MOL cell trace currently requires serial gas`: the
serial diagnostic input still contained five cell probes. Its failed record
and log are preserved. The ordinary rerun in
`build/pelef-mpi2-time12-activation200-ordinary-20261028/` removes only the
serial-only trace block, with outputs redirected to its own directory.
No numerical setting is changed. Notification metadata records command,
PID, logs and continuation; results remain pending until independently checked.

While this runs, `docs/pelec_mapping.md` now routes current 3D tree, MPI,
boundary, selected chemistry and persistence implementations separately from
older bounded-path tables. Every explicit new source path was checked live.
The roadmap now routes C2 case families to their existing contracts/checkers;
it does not silently finalize C0 or reduce the release matrix.


### 2026-10-28: larger-grid MPI trajectory independently verified

The ordinary MPI2 run completed 16 steps at 1.9999999999999994e-7 s in
2954.004 s. All 17 original native field comparisons were rerun, each over
32256 composite fluid cells; every gate passes. All 17 CSVs are byte-identical
to the recorded serial trajectory, all output clocks match, and the final CSV
matches step 16. The original inventory checker was rerun and passes, including
mass, elements, energy correction, ledger reconciliation and H2 consumption.

Reproduction: `python3 build/pelef-mpi2-time12-activation200-ordinary-20261028/verify_completion.py`.
That audit stores exact commands, Python version, audit HEAD, 64 artifact hashes
and per-step results in `verified.json`; fresh field reports are under
`independent-checks/`. Binary, input, initial state and reference-contract hashes
match the launch record. The frozen MPI binary SHA256 is
`16f7bab5b19c4b2bc3251452871a412bc29d34faf65e1fc8b78ce1ca8d4c6ad9`.
This uses the previously built explicit native-activation-compatible mechanism;
it is not a canonical-mechanism or current-source configuration-matrix claim.
The log includes display-authorization warnings before normal numerical
completion; no application failure is reported. The earlier serial-trace input
refusal remains recorded separately and is not counted as a successful run.

LiDryer provenance, project contract, Python syntax and diff whitespace checks
pass. No numerical kernel or acceptance threshold changed for this audit.
The n32/6.25-ns fixed-reference failures and maximum-norm convergence remain
open, as do the independent nonreacting AMR/MOL2 register/conservation issue
and C3–C5. This closes only the pending larger-grid MPI equivalence check.


### 2026-10-28: close the donor-to-halo difference algebraically

`tools/analyze_eb_limiter_conditioning.py` now decomposes each interpolated
conserved-state difference into the central value and three directional
factor/slope terms. It uses the explicit reference-anchored identity
`fc*sc-fr*sr = (fc-fr)*sr + fc*(sc-sr)`, multiplied by the child offset.
The cross term belongs to the slope term by convention; these are algebraic
contributions, not independent causal shares. Velocity differences are then
reconstructed as `dm/rho_candidate - velocity_reference*drho/rho_candidate`.
Both state and velocity closure residuals are retained rather than assumed zero.

Applied to the four existing stencils in
`build/pelec-time6-donor-reactor-20261028/conditioning-input.json`, the maximum
velocity closure residual is7.444e-16 m/s. In the ordinary first-step stencil,
the y-directed shared-factor contribution to interpolated y momentum is
6.105205977e-5 kg m^-2 s^-1; its slope-change term is-1.113e-16 and the central
change is-6.661e-16 in the same units. The complete transverse velocity
difference is5.880708070e-4 m/s; the density-denominator contribution is only
4.007e-12 m/s. This independently closes the already observed shared-factor
amplification rather than assigning it to a changed resolved momentum slope.

Results and source/input hashes are in `conditioning-decomposition.json` and
`decomposition-validation.json` in that directory. Three direct unit methods
pass, including an analytic case with simultaneous density, central-momentum,
factor and slope changes; the registered CTest passes1/1. LiDryer provenance,
project contract, Python syntax and diff whitespace checks pass. Existing
single-component substitution results are retained. No interpolation formula,
chemistry tolerance or original trajectory gate changed, and C2 remains open.


### 2026-10-28: separate observer outputs from held-state halo inputs

The previously reported `unit_eb_tree_mol_spatial_halo_3d` failure was reproduced
with the current-source library: `retained energy output preserves tree gas
rate`. The fixture compared a call without held-state input against one that
also supplies `internal_energy_state` and old held halos. The latter path
passes coarse geometry to halo filling and selects spatial coarse interpolation;
the former does not. It is therefore not an observer-only change of arguments.

The test now checks exact gas-rate/face-flux invariance in a separate call that
adds only internal-energy diagnostic outputs. The held-state case still checks
rate/halo equality, face-divergence conservation, direct halo filling, same-level
donors, the analytic half-increment from ancestor interpolation, local ghost
RHS, StateRedist inputs and transactional rollback. Direct halo/kernel calls
now explicitly use the same coarse geometry, preserved temperature seed and
total-energy primitive policy as the tree path. Exact assertions and existing
tolerances are retained. No numerical production source is changed.

The isolated revised fixture, linked to the same newly built library as the
failing original, passes. Source snapshots, task-only patch, compile/link
commands and both logs are in `build/eb-halo-observer-contract-20261028/`.
The first CTest attempt found no built executable and was Not Run, not a
reproduced numerical failure. A premature direct launch also found no executable
while the build was still running; the completed build was subsequently
confirmed before obtaining the actual baseline failure above.
The adopted test was rebuilt and the registered CTest passed1/1.
`validate.py` and `verified.json` record commands, timing and source/library/
binary hashes, which were independently rechecked after completion. Provenance,
Python syntax, project contract and diff whitespace checks pass. This is a
test-contract correction, not a trajectory, MPI or C2 qualification result.


### 2026-10-28: localize recursive non-MOL MPI timeout cost

A timing-only fixture linked to the frozen Debug archives of the failed suite
again reaches300 s without finishing. Both ranks complete the first five MPI
updates and the serial reference for layout3/step2. Serial references cost
22.1–22.3 s per step; completed MPI updates cost24.3–26.3 s. Thus the log shows
continued stage progress, not a fixed early collective deadlock. It does not
prove the final MPI update or rollback completes. CPU samples are retained,
but MPI polling can consume CPU, so CPU usage alone is not used as progress.
Artifacts: `build/mpi-nonmol-timing-20261028/record.json`, `run.log`, and
`timing-analysis.json`; the driver records timeout explicitly even though its
wrapper exits0. No timed-out test is counted as passed.

`check_recursive_registers` in `tests/unit/test_mpi_eb_tree_full_3d.F90` now
caches the two serial endpoints and their boundary/wall/chemistry/source,
limiter and substep-count results for `nonmol_rec` only. Each rank computes
those references independently during the first layout. The serial call has
no owner map or callback state; every later layout starts from the same seed.
The cache therefore removes four repeated serial advances per rank while
retaining all six MPI updates, three owner layouts, original field/budget/count
comparisons, and every layout's failure/rollback checks. Other modes and all
numerical production code are unchanged; no time limit or tolerance is raised.

The isolated modified fixture passes with two and three ranks in221.31 s and
215.47 s, sequentially, against the same frozen Debug archives and unchanged
300 s limits. All six layout/step comparisons and each final rejection/rollback
complete. The previously completed five budget residuals match the timing
baseline exactly; the sixth is4.2240294075782335e-16 for both rank counts.
`build/mpi-nonmol-reference-reuse-20261028/verified.json` records independent
checks of archive/source/binary identity, complete layout/step output, final
PASS markers and prior residual identity. This establishes that reference
recomputation cost was sufficient to cause the observed300 s timeout in this
fixture; no communication or numerical-kernel repair was needed for these runs.

A normal rebuild from current source and rerun of the original registered
MPI2/3 CTests is now active in `build/mpi-nonmol-current-reuse-20261028/`.
This step must be verified separately before calling the current-source
focused tests passed. The source manifest, build/CTest commands, logs and
notification state are retained. The frozen-library result is a test-cost
correction, not a current-source MPI release-matrix or C2 qualification claim.


The subsequent normal current-source build succeeds in70.38 s. Original
registered MPI2/3 recursive non-MOL CTests pass2/2 in220.61/214.30 s under
the unchanged300 s limits. All12 layout/step budget residuals exactly match
the isolated result, both final PASS markers are present, and every recorded
source/cache hash plus executable identity was independently rechecked.
`build/mpi-nonmol-current-reuse-20261028/verified.json` records that audit.
The earlier two timeouts are resolved for this current-source focused rerun;
this does not retroactively change old failures or qualify the full MPI matrix.


### 2026-10-28: remove intermediate overflow in falloff concentration derivative

A fresh Debug/GDB execution of the existing Godunov reactor-corrector fixture
locates SIGFPE in `reaction_effective_forward_rate`. At a trial temperature
2.154784396 K, reduced pressure is6.268637321e171, low-pressure rate is
3.029499544e184, high-pressure rate is5.570219452e10 and collider concentration
is0.01152589232. The coefficient is finite (5.569521552e10), but the derivative
expression squares `1+Pr`, overflowing before a representable derivative can
be computed. These are solver trial diagnostics, not a physically validated
2-K chemistry state or an endorsement of extrapolated thermodynamics.

For `Pr > 0.5*sqrt(huge)`, the implementation now uses the algebraically equal
identity `dk/dM = (k/M)*(1/(1+Pr) + d(log F)/d(log Pr))`. The ordinary-range
arithmetic is unchanged. No mechanism data, kinetic parameters, tolerance,
clipping or admission criterion changed. A synthetic Lindemann limit with
k0=1e180 and kinf=1e10 has Pr about1e170 and a representable derivative near
1e-160. The added regression checks finiteness and relative agreement with
the analytic high-pressure derivative. It fails on pre-change code and passes
after the correction; no absolute-error floor hides the small derivative.

Focused current-source checks pass the pressure-dependent kinetics and full
H2/O2 Jacobian tests (2 passes). Both original MPI reactor-corrector cases now
avoid SIGFPE but still fail `serial accepted held SDC iterate` (2 failures).
This fixes the intermediate arithmetic error; it does not claim the reactor
acceptance problem or integration tests are solved. Artifacts, GDB locals,
source snapshots, task-only patch and exact commands are under
`build/godunov-reactor-fpe-20261028/`; `verified.json` preserves that split.
Provenance, Python syntax, project contract and diff whitespace checks pass.
The next task is to localize the remaining reactor refusal without relaxing
its admission criteria or calling the unaccepted trial physical validation.


### 2026-10-28: exact forced interval for a structurally inactive mechanism

The next refusal is captured after42 adaptive substeps. The input is pure N2
at700 K, rhoY(N2)=1 kg m^-3, external N2 rate-6.777755426e7 kg m^-3 s^-1,
initial internal-energy density218595.256969 J m^-3, external energy rate
-3.182903270e13 J m^-3 s^-1, and duration1e-8 s. The existing pinned EOS
extrapolation admits the affine endpoint at approximately1.47e-11 K. These
extreme fixture values are a numerical admission control, not physical
low-temperature chemistry validation. Other species and their external rates
are exactly zero; full-mechanism kinetics unnecessarily evaluates cold rates
while chemistry is structurally absent.

`advance_forced_reactor_interval` now certifies this sufficient special case:
all reactions must pass the existing validity check, and every species appearing
in either stoichiometric direction must have exactly zero initial density and
zero external supply. Third-body-only bath species can remain present. No small
rate/density threshold is used. With no explicit initial-substep request, this
certificate selects the existing affine no-reaction interval, using the same
initial/endpoint EOS and physical admission checks. It does not remove reactions
from a mixed/active state, alter explicit substep requests, or permit negative
endpoints. Publication remains transactional.

The added unit regression compares the captured full-mechanism inert endpoint
against the empty-mechanism exact control, requires exact state/temperature
identity and one step, and checks negative-endpoint rollback. Relinking the
new test with the saved pre-change reactor fails the intended new assertion;
the revised implicit-reactor test passes. The existing forced signed-trial
suite also passes, including physical/signed Jacobian and stiff-trace checks.
Provenance, project contract and diff whitespace checks pass.

Both original MPI reactor-corrector fixtures now pass their Godunov1/2-iteration
reactor and gas/held inventory checks, then fail the MOL-side assertion
`serial reset ledger counts actual accepted forcing only`. They remain failed
integration tests; neither their old SIGFPE nor their initial reactor refusal
is reclassified as a complete pass. Inputs, isolated capture, pre-change
regression, commands, source patch and result hashes are under
`build/godunov-reactor-refusal-20261028/verified.json`. The remaining task is
to reconcile the MOL accepted-forcing/reset ledger without changing its limits.


### 2026-10-28: reconcile MOL wall-pressure ledger with actual gas pressure

Component-wise capture of the remaining MOL failure identifies x momentum,
not energy: inventory requires a wall impulse4.748534492e-8 kg m s^-1,
whereas the ledger reports3.323962425e-8. The missing1.424572067e-8 gives the
observed scaled residual0.4141197869. Mass, energy and other momentum components
close to roundoff. The fixture deliberately has approximately1000-K gas total
energy and700-K retained internal energy, exposing the pressure-policy mismatch.

The ordinary MOL hydro kernel already uses total-energy-derived pressure.
Both serial and MPI `wall_source` ledger helpers still unconditionally passed
the held scalar for MOL. They now pass it only for rotating or explicitly
selected separate-momentum MOL paths, matching the kernel's existing policy.
Godunov behavior is unchanged. All helper call sites were inspected: MOL uses
these values to integrate the wall ledger after constructing the gas rate;
the gas divergence remains computed by the existing kernel. No accepted
state, flux formula, EOS tolerance or test criterion is intentionally changed.

The failing component capture and source-only patches are in
`build/mol-wall-ledger-policy-20261028/`. The current-source Debug build and six
original MPI3 checks (standard held ledger, register closure, pinned EOS source,
temperature-seed refusal, and both reactor-corrector cases) pass 6/6 with
unchanged time limits (111.42 s total test time; 163.75 s build). Independent
verification confirms the five recorded source hashes and executable hash,
six distinct completed tests and six native PASS markers. The 20 printed MPI
inventory residuals have maximum1.048231883e-14, below the unchanged3e-12
criterion. The serial inventory checks and exact serial/MPI gas, temperature
and held-state comparisons also pass; serial residuals are asserted but not
printed separately. `verified.json` records the scope and log hashes.

These current-source results resolve the six focused integration checks,
including the previously failing reactor-corrector fixtures. Historical failed
runs above remain recorded as failures. This ledger correction does not close
fixed-PeleC trajectory, mesh/time convergence, or physical-validation gates.
Provenance, Python syntax, project contract and diff whitespace checks pass.
