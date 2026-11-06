# Completion roadmap

## Completion gates for the CPU/MPI release

The next release target is an explicitly scoped CPU/MPI alternative to the
pinned PeleC revision below. Algorithm compatibility and physical validation
are separate gates. GPU execution, new models beyond the pinned reference,
and fuels outside the selected mechanism contract require their own scope
decision; they must not silently enlarge this release or be called supported.
The following gates own the remaining work. The chronological detail below
is supporting evidence, not a list of new implementation tasks.

| Gate | Responsibility and existing starting point | Required evidence to close the gate |
| --- | --- | --- |
| C0: freeze scope and source | `references/pelec_baseline.json`, `docs/pelec_mapping.md`, current source including uncommitted/untracked numerical files | Reconcile stale mapping entries with current code; list required operators, application paths and explicit exclusions. Record source/dependency/input hashes and the observable, time, scale and tolerance for each acceptance case before comparison. |
| C1: reacting mesh integration | `cases/reactive_eb_tree_3d`, held-MOL native recovery, `check_eb_tree_reacting_restart.py` | Actual reacting H2/O2/N2 evolution through common/recursive steps, conservative regrid and changed-rank restart, independent gas/element/UEINT inventories, and transactional refusal. Nonzero-material partial recovery and prepared/routed sources still require distinct cases. |
| C2: fixed-reference trajectories | Existing Sod/HIT, `check_amr_tree_mol_reference.py`, gas transport and chemistry oracles | Matched gas/species time histories for required hydro, reaction, transport, boundary and AMR/EB combinations, plus mesh/time convergence. Preserve the original pinned solver as the authority; diagnostic variants are not replacement baselines. |
| C3: coupled particles and LES | `cases/spray_reference`, `cases/spray_jet_3d`, `cases/spray_a_reference_3d`, existing LES comparisons | Resolve gas and particle trajectory differences; progress from single-drop exchange to nonreacting breakup/AMR/walls and then reacting spray. Complete the declared duration and audit gas, liquid, wall and exterior budgets. |
| C4: research and resource envelope | Selected fuel mechanism and declared production case | Specify fuel/conditions/observables; establish physical comparison separately from parent-solver agreement. Measure runtime, peak memory, parallel scaling and output/restart cost at useful mesh/particle sizes. |
| C5: reproducible release | Existing serial/MPI, selected chemistry, optional CVODE and install workflows | Rerun the required configuration matrix and installed applications at one source identity, including rank-changing restart; report failures/skips and resolve licensing before distribution. |

C2's smallest matched comparisons accompany C1 fixes rather than waiting
until integration is finished. C3 follows its dependency cases; broad C5
qualification follows numerical stabilization. Existing implementation is
reused, not rebuilt merely because an older mapping row says "not started".


### C2 evidence routing and open comparisons

This table routes the existing cases to their authority; it is not a new,
smaller release matrix. C0 must still freeze required combinations. A checker
being implemented, an old pass, an isolated diagnostic and a current
source-bound trajectory are separate evidence classes.

| Comparison family | Contract / checker | Current boundary and required next evidence |
| --- | --- | --- |
| Regular hydro/HIT and two-level nonreacting MOL | `tools/check_hit_evolution_reference.py`, `references/pelec_amr_mol.json`, `tools/check_amr_tree_mol_reference.py` | Reuse frozen source/input hashes and all declared clock/iteration cases. The 2026-10-28 current-source serial four-case audit passes MOL1 with common/recursive clocks, but MOL2 fails both (maximum scaled difference about 3.4e-4; native mass drift about -1e-6 while PeleF conserves). See [the audit](validation/amr-mol-reference-20261019.md). An isolated native one-expression registration correction removes both differences to roundoff while preserving MOL1 outputs; the original reference remains unchanged and failed. This register/conservation conflict remains open and is separate from reacting EB. |
| Reacting regular AMR | `references/pelec_tgreact_amr.json`, `tools/check_tgreact_amr_reference.py` | Current-source serial MOL2 common/recursive clocks with imported/computed initial reaction history pass the original four-time field and inventory gates through 300 ns with transport off and with pinned Simple molecular transport; see [the dated comparison](validation/tgreact-amr-reference-20261019.md). The corrected ordinary old-state/pressure path passes the original n8 and n16 recursive transport gates through 300 ns; serial/MPI exact tests and rank-changing restart pass. The previous 13-gate n8 failure remains dated diagnostic evidence. A fresh GNU/OpenMPI n8 recursive probe passes the original reference gates for MPI2/MPI3 and imported/computed histories, with serial byte identity. Common-clock MPI3, recursive n16 and n8/dt25-ns comparisons also pass. Stationary adiabatic slip/no-slip x-wall serial/MPI controls pass original field and boundary-ledger gates; n8 no-slip covers both clocks. First-order extrapolation x boundaries pass n4 recursive serial/MPI3 field and open-inventory checks with transport off/on, plus n8 transport with common/recursive clocks. Stationary 500 K isothermal no-slip x walls pass n4/n16 recursive and n8 common/recursive serial/MPI3 original gates through 300 ns, with exact serial/MPI state/budget CSVs. Temporal self-differences decrease at n8, but n4/n8/n16 energy/temperature evolution differences increase in both solvers and peak at the wall; spatial convergence remains open. A 2050 K wall control separates the remaining interface error: first-predictor pressure sensitivity traces to density-extremum shared interpolation limiting; 96 same-input slopes match the pinned AMReX kernel exactly. This bounded shared-policy finding does not close spatial convergence or authorize changing the reference method. The subsequent n32/2050 K failure is resolved by matching total-energy transport primitives and no-slip tangential ghost reflection. Corrected n8/n16/n32 comparisons pass original trajectory gates through300 ns; n16 serial/MPI and n32 MPI/isolated serial outputs are exact. The wall-transport MOL rank-changing restart also passes with explicit covered-energy-reset checkpoint validation. Corrected z-momentum increment Linf remains nonmonotone at100/200 ns in both solvers; the matched dt50-ns three-grid study also passes reference gates but retains that nonmonotonicity. Paired stage/face-pressure and native-donor diagnostics identify shared component limiting of ghost normal velocity. An isolated independent-momentum control has decreasing all-field state/increment L1/Linf grid differences and passing conservation through300ns, but fails original reference limits and is not adopted. The shared default remains unchanged; extrapolation n4/n8/n16 passes trajectory and candidate open-ledger gates, but density/energy/temperature increment maxima remain nonmonotone at physical boundary cells. A matched zero-initial-velocity control passes reference gates and greatly reduces those maxima, but leaves nonmonotone interface density/energy/temperature maxima. Nonreacting controls and exact stage/flux replay isolate a species-transport pressure response under the held-energy reset. An initial-symmetry transport-ghost intervention nearly removes the first-predictor pressure mismatch but reduces the100-ns momentum difference by only28%; it is diagnostic and unadopted. The n32/dt50ns outflow trajectory and MPI3 comparison pass original gates with14 exact serial CSVs;8/16/32 maximum norms decrease at250/300ns but retain early-time nonmonotonicity. Earlier failures remain dated evidence. Explicit isothermal slip now passes bounded n4/500 K common/recursive serial/MPI3 reference gates and fixed/selected application/restart checks. Its n8/n16 reference gates also pass with exact serial/MPI outputs, but the500 K wall spatial audit remains nonconvergent. Convergence evidence remains open, alongside extrapolation and other thermal-wall convergence, other open prescriptions and moving walls. EB agreement cannot substitute for this family. |
| Reacting planar EB-AMR | `references/pelec_tgreact_eb_amr.json`, `tools/check_tgreact_eb_amr_reference.py`; [application evidence](validation/eb-full-cvode-application-20261024.md) | Distinguish the reference probe from the actual application, and canonical from native-activation compatibility. Serial n32/12.5-ns/200-ns evidence is bounded; n32/6.25-ns field failures and maximum-norm convergence remain open. |
| CPU/MPI trajectory equivalence | Shared EB-tree frontend, selected mechanism identity and the same native field/inventory authority | Focused n8 serial/MPI/rank-changing restart exists. The frozen activation-compatible n32 MPI2/12.5-ns/200-ns run passes all 17 original native field comparisons, serial byte identity, clocks and inventory; independently rechecked in `build/pelef-mpi2-time12-activation200-ordinary-20261028/verified.json`. This is not a fresh current-source release matrix, and the short-step/convergence gates remain open. |
| Temporal/spatial convergence | `tools/analyze_eb_temporal_convergence.py`, `tools/compare_tgreact_eb_meshes.py`, frozen parent references | Declining weighted norms and a localized shared-limiter sensitivity do not close nonmonotone maximum norms or failed original trajectory limits. Diagnostic native variants never replace the pinned authority. |

The newest pressure/source/halo/flux/chemistry controls are maintained only in
the dated application record. Do not repeat arbitrary tolerance/rounding
sweeps solely to obtain a favorable member of the n32 family. Progress on
independent C2 obligations does not remove that family's unresolved gates.

The regular serial/MPI application now independently initializes fresh MOL
reaction history with level clocks and preserves restored checkpoint history.
The computed-history TGReact comparison and selected LiDryer continuation pass
([initialization record](validation/mol-initial-reactions-20261019.md)).
The gas EB startup connection also passes fluid-mask/rollback and recursive
reacting regrid/changed-rank continuation checks
([EB initialization](validation/eb-initial-reactions-20261019.md)).
Remaining C1/C2 work includes matched spatial reaction/transport evolution
with refinement/convergence; the nonreacting MOL2
register/conservation difference remains open.
The ordinary nonrotating EB-AMR comparison now passes the original pinned
300 ns field/budget gates after consistent total-energy pressure recovery for
faces and walls. Focused serial/MPI/restart tests pass; fixed-reference
mesh/time convergence remains open
([current evidence](validation/eb-default-localization-20261023.md)).
The earlier diagnostic progression below records why intermediate corrections
alone did not close those field gates.
Ordinary held serial/MPI closure now connects the missing flux re-reflux
redistribution even with NoRedist. At 300 ns this reduces density error by
about 72 times and x-momentum error by about 10 times, but the field gate
still fails ([re-reflux correction](validation/eb-no-redist-rereflux-20261019.md)).
Connected held-MOL coarse/fine halo interpolation now matches initial native
conserved values to roundoff in the diagnostic; endpoint maxima are nearly
unchanged ([halo correction](validation/eb-mol-linear-halo-20261019.md)).
Ordinary held transport now uses the full spatially interpolated neighborhood,
with serial/MPI restart checks passing, but endpoint field gates remain open.
Primitive energy substitution alone did not improve the coupled trajectory.
Transport centroid reconstruction with conservative species/enthalpy limiting
and two-sided interface viscous gradients are now implemented and checked.
Fine-level total/held-energy reset timing after parent re-reflux is now
corrected in ordinary reacting NoRedist MOL, with actual held EOS and repair
fallback preserved. First-step field limits pass; the 300 ns temperature
error falls to `1.10e-4 K`, but twelve later density/momentum/energy limits
still fail ([evidence](validation/eb-energy-sync-20261020.md)). Importing the
native 100 ns state and reaction history reproduces nearly the same next-step
error; accumulated startup/history error is not its main cause. Covered coarse donor
energy is now repaired after restriction and retained through driver entry/
exit; its gap is at roundoff, but the twelve field failures remain nearly
unchanged ([record](validation/eb-covered-energy-20261020.md)). A native diagnostic
using UEINT consistently for spatial pressure and transport temperature now
satisfies every original field limit through 300 ns. This identifies the
dominant remaining gap as the spatial energy-policy difference
([record](validation/eb-held-policy-20261020.md)); it does not close the
unmodified pinned-reference gate. Keep that gate open and preserve the
held-energy cancellation contract rather than making further flux/reset
changes merely to erase this policy difference. A completed 100/50/25 ns
study finds observed orders near 1.1--1.2 in both implementations; the 1.2
microsecond trajectory retains conservation but fails O2/H2O field limits
from 1.0 microsecond ([record](validation/eb-temporal-study-20261020.md)).
The native activation-energy conversion is now isolated: an explicitly
generated compatibility mechanism removes all six longer O2/H2O exceedances
against held-EOS native through 1.2 microseconds, without altering the
canonical mechanism ([record](validation/eb-native-activation-20261020.md)).
The variable-grid 16-cubed pilot passes through 100 ns with 50 ns root steps
([record](validation/eb-mesh-pilot-20261021.md)). The common-50-ns 8/16 pair and
conservative fluid-volume restriction now separate initial sampling from
evolution differences ([record](validation/eb-mesh-restriction-20261021.md)).
The 8/16/32 diagnostic now passes all field and budget gates through 100 ns
with explicit tighter chemistry and caloric inverse controls. Forced-CVODE
admission and temperature-residual amplification have focused corrections
([EOS record](validation/eb-eos-convergence-20261021.md),
[forced-reactor record](validation/eb-forced-cvode-20261021.md)).
Observed evolution orders are not uniformly second order. The fixed 33-point
profile explains its first-order initial-temperature differences. The separate
129-point cosine-profile study passes all three grids and restores near-second
order endpoint temperature/density/H2O, but gas evolution-increment orders
remain about 1.1--1.5 in both solvers. Regional saved-field analysis does not
localize the discrepancy solely to the wall. Separate temporal, limiter and
refinement-interface effects before assigning the remaining order loss,
keeping the frozen-profile result intact
([three-grid record](validation/eb-three-grid-20261021.md)).
Fixed 8/16 meshes now also pass 50/25/12.5 ns time refinement. Observed
temporal orders remain near one in both solvers, while quartering the step
changes density/energy/temperature spatial differences by less than 1%.
Initial-flow, limiter and interface effects therefore remain the next
spatial diagnostics ([time record](validation/eb-smooth-time-refinement-20261021.md)).
The initially stationary control passes 8/16 but exposes 19 field-limit
failures at 32. Saved coarse stencils show shared-limiter sensitivity to
near-zero transverse momentum.
Resolve this case before claiming robust spatial agreement
([stationary record](validation/eb-stationary-control-20261021.md)).
Live-stage intervention and a matched diagnostic ablation now identify
shared momentum limiting as the stationary first-step amplification mechanism.
An opt-in [bounded momentum kernel](design_decisions/eb-halo-momentum-limiter.md)
now reaches common/recursive held MOL. Kernel/operator, serial/MPI retry and
full-history replay checks pass, including policy mismatch refusal. Application
input/checkpoint dispatch and fixed/selected common/recursive lifecycle checks
pass, including changed-rank restart and policy-mismatch refusal
([application evidence](validation/eb-momentum-application-20261021.md)).
The stationary matched-candidate 8/16/32 field comparison passes, but
evolved temperature/energy orders are only 0.26/0.46; endpoint near-second
order is dominated by initial restriction differences. Most evolved error
lies near coarse-fine interfaces and cut cells
([spatial diagnostic](validation/eb-momentum-stationary-mesh-20261021.md)).
Candidate convergence and adoption remain open. Halving the step changes
8/16 evolved temperature/energy differences by less than1%; further blind
grid sweeps are not the next action. The next focused C2 discriminator is
an initial spatial-operator audit: capture coarse-fine ghost pressure and
normal mass/momentum/energy fluxes before chemistry, correlate them with the
covered-coarse restriction pressure difference, and compare periodic interface
and cut-cell controls. Preserve conservative restriction and the original
reference path; change numerical policy only after this isolates a cause.
The initial production-operator audit is now implemented and the eight
8-grid controls complete. Nonuniform interfaces have a 0.275 Pa halo pressure
offset before chemistry, with mass flux strongly correlated with pressure
jump; the disturbance persists without cut cells and vanishes for uniform
states. Initial shared/separate momentum results coincide. Covered-coarse
held-EOS pressure and stored-temperature pressure differ, so the next step is
a conservative restriction/interpolation stencil audit using the actual EOS,
not a pressure reset ([initial audit](validation/eb-initial-interface-20261021.md)).
The independent stencil reconstruction now reproduces all adjacent gas/held
halos exactly, with child-mean and fluid-volume restriction errors below
2.5e-16. The worst interpolated EOS pressure falls below all donor/neighbor
pressures; EOS-consistent temperature guesses do not change it. A regular
control with recomputed conservative restriction retains the pressure deficit.
The next numerical candidate must therefore address thermodynamic pressure
reconstruction while retaining conservation and spatial accuracy, rather than
repairing donor indexing ([stencil audit](validation/eb-initial-stencil-20261021.md)).
A common conservative pressure-envelope slope contraction was also tested:
it suppresses 48/128 planar linear stencils nearly to PCM (minimum factor
1.81e-8). It is not adopted on the strength of pressure boundedness alone;
retain resolved spatial gradients in the next thermodynamic reconstruction.
A paired spatial-pressure diagnostic now retains interpolated density,
momentum and species, but derives covered coarse pressure from volume-averaged
fine EOS pressure and reconstructs auxiliary energies on both sides of the
interface. The planar initial energy-flux mismatch falls from64.98 to9.52e-8
W/m2; fine-only pressure repair instead worsens it to203.57 W/m2.
This is an operator-local policy experiment, not conservative energy
prolongation or application adoption. Qualify time evolution, reflux and
held-energy work before threading it through serial/MPI applications
([pressure representation](validation/eb-pressure-ghost-20261021.md)).
The paired candidate now also passes a private 50 ns forward-Euler hydro
increment with production gas/scalar reflux on stationary and three moving
controls. Composite mass/energy/species errors stay below6.5e-15 and held
energy below2.3e-15. Stationary pressure remains at roundoff; moving pressure
offsets improve fromabout0.0196 Pa to0.000327/0.000104 Pa. This advances the
one-step conservation discriminator only; repeated staged hydro and actual
held-MOL/reacting closure remain next
([one-step record](validation/eb-pressure-step-20261021.md)).
The paired diagnostic now reevaluates both stages of synchronous SSPRK2.
Stationary equilibrium holds through200 ns and moving pressure error is
0.001247 Pa versus0.061235 Pa originally, with closed composite inventories.
However, local E-H-KE gaps accumulate and largely cancel in the global sum;
the completed50/25 ns equal-endpoint comparison changes pressure by0.0053%
and local energy gap by-0.56%, so time-step reduction does not remove it. Route this through actual
held closure with explicit numerical-energy accounting before application
adoption ([staged diagnostic](validation/eb-pressure-evolution-20261022.md)).
The existing cell dual-energy reset is now tested after accepted diagnostic
RK steps with signed and absolute correction ledgers. It closes the local
E-H-KE identity and passes correction-accounted budgets, while moving pressure
error stays0.001247 Pa. Residual pressure error is therefore not cured by
reset alone. Full prepared-MOL/re-reflux and spatial-accuracy qualification
remain open ([reset diagnostic](validation/eb-pressure-reset-20261022.md)).
Cell-history localization now places the paired moving residual on the coarse
level, with nearly equal cut/full-cell maxima. A regular geometry ablation
retains the error; coarse advection versus interface synchronization on a
consistent smooth initializer is the next discriminator
([localization evidence](validation/eb-pressure-location-20261022.md)).
A smooth 900 K H2/N2 contact now passes initial restriction to4.41e-16.
It preserves pressure through the first Euler stage but develops a fine
pressure error during the second spatial update, followed by a coarse
interface error during synchronization. Initial paired pressure alone is
therefore insufficient: audit stage-updated material/held flux consistency
before production adoption
([stage synchronization record](validation/eb-pressure-sync-20261022.md)).
The paired pressure representation is now wired into ordinary serial/MPI
common and recursive held MOL with explicit policy/restart admission
([application integration](validation/eb-pressure-application-20261022.md)).
Time-consistent recursive ancestor history is implemented; see the
[recursive/native MOL evidence](validation/eb-pressure-recursive-20261022.md).
The bounded nonisothermal pressure-refinement check is now complete:
8/16/32 pressure errors decrease from 0.001271 to 0.0001390 to 0.00001456 Pa
with closed correction-accounted budgets
([mesh evidence](validation/eb-pressure-mesh-20261022.md)).
The reacting reference comparison is now executed and fails the unchanged
field gates; paired temperature difference is 0.239 K versus 0.000110 K
for the default. The four-policy ablation places the dominant departure at separate momentum
(0.241 K), with consistent velocity and paired pressure slightly reducing the
temperature difference. Assess variable-velocity accuracy of that limiter
before choosing a numerical correction
([reacting evidence](validation/eb-pressure-reacting-20261022.md)).
Recorded-stencil replay now confirms the limiter formulas to 1.12e-16;
48/384 scalar factors change when resolved momentum is excluded, while
initial coarse/fine flux mismatch decreases
([attribution](validation/eb-momentum-attribution-20261022.md)).
Keep original-policy parity and candidate accuracy as separate gates. The
analytic shear 8/16/32 four-policy comparison is now complete: most integrated
errors improve with momentum separation, but z-momentum L2 and transverse
Linf behavior leave all-field convergence open
([analytic results](validation/eb-shear-analytic-20261022.md)).
The 16-grid localization/time-step/regular comparison is complete: maxima
remain coarse/fine-adjacent, paired cut-case maxima change less than 0.3% on
time halving, and regular geometry reduces candidate maxima without removing
them. Next separate coarse/fine reconstruction/flux effects from wall coupling
using operator-local evidence; the regular ablation also changes fluid volume
and wall position. The completed initial operator audit finds zero transverse
rates and roundoff-level halo pressure errors for all eight cases, but nonzero
streamwise/density/energy rates. The first-step trace now brackets resolved pressure error between the first
fine spatial exit and the next predictor spatial entry, before closure/reflux.
The sampled cellwise predictor balance now closes to roundoff and attributes
99.4–99.9% of representative pressure drops to density change before recover;
the reset leaves rho/H/pressure unchanged. The saved six-face balance reproduces the density change entirely through
the x coarse/fine face; cut ghosts retain the parent mean because the pinned
3x3x3 connectivity guard disables tangential interpolation near the wall.
The isolated connected-direction candidate completed one root step: all
nonzero cut-field norms improve with momentum separation, regular outputs
match and budgets pass, but default y-momentum L2 worsens 5.47%. The completed 200 ns 8/16/32 comparison confirms decreasing separated-policy
L1/L2 errors, but y-momentum Linf is nonmonotone and z-momentum Linf convergence
is weak. The optional planner kernel passes its focused test. Application dispatch and lifecycle/policy identity are implemented behind
a false-default selector requiring separate momentum; all nine focused integration
checks pass, as recorded in the [application record](validation/eb-connected-application-20261022.md). The planner now rejects disconnected retained diagonals; sphere
and corner donor-enumeration tests pass, but general geometry flow accuracy
and volume-weighted reconstruction conservation remain unqualified;
The production serial recursive 16-cell/50 ns cut/regular comparison matches
all 32 diagnostic tables against the candidate (enabled) or original (default),
with budgets below 6.06e-14 ([evidence](validation/eb-connected-shear-integration-20261023.md));
local-error convergence remains an explicit open gate. The saved maxima
separate cut-plane y errors from full-volume x/z-interface z errors; next
compare directional pressure/flux differences at those two locations before
changing the limiter ([audit](validation/eb-shear-local-convergence-20261023.md)). Do not tune the limiter solely to match the pinned policy.
The smooth-contact defect is now traced to independently limited density and
momentum: fine halos violate constant velocity by0.0888 m/s. A further opt-in
interpolator couples advected momentum to the limited density slope and limits
relative momentum. In the paired diagnostic it reduces the200 ns pressure
offset from0.0015503 Pa to2.91e-11 Pa, with conserved child averages and closed
composite budgets. Existing policy histories remain unchanged. Carry this
candidate through full held-MOL and nonuniform-temperature/cut controls before
application exposure ([velocity correction](validation/eb-contact-velocity-20261022.md)).
Planar cut-cell controls now retain the smooth isothermal contact to2.91e-11 Pa
with the velocity candidate. The imported nonuniform-temperature control
remains at0.001247 Pa with or without it, despite closed budgets. Separate
thermal reconstruction/mixing error from the solved constant-velocity defect
before treating the candidate as a general interface correction
([cut-cell controls](validation/eb-velocity-cut-20261022.md)).
A flux-free local EOS discriminator now shows that conservative mixing of
nonisothermal donors alone produces0.001004 Pa pressure departure; equalizing
composition retains it, while equalizing temperature reduces it to roundoff.
This is not a decomposition of the MOL trajectory and does not close its field
gate. Preserve the explicit energy budget when qualifying thermal pressure
behavior; do not equate all nonzero pressure error with faulty synchronization
([mixing discriminator](validation/eb-thermal-mixing-20261022.md)).
The velocity-consistent candidate is now threaded through ordinary serial/MPI
held-MOL hydro/transport halos, common/recursive stages, retry and replay
policy identity, and application restart metadata. Fixed/selected common and recursive lifecycle, selected hot reacting,
regrid and2-to-3-rank restart checks pass. Pressure reconstruction remains
separate and diagnostic-only; physical/reference pressure convergence is
still open ([application integration](validation/eb-velocity-application-20261022.md)).
Default reference interpolation is unchanged.
The opt-in caloric mode is also exposed in the ordinary EB-tree namelist;
selected-mechanism serial/MPI regrid/restart and policy-mismatch refusal pass
([application record](validation/eb-eos-application-20261021.md)).
Establish asymptotic time/spatial
convergence with explicit energy and mechanism conventions; the current time
trend is not second order. Keep the original-reference result separately
visible. Extend inlet and
missing same-level spatial stencils beyond their current PCM fallback separately.

The first new application gate uses a reactive H2/O2/N2 mixture rather than
the pure-N2 chemistry-path control. It must demonstrate a nonzero chemical
change and a reaction-off intervention, not just enable a chemistry flag.
Its small mesh and short duration do not qualify ignition, production AMR,
native partial-source recovery, or fixed-PeleC trajectories.

For each numerical increment, retain one reproducible failure or missing
capability, the smallest necessary change, and its acceptance result. Report
which gate advanced. A new low-level recovery branch must identify the
required application case it enables. Conservation-preserving deviations
from upstream cleanup need their own activation/convergence evidence; do not
silently treat them as exact algorithm parity or weaken their checks.

Long jobs retain a process identity, immutable input/source identity, logs
and a completion/resume condition. Work on independent checks while they
run; when only waiting remains, hand back the job record without repeated
polling. Elapsed days do not reset the completion plan.

## Definition and frozen reference

A selected PeleC responsibility is complete only when it has an independent
Fortran implementation and automated numerical evidence, or is explicitly
excluded with a reason and a usable alternative. A scaffold, test count or
clean compilation is not sufficient.

The reference is PeleC `development`
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; recursive dependencies are pinned in
[the manifest](../references/pelec_baseline.json). Changing that baseline must
identify affected numerical signatures and acceptance criteria.

Ordinary held MOL has [full-root replay evidence](validation/native-mol-root-replay-20261011.md),
[joint gas/UEINT spatial-weight evidence](validation/native-mol-weighted-root-replay-20261011.md)
and [sub-parent replay evidence](validation/native-mol-subparent-replay-20261012.md).
Explicit observed windows now restore ancestor clocks/endpoints and open
parent families. Native current-stock admission now resolves the observed
future-fine root refusal and qualifies actual nonzero gas/UEINT receipt payment;
see the [funding record](validation/native-mol-future-source-funding-20261012.md).
Explicit frozen root/first-fine replay with nonzero material pending receipts
now has [full-mechanism evidence](validation/native-mol-material-replay-20261013.md).
Automatic weighted recovery of completely observed ordinary held MOL parent
windows now has [natural closure evidence](validation/native-mol-automatic-parent-20261014.md).
Selected common-clock [observed-source capture](validation/native-mol-source-prefix-20261014.md)
now retains the actual gas/UEINT packets without future entries. Explicit
[common-clock consumption](validation/native-mol-prefix-consumer-20261015.md)
rebuilds spatial operators and evaluates newly reached stages once. Recursive
common-clock retries now have [retry and rollback evidence](validation/native-mol-recursive-prefix-retry-20261015.md).
Successful partial single-root recovery now has
[serial/MPI inventory and rollback evidence](validation/native-mol-partial-root-recovery-20261016.md).
Partial root recovery now completes descendants with
[composite inventories and late-child rollback](validation/native-mol-partial-descendants-20261016.md).
Partial leaf fine-clock recovery under an open root now has
[mixed-leaf serial/MPI and explicit-prefix evidence](validation/native-mol-partial-fine-20261016.md).
Partial sub-parent clocks with descendants and continuation from a
complete-clock prefix now have [focused evidence](validation/native-mol-window-prefix-20261018.md).
The next gates are broader partial-window replacement and accepted partial
recovery with nonzero pending material funding. Automatic nonzero-material
sub-parent recovery now has [controlled three-level evidence](validation/native-mol-automatic-material-20261014.md).
Prepared/routed consumers, broader material debt topologies and pinned solver
trajectories remain separate completion gates.

The pinned reference rejects non-regular EB cells for Smagorinsky, Dynamic,
WALE and Vreman LES. Cut-cell LES is therefore a future model extension, not a
missing reference operator; the implemented regular-grid AMR LES path remains
the usable alternative. Independent serial/MPI SSPRK2 LES level stepping,
recursive intervals and flux synchronization now reach the MPI spray application
([coupling evidence](validation/les-mpi-spray-subcycle-20260925.md)). This does not
establish parity with the reference time integrator or physical validation.
See the [dated source audit](validation/les-reference-scope-20260925.md).

Watson latent heat and the PeleMP pressure-dependent boiling calculation now
have independent library APIs with distinct reference and boiling temperatures
([kernel evidence](validation/watson-properties-20260925.md)). The optional
`pelemp-cap` mode now connects distinct liquid inputs, the vapor-state cap and
restart identity, including ordinary EOS enthalpy-derived saturation latent heat
([integration evidence](validation/pelemp-cap-integration-20260925.md)).
The constant-latent `surface-cap` remains available and `reject` remains the
default. The ordinary EOS vapor-pressure path uses
enthalpy-derived latent heat, so replacing every latent heat with Watson would
not match the reference. Isolated property and fixed/selected regular
conservation/restart checks are complete. EB tree input and changed-rank
continuation now also have [focused evidence](validation/pelemp-eb-tree-input-20260925.md).
MPI AMR execution of this option, broader transfer-integrator parity and physical
qualification remain separate. Shared gas-species mapping now accumulates
component vapor and sources conservatively
([mapping evidence](validation/shared-vapor-mapping-20260925.md)). Independent
liquid molecular-weight inputs now feed liquid mole fractions, vapor pressure,
boiling and restart identity
([integration evidence](validation/liquid-molecular-weights-20260925.md)).
Single/mixed vapor-state conversion now uses the complete ambient mean molecular
weight and reference activity ordering
([skin-state evidence](validation/ambient-vapor-skin-20260925.md)).
GCM is excluded from this pinned PeleC port because upstream restricts it to
PeleLMeX; PeleMP is the available alternative. See the
[property source audit](validation/spray-property-scope-20260925.md).

## Rotating-frame integration

The momentum source, reversible relative-energy conversion and MOL interface
energy corrections now have [focused evidence](validation/rotating-frame-kernels-20260927.md).
Cell EOS, wall-pressure recovery and single-patch flux divergence are connected
([bounded patch evidence and larger-speed failure](validation/rotating-mol-patch-20260927.md)).
Serial tree/halo recovery, conservative PCM/restriction and inlet energy adapters
are connected ([transfer evidence](validation/rotating-eb-tree-transfer-20260927.md)).
Serial inviscid, nonreacting old/new sources, time integration and reflux recovery
are connected ([SDC evidence](validation/rotating-sdc-20260927.md)).
Forced chemistry now carries rotating energy and density forcing
([reaction evidence](validation/rotating-reaction-20260927.md)).
Serial transport and embedded-wall energy adapters are connected
([transport evidence](validation/rotating-transport-20260927.md)).
MPI owner-local transfer/stepping, conservative remapping, timestep selection,
boundary/initial states, budgets and checkpoint identity now carry the frame in
the serial and MPI gas applications. Frame axis, origin and angular speed must
agree across restarts. MPI EB spray now converts particle/film sampling and
deposition, LES states, inventories, injection and outflow to rotating energy.
The pinned [particle setup](https://github.com/Pele-Suite/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/Particle.cpp)
passes constant external forcing to the pinned
[droplet update](https://github.com/Pele-Suite/PelePhysics/blob/6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b/Source/Spray/Drag.H);
the inspected path has no position-dependent rotating particle force. Compare
particle trajectories and complete boundary/force ledgers against that fixed
behavior before deciding whether any extension is needed. Source/flux
kernel agreement alone does not satisfy this responsibility.

## Consolidation versus solver completion

The source-provenance repair, documentation ownership, CMake module split,
archived PR11 supersession review, and hosted serial/MPI/CVODE/clean-install
workflows are implemented. [PR #102](https://github.com/GTT10/Fortranslate/pull/102)
records the actual integration status and completed run identities; see
[current status](current_status.md) and [the layout evidence](validation/build-layout-20260909.md).
A targeted NASA7 failure-status fix and its red/green tests are recorded
[separately](validation/nasa7-conversion-status-20260909.md).

Do not redo this housekeeping as a substitute for the remaining physical
capabilities. Conversely, do not label the following capabilities implemented
merely because the repository has been consolidated.

## 0.246.0 functional increment

The opt-in regular 3D path now connects coexisting single-component liquid populations selected per port, conservative
mass/momentum/total-energy exchange, drag/heating/evaporation, timed cone
injection from point or finite circular solid/hollow outlets with swirl with up to 16 independently configured ports, shared or independent constant/tabulated-rate schedules with prescribed or discharge-based speeds and optional per-port mass- or number-CDF Weibull/normal/lognormal/ChiSquared diameters, deviatoric Smagorinsky and coupled checkpoint/restart to the existing
gas solver. Cell-local CVODE is connected to this CFD path, with a pinned
FFCM1_Red liquid-methanol burn fixture and independent Cantera gas-only oracles.
Chemical labels, single-interval NASA7 and classified Lindemann/SRI/Tsang, chemically activated, PLOG and Chebyshev ingestion are
supported. Selected mechanism capacity is now 128 species; a 53-species GRI-Mech
0D and coupled-spray case has focused evidence. Large-mechanism MPI/scaling and
production diesel chemistry qualification remain open. Lu n-dodecane (53 species,
268 reactions) now imports with extended Troe coefficients and completes a
lean constant-volume ignition comparison; Newton trial-state handling preserves
accepted-state checks. See [focused evidence](validation/dodecane-20260915.md). See [scope and commands](spray.md) and [acceptance gates](validation/0.246.0.md).
These are bounded implemented capabilities. Multicomponent liquids and MPI/AMR
particle coupling are implemented; detailed-fuel physical qualification,
mixture boiling, general dispersion/collision and near-wall turbulence
qualification remain open.
Regular wall reflection and outflow now include gas/particle surface budgets
and restart persistence; isothermal/tangentially moving gas walls include molecular
and SGS heat/work budgets. The Ahamed impact-regime, splash-velocity, mass-distribution, four-child construction and deposited-film geometry kernels are connected by a transactional impact event with explicit wall impulse/energy budgets; regular coupled-step boundary handling now retains deposition patches, advances splash children and evolves stationary film temperature/mass with per-material constant or cubic conductivity; the CLI and coupled checkpoint retain these states, while physical qualification remains open. Lateral film transport is a future model extension beyond the pinned stationary-film path, not an omitted reference transport operator ([source audit](validation/film-reference-scope-20260918.md)). ETAB and KH-RT representative-drop
size/deformation, bounded computational parcel splitting and restart are implemented;
KH shed-mass accumulation and conservative four-child transverse release now
connect to CFD and restart, including computational splitting with partitioned
KH inventories and stationary wall impact processed before KH release on a
crossing substep. General turbulent dispersion, wet-impact correlation limits
and physical qualification remain open. Uniform gas/particle acceleration
now retains external impulse/work separately through coupled restart. Single-component
Antoine saturation-pressure fits are connected with explicit temperature bounds
and restart identity. Polynomial liquid density now drives injection, drag,
evaporation and breakup consistently; inverse-temperature liquid viscosity now
drives ETAB/KH-RT. The pinned piecewise drag law and ETAB deformation multiplier are selectable alongside Schiller-Naumann drag, with restart identity. An optional near-boiling vapor-state temperature cap follows
the pinned reference convention. Port material identity now routes evaporation to each liquid's gas species and persists through restart; density/pressure fits, vapor-state caps and breakup viscosity/surface tension can be set independently per port. Polynomial liquid heat capacity now supplies the thermal update and integrated liquid energy, including per-port fits. Mixed-liquid density, heat capacity, stored energy and fixed-composition heating
now have a shared thermodynamic API. Raoult surface partitioning and supplied
component-mass evaporation exchange, diffusion-weighted rates and a conservative
frozen-gas transport update are implemented. Regular-grid feedback, persistent
component masses, boundary budgets, representative breakup/splitting and composition
checkpoint APIs are connected. Combined hydro/transport/LES/chemistry stepping
now retains mixture composition through accepted updates and rollback. Physical mixture boiling rates remain open.
The mixture impact event kernel now partitions component masses into retained
parent, splash and film with wall budgets. Supplied film evaporation/heat exchange
now updates composition, temperature and height. Mixed film heat/evaporation rates
and a frozen-rate update are connected to grid feedback. Planar
crossing dispatch now computes mixed impact parameters and retains component budgets.
The checkpoint API retains mixed-film component masses. Stationary mixed-film grid
feedback now updates gas and film with wall budgets. Particle crossing, splash and
remaining-interval film evolution now share one transaction. Reactive stepping
retains both component matrices with full rollback. CLI wall inputs, component-film
output and restart are now connected. Shared-composition CLI inputs, initial/scheduled injection,
component output and restart are now connected. Mixed KH shed children preserve
composition and conservative transverse energy through coupling, splitting and restart.
Component density/heat-capacity fits, Antoine pressure and Clausius surface caps
are connected to CLI inputs and restart identity. Per-port compositions now feed
initial/scheduled injection, component budgets and density-based injection speeds.

## Ordered implementation work

EB hydro currently uses cell-centered Riemann states. The pinned PLM covered-cell
slope policy (4th/2nd/zero fallback) is now available through the normal PLM
reconstruction API ([checks](validation/eb-plm-slope-20260925.md)). An aperture-aware
three-direction normal predictor now returns primitive states and independently
traced internal energy ([checks](validation/eb-plm-normal-20260925.md)). A three-layer
halo API now constructs normal Riemann fluxes and Godunov pressure/velocity
([checks](validation/eb-plm-normal-flux-20260925.md)). Cell-local single/double
transverse correction now includes the joint aperture gate, pressure work and
primitive sources ([checks](validation/eb-transverse-20260925.md)). The twelve-solve
CTU flux pipeline is connected, including corrected Godunov pressure/velocity
([checks](validation/eb-ctu-flux-20260925.md)). A connectivity-gated face-centroid
flux interpolation API now follows the reference regular-area/aperture ordering
([checks](validation/eb-centroid-flux-20260925.md)). The prepared-halo hydro API
now assembles CTU, viscosity/normalization, centroid fluxes and conservative,
FluxRedist or StateRedist updates with integrated flux and inventory publication
([checks](validation/eb-ctu-update-20260925.md)). A cut-aware, multidirectional
flattening preparation API now feeds the prepared-halo update
([checks](validation/eb-flattening-20260925.md)). Physical wall/outflow Riemann
corrections and interior physical-wall viscosity masks now reach the complete
prepared CTU update ([checks](validation/eb-physical-flux-20260925.md)).
Configured plane/extruded/nodal geometry now supports domain-face extension,
and an outflow entry prepares state halos and flattening before the CTU update
([checks](validation/eb-extended-geometry-20260925.md)). Fixed/selected single-level applications now select `plm_ctu`, compose it with
reaction/transport and bind checkpoint continuation to the hydro method
([checks](validation/eb-ctu-application-20260925.md)). Serial/MPI patch-tree
halo APIs now accept six-layer target-resolution metrics with peer/ancestor
exchange and temporal interpolation ([checks](validation/eb-wide-halo-20260925.md)).
A conservative extended-geometry builder and serial common-interval CTU
patch stepping now connect peer synchronization and coarse/fine reflux
([checks](validation/eb-tree-ctu-20260925.md)). Serial recursive CTU stepping
and coupled reaction/transport composition now also connect time-interpolated
ancestor halos and accumulated fine fluxes
([checks](validation/eb-tree-ctu-subcycle-20260925.md)). MPI common CTU steps
now match serial under distributed and idle-rank ownership
([checks](validation/mpi-eb-ctu-common-20260925.md)). Recursive MPI CTU and
coupled common/recursive dispatch now also match serial
([checks](validation/mpi-eb-ctu-recursive-20260925.md)). Fixed/selected serial/MPI
tree application selection, scheduled regridding and changed-rank restart are
connected ([checks](validation/eb-tree-ctu-application-20260925.md)). Other physical
boundary halo preparation remains before claiming EB Godunov parity. MPI EB
spray now forwards CTU through coupled gas stepping and checkpoint identity
([checks](validation/mpi-eb-spray-ctu-20260925.md)). The serial single-level spray
entry and fixed/selected restart now also connect CTU
([checks](validation/serial-eb-spray-ctu-20260925.md)); broader physical-boundary
preparation and complete reference parity remain open. A single-level library
entry now prepares six-layer periodic/fixed-wall/prescribed-inflow halos for CTU
([checks](validation/eb-ctu-boundary-halos-20260925.md)); coupled transport,
application and tree/MPI boundary selection still need connection. The redistributed
transport Euler stage now accepts prepared donors and returns transactional
integrated face fluxes ([checks](validation/eb-transport-redistributed-halos-20260925.md));
SSPRK2 now invokes a caller-supplied donor provider at both stages and returns
weighted integrated fluxes ([checks](validation/eb-transport-stage-halos-20260925.md)).
Single-level coupled CTU now forwards inflow/outflow labels and a transport donor
provider with combined hydro/transport boundary accounting
([checks](validation/eb-coupled-inflow-boundaries-20260925.md)). Built-in inflow/outflow donors now reuse the regular boundary filler at every
transport stage ([checks](validation/eb-built-in-inflow-donors-20260925.md)).
Fixed/selected single-level gas applications now connect inflow/outflow input,
open-face donor timestep bounds and boundary-bound restart
([checks](validation/eb-inflow-application-20260925.md)). Fixed-wall transport policies
and isothermal heat flux now connect to single-level coupled CTU
([checks](validation/eb-fixed-wall-transport-20260925.md)). An aperture/fluid-volume
weighted isothermal timestep bound is available
([checks](validation/eb-isothermal-timestep-20260925.md)). Reflected configured
geometry and fixed/selected single-level wall input, thermal timestep and restart
binding are connected ([checks](validation/eb-wall-application-20260925.md)).
The periodic transport operator now wraps gradients and shares aperture-weighted
seam fluxes and species limiter coefficients
([checks](validation/eb-periodic-transport-20260925.md)). Periodic configured geometry,
coupled gas dispatch and fixed/selected application restart are now connected
([checks](validation/eb-periodic-application-20260925.md)). Fixed/selected single-level gas applications now accept independent inlet mixtures;
single-level spray state generation and restart now forward those mixtures. Serial single-level spray boundary
input now covers inflow/outflow, periodic and domain-wall gas conditions with
specular particle reflection ([checks](validation/eb-spray-wall-20260926.md)).
The serial EB tree halo API and common/recursive time updates now accept
prescribed inflow states. MPI gas stepping and fixed/selected application
dispatch now include changed-rank restart
([checks](validation/eb-tree-inflow-application-20260926.md)). Tree applications
also accept independent named mixtures per inlet face
([checks](validation/eb-tree-inlet-composition-20260926.md)). Periodic gas and MPI
spray boundaries are connected; adiabatic tree gas wall inputs and restart are connected
([checks](validation/eb-tree-wall-application-20260926.md)). Isothermal tree gas and MPI spray input, stepping and restart are connected
([checks](validation/mpi-eb-isothermal-spray-20260926.md)); MPI spray specular wall routing and changed-rank application
restart are connected
([checks](validation/mpi-eb-wall-spray-application-20260926.md)).

The 3D reactive flow path still lacks complete characteristic PPM CTU tracing.
A periodic three-axis normal PPM face-flux predictor now reuses the existing
1D/2D characteristic reconstruction, with analytic-flux, rotational and periodic
conservation [checks](validation/ppm-normal-3d-20260925.md). Twelve-solve
transverse/corner flux construction is now implemented with PeleF conservative
EOS recovery ([focused checks](validation/ppm-ctu-3d-20260925.md)); reference
primitive-pressure parity remains open. The periodic conservative CTU time
update and full-split library dispatch now have multistep conservation and
rollback [checks](validation/ppm-step-3d-20260925.md). Fixed/selected regular
gas application input now selects periodic PPM
([input checks](validation/ppm-input-3d-20260925.md)). Four-layer CTU halos and
physical-boundary library stepping now have exact subdomain-flux and boundary
ledger [checks](validation/ppm-halo-3d-20260925.md). Fixed/selected regular gas
applications now also select physical-boundary PPM with multistep molecular
transport [checks](validation/ppm-boundary-input-20260925.md). Serial common-step
AMR tree hydro now connects PPM and reflux
([checks](validation/ppm-amr-common-20260925.md)). Serial recursive level
stepping also supports PPM with integrated fluxes and 1/2/6-step synchronization
([checks](validation/ppm-amr-recursive-20260925.md)). MPI common-step PPM now
matches serial fields and diagnostics, including empty ranks
([checks](validation/ppm-mpi-amr-common-20260925.md)). Recursive MPI PPM now
also matches serial fields, diagnostics and 1/2/6 step counts
([checks](validation/ppm-mpi-amr-recursive-20260925.md)). Fixed/selected
serial/MPI AMR tree gas application input now enables PPM, with reacting
transport and fixed changed-rank continuation
[checks](validation/ppm-amr-application-20260925.md). MPI AMR spray input now also enables PPM,
with fixed 2-to-4-rank particle/gas continuation
[checks](validation/ppm-amr-spray-20260925.md). Regular serial spray also exposes PPM with Mach-boundary/transport/Smagorinsky
restart [checks](validation/ppm-regular-spray-20260925.md). Serial two-level AMR
hydro and reaction splitting now support PPM, including periodic/physical
boundaries, reflux and exact complete-refinement parity
([checks](validation/ppm-two-level-hydro-20260925.md)). Serial fixed/selected
two-level application input, conservative regrid, transport and checkpoint
selection now support PPM, with fixed reacting restart and exact selected output
[checks](validation/ppm-two-level-application-20260925.md). Sparse two-level MPI coarse PPM and four-layer halos now match serial even with
one owned cell per rank ([checks](validation/ppm-mpi-two-level-coarse-20260925.md)).
Fine-level MPI PPM now supports subcycling, reflux, reaction/transport splitting
and exact serial hierarchy parity
([checks](validation/ppm-mpi-two-level-hierarchy-20260925.md)). Fixed/selected MPI
two-level input and regridding also enable PPM; fixed four-to-two-rank restart
matches serial output ([checks](validation/ppm-mpi-two-level-application-20260925.md)).
A primitive-state Riemann entry now preserves traced pressure and uses the pinned
PeleC contact/regularization policy
([checks](validation/ppm-primitive-riemann-20260925.md)). An optional CTU API
now applies separate single/double transverse pressure work using Godunov face
pressure and velocity ([checks](validation/ppm-transverse-pressure-20260925.md)).
The pressure-traced path now uses the minimum shock flattening across all three
axes ([checks](validation/ppm-multidirectional-flattening-20260925.md)). Pinned
normal tracing recomputes internal energy from pressure through the EOS, as does
PeleF; independent normal internal-energy tracing is not an additional missing
reference step. The existing `characteristic_ppm`/`pelec` application combination
now selects this pressure policy through the common periodic/halo CTU entries
([checks](validation/ppm-pressure-dispatch-20260925.md)). Scalar reconstruction
and profile integration now agree with compiled pinned C++ bodies in a bounded
675-case comparison ([checks](validation/ppm-pinned-kernels-20260925.md)).
The primitive Riemann entry also agrees with the pinned compiled kernel using a
NASA7 EOS adapter in 56 cases ([checks](validation/riemann-pinned-kernel-20260925.md)).
Double-transverse pressure subtraction now follows the pinned spatial order
([checks](validation/ppm-transverse-order-20260925.md)). Regular-grid inflow
faces now receive final fourth-order-slope PLM flux replacement
([checks](validation/ppm-inflow-plm-20260925.md)). Serial/MPI two-level AMR now
forwards physical inflow masks on coarse/fine ownership boundaries and normalizes
PLM face composition before EOS evaluation
([checks](validation/ppm-two-level-inflow-20260925.md)). Serial/MPI patch-tree
common and recursive stepping now also forwards physical inflow masks
([checks](validation/ppm-tree-inflow-20260925.md)). Regular MPI spray hydro now
forwards domain-end inflow masks too
([checks](validation/ppm-regular-mpi-inflow-20260925.md)). External
normal primitive projection now preserves independent species weights and agrees
in a bounded 315-case compiled-kernel comparison
([checks](validation/ppm-mixture-projection-20260926.md)). Regular zero-source
single/double transverse kernels also agree in 144 compiled-reference cases
([checks](validation/ppm-pinned-transverse-20260926.md)). Godunov contact-speed
species upwinding and zero-viscosity final normalization now connect through
regular/EB CTU ([checks](validation/ppm-passive-flux-20260926.md)). The local
double-transverse kernel now accepts primitive sources and agrees in 216 compiled
reference cases ([checks](validation/ppm-transverse-source-20260926.md)); full
periodic/halo CTU source routing now also agrees in 28 full-assembly cases
([checks](validation/ppm-ctu-source-20260926.md)). Cell-local conserved-source conversion
is implemented ([checks](validation/ppm-source-conversion-20260926.md)). Grid
conversion now feeds periodic/halo CTU directly
([checks](validation/ppm-conserved-source-20260926.md)). Periodic CTU stepping
now adds supplied frozen conserved rates with inventory and rollback checks
([evidence](validation/ppm-source-step-20260926.md)). Boundary stepping also
accepts caller-supplied source halos and separates boundary/source inventories
([checks](validation/ppm-boundary-source-step-20260926.md)). The full library
split now forwards these sources with exact transport-composition checks
([evidence](validation/ppm-source-split-20260926.md)). A regular-cell external
force/work constructor now follows the pinned force-per-volume policy
([checks](validation/external-force-source-20260926.md)). Old/new external work
correction and a separate volume-source ledger now connect through full splitting
([checks](validation/external-force-step-20260926.md)). Fixed/selected regular
gas input and source-accounted diagnostics now connect external forcing
([checks](validation/external-force-app-20260926.md)); broader application routing
and full SDC integration remain open. The reference-ordered SDC candidate
assembler is now implemented and checked for source/hydro/reaction component
mapping ([checks](validation/sdc-state-assembly-20260926.md)). Nonreacting periodic
iteration now connects source callbacks, CTU and repeated state assembly
([checks](validation/sdc-periodic-iteration-20260926.md)). Fixed/selected regular
gas input selects this nonreacting periodic path with external force and source
ledgers ([checks](validation/sdc-app-20260926.md)). Molecular transport now
participates in old/new source evaluation with application inventory checks
([evidence](validation/sdc-transport-20260926.md)); forced chemistry and broader
boundary/AMR/MPI support remain open. The forced chemistry RHS now supports
variable density and externally changing internal energy
([checks](validation/forced-reactor-rhs-20260926.md)). An adaptive implicit
substep now integrates the forced reactor
([checks](validation/forced-reactor-integrator-20260926.md)). Full-interval
subcycling now preserves caller state on failure
([checks](validation/forced-reactor-interval-20260926.md)). Periodic SDC library
reaction coupling is implemented with a uniform reacting check
([evidence](validation/sdc-reactions-20260926.md)). Fixed/selected application
routing now retains reaction sources across steps
([checks](validation/sdc-reacting-app-20260926.md)); nonuniform reacting
qualification and broader boundary/AMR/MPI routing remain open. The library
now handles zero-gradient outflow with a separate boundary inventory
([checks](validation/sdc-outflow-20260926.md)). Fixed/selected dispatch now
supports reacting outflow with separate boundary/source accounting
([checks](validation/sdc-outflow-app-20260926.md)). Outflow molecular transport
now participates in old/new source evaluation and boundary accounting
([checks](validation/sdc-outflow-transport-20260926.md)); other physical
boundary types remain open.
SDC prescribed outlet pressure now reaches hydro and transport
([checks](validation/sdc-pressure-outlet-20260926.md)). Specified-state and Mach
inflow now share existing ghost and transport treatment in SDC
([checks](validation/sdc-inflow-20260926.md)). Slip, no-slip and isothermal
wall dispatch now reaches SDC hydro and transport
([checks](validation/sdc-walls-20260926.md)). Full-face hydro/transport fluxes,
coarse/fine halo interpolation and separate reaction/transport stages now feed
a common-step serial AMR SDC loop with final reflux/restriction
([checks](validation/sdc-tree-loop-20260926.md)). Regular and EB tree application routing, persistent reaction-source restart
and serial/MPI recursive gas SDC are now connected; see the
[current scope](current_status.md) and [MPI EB recursive checks](validation/sdc-mpi-eb-subcycle-20260927.md).
Remaining work includes tree small-cell redistribution, spray SDC integration
and broader qualification. Physical
wall/outflow Riemann-state correction now reaches regular serial/MPI AMR and
MPI spray dispatch ([checks](validation/ppm-physical-riemann-20260926.md)).
The full periodic zero-source twelve-solve assembly now agrees in six bounded
compiled-reference cases ([checks](validation/ppm-pinned-ctu-20260926.md)).
The expanded ten-case oracle also covers active flattening and artificial
viscosity ([checks](validation/ppm-pinned-shock-20260926.md)). Eighteen cases now
include inflow/outflow/slip/no-slip boundary assembly; inflow PLM replacement
preserves independent species weights
([checks](validation/ppm-pinned-boundary-20260926.md)). Broader qualification
remains open.
The reference also adjusts final fluxes with compressive artificial viscosity
(default `difmag=0.1`). A regular-cell nodal-divergence/flux-correction API now
implements this step with no-slip tangential exclusions and species-flux
normalization ([checks](validation/artificial-viscosity-20260925.md)). Periodic and
halo CTU now apply coefficient 0.1 by default with the `pelec` pressure policy;
regular, AMR, tree and MPI spray callers provide physical no-slip masks
([checks](validation/ppm-viscosity-dispatch-20260925.md)). The CTU API accepts an
explicit coefficient. Fixed/selected regular gas namelist input now exposes it
([checks](validation/difmag-regular-input-20260925.md)); those executables have
no restart interface. Serial fixed/selected regular spray now forwards the
coefficient through coupled stepping and binds its resolved value to restart
([checks](validation/difmag-spray-input-20260925.md)). Regular MPI spray now forwards
collectively agreed coefficients and binds restart context
([checks](validation/difmag-mpi-spray-20260925.md)). Serial/MPI patch-tree common,
level and recursive hydro now forward explicit coefficients with collective
agreement and rollback ([checks](validation/difmag-amr-tree-hydro-20260925.md)).
Serial/MPI reaction/transport composition and the AMR spray gas call now forward
the coefficient ([checks](validation/difmag-amr-split-20260925.md)). Serial/MPI
tree gas and MPI AMR spray now expose coefficient input and bind
resolved restart identity ([checks](validation/difmag-amr-input-20260925.md)).
The older two-level AMR hydro, reaction and transport composition now also
forward explicit coefficients through coarse/fine stepping and reflux
([checks](validation/difmag-two-level-20260925.md)). Fixed/selected serial/MPI
two-level applications now expose input and bind the resolved coefficient to
restart identity ([checks](validation/difmag-two-level-input-20260925.md)).
Regular serial PCM/PLM library stepping now applies explicit positive coefficients
at both SSPRK2 stages, including mean face fluxes and boundary accounting
([checks](validation/pcm-plm-viscosity-20260925.md)). Explicit zero retains the unadjusted PCM/PLM calculation. Serial/MPI tree common and level hydro now apply
the coefficient before flux integration and reflux, including recursive
subcycling ([checks](validation/pcm-plm-tree-viscosity-20260925.md)). PCM/PLM
application defaults now resolve to the reference 0.1 through shared regular
gas configuration, including AMR checkpoint identity
([checks](validation/regular-gas-default-difmag-20260925.md)). Serial/MPI regular spray
now also defaults to 0.1 with changed-rank restart evidence
([checks](validation/spray-default-difmag-20260925.md)). EB application inputs also now default to 0.1 with zero-to-default restart
rejection ([checks](validation/eb-default-difmag-20260925.md)). Regular boundary/full-composition and PPM
flux APIs now default to 0.1 independently of solver/reconstruction
([checks](validation/regular-library-default-difmag-20260925.md)). Periodic
PCM/PLM direct stepping and the older Strang wrapper now accept the coefficient
and default to 0.1, preserving explicit-zero dispatch
([checks](validation/periodic-api-difmag-20260925.md)). Serial/MPI AMR tree
common/level hydro now resolves omission to 0.1, with common/recursive
serial-versus-MPI equivalence checks
([evidence](validation/amr-tree-default-difmag-20260925.md)). Serial/MPI
older two-level hydro now also resolves and forwards 0.1 on every coarse/fine
stage, with eight-rank default-versus-explicit checks
([evidence](validation/two-level-default-difmag-20260925.md)). Regular MPI
spray hydro/coupling library defaults now also resolve to 0.1, with four-rank
PLM/PPM omitted-versus-explicit checks
([evidence](validation/mpi-spray-library-default-20260925.md)). EB hydro and
hierarchy/coupling libraries now also default to 0.1; transport-only shared
operators retain zero
([evidence](validation/eb-library-default-20260925.md)). The viscosity kernel now accepts
EB volume fractions and open-face masks, suppressing nodal divergence near
cut/covered cells using the reference cutoff
([checks](validation/eb-viscosity-kernel-20260925.md)). The EB face-flux API and
plain/flux-redistributed/state-redistributed Euler APIs now accept explicit
coefficients. Outflow halos are built locally; external positive-coefficient
halos must include corner states and matching volume fractions
([checks](validation/eb-viscosity-dispatch-20260925.md)). Serial common-step EB
tree hydro now requests complete donor halos and forwards explicit coefficients
([checks](validation/eb-tree-viscosity-20260925.md)). Extracted patches now retain same-resolution exterior volume fractions through
serial/MPI halo fill and checkpoint restart
([checks](validation/eb-fine-geometry-halo-20260925.md)). Geometry without that
context still uses donor fractions; external reference geometry fidelity remains open,
as does the uniform default policy. MPI common-step hydro now
fetches complete owner-local halos and applies agreed coefficients, with
three-rank serial parity
([checks](validation/mpi-eb-tree-viscosity-20260925.md)). Serial recursive hydro and full
reaction/transport composition now forward the coefficient and request
time-interpolated complete halos
([checks](validation/eb-tree-viscosity-split-20260925.md)). MPI recursive and
coupled hydro now do likewise with collective coefficient agreement and
three-rank reacting/transport parity
([checks](validation/mpi-eb-viscosity-split-20260925.md)). Single-domain EB gas
input now accepts finite nonnegative coefficients and binds positive values
to checkpoint identity; omitted/zero retains legacy identity
([checks](validation/eb-difmag-input-20260925.md)). EB tree gas input now forwards
the coefficient through serial/MPI composition and patch checkpoint identity,
with a focused serial recursive restart check
([checks](validation/eb-tree-difmag-input-20260925.md)). Single-domain EB spray
now forwards explicit coefficients through particle/gas composition and binds
them to restart context while retaining the legacy zero component sequence
([checks](validation/eb-spray-difmag-20260925.md)). MPI EB spray now forwards
agreed coefficients before gas stepping, with 2/4-rank application parity
and changed-rank restart checks
([checks](validation/mpi-eb-spray-difmag-20260925.md)). Older two-level EB
common/recursive hydro now accepts explicit coefficients with complete patch
halos and conservative reflux
([checks](validation/eb-two-level-viscosity-20260925.md)). Its reaction/transport
split wrappers now forward the coefficient through common/recursive hydro;
see [coupled checks](validation/eb-two-level-viscosity-full-20260925.md).
Multipatch hydro/full composition now forwards coefficients through complete
peer halos, flux reconciliation and reflux
([checks](validation/eb-multipatch-viscosity-20260925.md)). The older
two-level/multipatch application now accepts explicit coefficients and binds
patch restart identity
([checks](validation/eb-two-level-difmag-input-20260925.md)). Regular MPI spray PLM hydro now
applies explicit coefficients with short-slab serial parity
([checks](validation/plm-mpi-spray-viscosity-20260925.md)); serial spray and MPI
mixture wrappers accept those coefficients. Serial two-level AMR now applies
explicit PCM/PLM coefficients on coarse and fine stages before flux averaging
and reflux ([checks](validation/pcm-plm-two-level-viscosity-20260925.md)). MPI
two-level AMR now applies the same corrections using local coarse/fine halos,
with eight-rank serial parity
([checks](validation/mpi-two-level-viscosity-20260925.md)). Shared gas configuration
now accepts explicit PCM/PLM coefficients, and two-level/tree checkpoints bind
them while preserving legacy zero identity
([checks](validation/nonppm-difmag-input-20260925.md)). Regular serial/MPI spray
now accepts explicit PLM coefficients and binds restart identity, with fixed
coupled and changed-rank checks
([evidence](validation/spray-plm-difmag-input-20260925.md)). Regular MPI spray
now exposes PPM through particle coupling and checkpoint context, with fixed
changed-rank continuation [checks](validation/ppm-regular-mpi-application-20260925.md). The regular MPI spray hydro library now accepts PPM
with four-layer exchange across two-cell slabs and serial state parity
[checks](validation/ppm-regular-mpi-hydro-20260925.md); the existing
PCM/PLM solver remains available.

| Order | Work unit | Measurable exit gate |
| --- | --- | --- |
| 1 | Validate existing physics before expanding it | Pin external cases, states, observables and tolerances; compare H2/O2 trajectories/rates with Cantera and flow fields with frozen PeleC; characterize convergence and conservation |
| 2 | Production chemistry path | Specify the target fuel mechanism and unsupported forms; validate ingestion/0D ignition; qualify stiff CFD integration with rollback and conservation before detailed-fuel claims |
| 3 | General 3D flow and mesh | Physical boundaries, then dynamic/multipatch/multilevel AMR and general EB in separate increments; require conservative regrid/reflux, changed-rank restart and convergence |
| 4 | LES and Lagrangian spray | Validate closures and particle submodels, then two-way mass/momentum/energy coupling; qualify evaporation and breakup before penetration/ignition comparisons |
| 5 | Practical release qualification | Reproducible clean install, documented supported inputs, distributed I/O/scaling measurements, physical benchmarks and owner-resolved licensing |

These are dependencies, not permission to bundle several physical models into
one PR. Existing regular 3D and static two-level AMR foundations must not be
reimplemented. Their restrictions remain explicit in current status.

Regular 3D hydro/chemistry/molecular transport now accepts per-face prescribed reservoir inflow, outflow, slip-wall and adiabatic/isothermal no-slip conditions with prescribed tangential wall velocity
with integrated boundary-flux accounting; see the [case guide](../cases/reactive_boundaries_3d/README.md).
Mach-state inflow is connected for regular and patch-tree gas/spray paths;
broader characteristic-boundary qualification and moving geometry remain open.
EB tree halos and serial common-step hydro now also accept Mach inlets
([focused evidence](validation/eb-mach-hydro-20260925.md)). Serial common/recursive
transport, recursive hydro and full splitting now forward those controls
([split evidence](validation/eb-mach-split-20260925.md)). MPI halos, common and
recursive hydro/transport and full splitting now also forward these controls
([MPI evidence](validation/eb-mach-mpi-20260925.md)). Serial/MPI timestep selection
now includes inlet signal speeds and diffusivity
([bound checks](validation/eb-mach-timestep-20260925.md)). Fixed/selected gas
application input and checkpoint binding are connected; fixed serial/2/3-rank
execution and rank-changing restart have
[focused evidence](validation/eb-mach-application-20260925.md). MPI EB spray now
forwards the same controls, with 2/4-rank conservation/parity and 2-to-4-rank
restart [checks](validation/eb-mach-spray-20260925.md). These bounded tests do
not establish nozzle validation or deposition/splash/film inlet qualification.
AMR-particle boundary coupling and stationary level-set EB applications are
implemented within the documented subsets. A three-axis extrusion API now maps
validated 2D level-set sections into 3D volume/face/wall metrics, including faceted
normal integrals and uniform-pressure balance. Nonplanar section-file geometry now connects
to fixed/selected application inputs and geometry-bound restart. Outward-neighbor StateRedist connects
nonplanar hydro/transport/reaction updates in a focused extruded-section case.
A 3D nodal level-set builder now supplies tetrahedrally reconstructed cell/face
metrics and pressure balance. Nodal-file geometry now connects to fixed/selected
application/restart and coupled reaction/transport/hydro updates. MPI EB spray
regridding and changed-rank restart are connected
([dated evidence](validation/spray-mpi-eb-regrid-application-20260918.md)); broader
small-cell stabilization, convergence and physical qualification remain open.

Serial/MPI 3D temperature-tagged patch movement/resizing is implemented with conservative
remapping, overlap reuse and geometry-aware changed-rank restart, including
physical-boundary tagging and boundary-touching patches. Remapping now runs
within fixed coarse slabs without full-level root aggregation. General
patch-tree evolution and integration are implemented as described below;
measured scaling and broader distributed I/O qualification remain open.
Two-level MPI AMR CSV output and final diagnostics now
operate on local slabs, as do analytic coarse/fine initialization and checkpoint read/write; see the [AMR case guide](../cases/amr_reactive_3d/README.md).

The regular 3D patch-tree geometry/state foundation now supports multiple children,
recursive PCM initialization, bottom-up restriction and composite inventories.
Snapshot peer/ancestor ghost filling and physical boundaries are connected.
Per-donor-level temporal interpolation is available for coarse ghosts.
Tree flux synchronization/reflux now handles nested and adjacent boxes, periodic
seams and physical boundary-ledger replacement. Synchronous SSPRK2 hydro now
connects flux production and temperature recovery, using a common timestep
bounded by all node CFL limits. Limited coarse-donor interpolation is connected
for PLM hydro and transport. Single-depth RK2 stepping now exposes integrated
fluxes and coarse-time endpoints. Interval-scoped reflux can return corrected
fluxes while preserving coarser endpoints. Recursive hydro subcycling and scaled
coarse CFL selection are connected. Recursive transport uses ratio-squared substeps
and scaled diffusion/wall bounds. Coupled/application dispatch and mode-aware
checkpoint continuation are now connected. MPI patch-owner scatter/gather is
implemented, and sparse owner-cell queries provide the payload exchange primitive;
synchronous halo donor selection, PCM/PLM interpolation and physical transforms
are connected. Distributed bottom-up restriction with depth limits is implemented.
MPI shared/periodic face flux synchronization is implemented. MPI application subcycling dispatch and mode-aware checkpoint identity are now connected.
Focused hydro-only scheduled/adaptive recursive execution and changed-rank restart
are checked. Optional deterministic cell-step-weighted whole-patch partitioning is
connected at initialization/regrid/restart. Final CSV now uses owner-local MPI-IO;
initial inventory and restart reference scales also avoid field gathers.
Checkpoint saves now use owner-local MPI-IO.
Checkpoint loading now streams an index with bounded root memory and reads node
ranges on their owners; optional indexed V2 eliminates the serial payload scan.
Direct same-geometry owner migration now transfers U/T without interpolation or
a root gather. A cost-driven rebalance API now accepts owner-reported patch costs
and migrates only for a sufficient predicted peak-load reduction. Scheduled recursive full-split execution now has a focused combined regrid,
rebalance and V2 changed-rank restart check. Temperature-triggered full-split
refinement also has this combined check. Particle/AMR application integration is
connected; large-tree scaling and broader EB/particle qualification remain open.
Owner-local kernel timing
and opt-in application rebalance cadence are now connected. Recursive MPI
reaction/transport/hydro composition is connected. MPI recursive transport now connects
ratio-squared steps, interval reflux and scaled diffusion/wall bounds. MPI recursive
hydro now connects ratio-based steps, interval reflux and scaled coarse CFL.
MPI transport level stepping with temporal donors, physical walls and integrated
fluxes is now connected. MPI hydro level stepping
with coarse-time donors and integrated face fluxes is connected. Sparse coarse-time ghost endpoints
and per-depth interpolation fractions are now connected. MPI temperature-adaptation application dispatch
and adaptation-aware changed-rank continuation are connected.
MPI temperature tagging now generates deterministic trees from owner-local hot
exposed cells; fixed/selected application controls and restart policy are connected.
Scheduled MPI tree remapping and event-aware checkpoint continuation are connected.
Sparse conservative MPI tree remapping now supports changed geometry/ratios/owners
with collective rollback and reactive temperature recovery.
Fixed/selected static MPI tree restart now supports changed rank counts and serial interchange.
Fixed/selected-mechanism static MPI tree application dispatch and final CSV output are connected. Distributed
transport/wall timestep selection and the full common-step operator composition
are connected. Synchronous molecular-transport RK2 is connected with moving/isothermal
physical walls and serial PCM/PLM donor parity. Exposed-cell
MPI chemistry and reaction/hydro Strang composition are connected. Synchronous
distributed SSPRK2 hydro is connected with serial PCM/PLM parity. Owner-local EOS recovery
and global minimum CFL selection are implemented. Parent flux
averaging, conservative cell correction, physical boundary accounting and depth-limited
distributed reflux are connected. Broader qualification remains. Exposed-cell chemistry and reaction/hydro splitting are connected;
synchronous molecular transport now also connects nested reflux and physical
walls. The transactional reaction/transport/hydro split API is connected.
Conservative overlap remapping now handles topology/refinement-ratio changes.
The serial application now connects static topology input, analytic initialization,
split evolution, scheduled conservative regridding and exposed-cell output.
Configure-time selected-mechanism dispatch now uses an isolated tree runtime.
Tree snapshot I/O, application stop/restart and schedule-position checks are implemented. Temperature-driven tiled topology generation, successful-step application cadence
and automatic-geometry restart are connected for serial and MPI ownership. See the
[tree design](design_decisions/amr-patch-tree-3d.md).

Particle Cartesian MPI migration preserves parcel/component payloads and local
ID ordering, with collective input/capacity rejection. Distributed gas feedback,
thermodynamic stepping, scheduled injection, breakup and rank-changing restart
are connected to applications. AMR film ownership, regridding/rebalancing and
restart are also connected
([dated film evidence](validation/spray-amr-film-application-20260918.md)). These
implementation results do not establish experimental spray/film validation,
production scaling or complete reference-solver parity.

## Contract for every increment

Before modifying code, specify the failure/capability, supported dimensions,
physics/parallel mode and a measurable exit test. Include a focused unit test,
an application/conservation or external-parity test, and a precise scope
statement. Separate numerical changes from structural refactoring.

For comparisons, record reference revision, mechanism/source hash, initial
state, observable, tolerance, comparison time and mesh/time-step convergence.
Do not loosen tolerances or regenerate golden outputs merely to obtain green
CI. Failed, interrupted, missing and completed checks are distinct states.
New work starts from the integrated main; a deliberate PR dependency must
state its base and merge order rather than recreating an accidental stack.

## Historical detail

The old milestone-by-milestone roadmap remains at the
[immutable checkpoint](https://github.com/GTT10/Fortranslate/blob/509d07929d52d9e570a6af1f6ec7e5836bd4fe5d/docs/completion_roadmap.md).
Detailed rationale remains in [design decisions](design_decisions/), the
[historical implementation inventory](implementation_status.md) and dated
[validation records](validation/). [Issue #101](https://github.com/GTT10/Fortranslate/issues/101)
tracks completion and actual tooling limitations. The unexecuted ancestor-only
branch-deletion audit is an administrative remainder, not a numerical feature.
