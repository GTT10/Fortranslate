# Current status: PeleF 0.246.0

This page separates implemented capability, dated reported results, and
independent checks. The preceding integration is
recorded in [PR #102](https://github.com/GTT10/Fortranslate/pull/102). The new
opt-in spray scope and qualification boundary are in
[the spray guide](spray.md) and [0.246.0 validation](validation/0.246.0.md).
Older milestone prose remains historical detail, not a current global status.

The validated 2026-11-06 parallel integration changes have been applied to this
checkout after user approval. Results and remaining qualification boundaries are
recorded in [the integration evidence](validation/parallel-integration-20261106.md).


The 2026-11-04 elementary-kinetics library exposes an explicit equilibrium
reference-concentration constant, preserving the default convention. Clean
OpenMPI/CVODE focused tests and n8 serial/MPI3 reference checks pass; eight
CSV pairs and42 frozen reactor outputs (21 per coordinate formulation) match
their respective comparison baselines exactly. The isolated constant change
does not resolve the captured reactor discrepancy. It is not an application
or restart option and does not close C2; see the [dated evidence](validation/tgreact-amr-reference-20261019.md).

The full-coordinate CVODE library also accepts explicit species/temperature
absolute-error weights capped by the caller tolerance. Default application
outputs remain byte-identical in n8 serial/MPI3 checks. This is a library-only
control. A separate first-MOL source-assembly diagnostic now passes the
original n64 field/budget gates through50ns with explicit compatible chemistry;
maximum x-momentum difference is9.6030e-7. The isolated MPI equality claim was invalidated: its probe lacked the MPI
compile definition and dispatched serial code. The fresh standard-source build
passed11 tests and both n8 reference gates, but exposed the missing corresponding
MPI source-assembly change. Both serial and MPI paths are now updated: a fresh build passes11 focused
tests and standard n8 reference gates, with8 exact CSV pairs and verified MPI
dispatch. Reacting MOL fixed/selected common/recursive rank-neutral restart
also passes after correcting the MPI reflux diagnostic scope and checkpoint
test helper. Current-API compatible n8 serial/MPI3 checks pass with8 exact pairs;
the current-API n64/50ns comparison also passes all original gates with6 exact
serial/MPI3 CSV pairs. The16-grid compatible300ns trajectory now passes all13 original snapshot
gates. The16/32/64 fixed-step audit passes all original gates and all256 positive-time
norm comparisons decrease, but maximum observed orders remain low and closely
match native behavior. The16-grid100/50/25ns temporal audit passes all trajectory gates and all192
norm comparisons decrease, with sub-second-order maximum ratios also seen in
native. The12.5ns extension also passes all25 snapshots; temporal maxima
localize near the coarse/fine boundary, without proving an interface-only cause.
The32-grid300ns compatible serial trajectory completed on2026-11-06 with all13 snapshot gates passing; this does not qualify the canonical mechanism.
Longer trajectories and convergence remain open;
this is not C2 completion.

The 2026-10-29 bounded evaporation comparison identifies and fixes two CTU
issues: incomplete tangential transport neighborhoods and retained SDC
inverse-temperature seeds in the nonreacting pinned-EOS path. A16-cubed
trace-species failure additionally exposed gradient-product underflow in the
PLM limiter; direct sign comparison now preserves the nonnegative-face guard.
Focused regressions pass and the formerly rejected32-step case completes
with exact MPI1/MPI2 state, particle and budget agreement. A16-cubed MPI2-to-MPI3
checkpoint continuation also matches the uninterrupted run after aligning
the context assembly capacity with the existing reader; changed geometry
still rejects restart. Fixed-reference
errors and spatial convergence remain open; see the
[evaporation evidence](validation/evaporation-full-reference-20260927.md).

The 2026-10-30 regular-AMR gas tree path now supports explicitly enabled
isothermal slip faces with SDC/MOL, common/recursive timesteps and restart
binding. Fixed/selected application checks each pass68 cases, including
rank-changing continuation and changed-boundary refusal; four subcycle and
checkpoint tests pass. A pre-existing nonrotating SDC checkpoint issue was
resolved by explicitly validating the same covered-cell energy reset used
by MOL. Bounded n4/500 K fixed-PeleC comparisons pass with exact serial/MPI3
outputs. The n8 common/recursive and n16 recursive comparisons also pass,
with24 exact serial/MPI CSV pairs. The500 K wall spatial audit remains
nonconvergent, with energy/temperature maxima at the physical wall. A2050 K
diagnostic passes n32 serial/MPI3 reference gates, but8/16/32 z-momentum
maximum differences still grow at100/200ns near the coarse/fine interface.
Halving the timestep does not remove this spatial nonmonotonicity. Common
and recursive stepping both retain approximately first-order maximum temporal
differences; stage diagnostics identify delayed synchronization pressure
feedback as a candidate cause, without an adopted numerical correction. The
2026-10-31 coupled16/32/64-grid,100/50/25ns diagnostic has decreasing L1/Linf
differences in all224 rows for both solvers, but z-momentum maximum orders
remain low (0.21--0.61) and the n64 trajectory fails17 original reference
field gates. Paired stage observations identify amplification of tiny H2O2
neighbor differences in a coarse/fine ghost; same-input interpolation slopes
match the fixed AMReX kernel. Identical-input reactor replay explains the
leading H2O2 discrepancy. An isolated full-coordinate integrator improves
the first local ghost error but fails n64/50ns original momentum and energy
limits at a different interface location; it is not adopted. A separate
full-coordinate CPU policy correction now matches the pinned 100-step Jacobian
refresh interval instead of the CVODE default of51. Three focused tests and
15 expected serial/MPI lifecycle outcomes pass, including exact changed-rank
restart states. Post-correction n8 serial/MPI3 reference gates pass, but
n64/50ns still fails x-momentum (1.306e-6 versus1e-6); the two prior energy
failures now pass. A separate native-compatibility diagnostic reduces the n64/50ns
x-momentum maximum to1.005e-6, still above1e-6; it is not adopted. Same-input
interpolation slopes match AMReX, while nearly flat trace-species gradients
amplify small evolved differences. A subsequent production primitive-conversion change preserves admitted rhoY/rho
for nonrotating regular MOL, matching the native conversion policy. Two focused
tests, four restart regressions and fixed/selected thermal-slip checks
(68 outcomes each) pass. The reduced-coordinate n64/50ns comparison still
fails x-momentum (1.309e-6 versus1e-6, slightly above its prior1.293e-6);
the full-coordinate comparison also fails x-momentum (1.983e-6) plus density
and both energy gates. Same-input pressure improvement has not translated
into evolved-trajectory improvement. Paired stage/flux capture identifies
shared-limiter amplification and predominantly interface-normal viscous flux
differences. Same-input reactor replay attributes most early temperature
difference to the pinned-native activation convention. Explicit activation-only
and caller-capped native-compatibility diagnostics still fail x-momentum
(1.764e-6 and1.100e-6 versus1e-6); neither qualifies canonical C2.
The remaining four interface-cell exceedances are traced to shared-limiter
amplification of small evolved differences. A native-only diagnostic changing
pressure multiplication association gives50ns momentum differences8.04e-11,
9.49e-10 and1.052e-6 on16/32/64 grids at fixed25ns steps. This demonstrates
roundoff sensitivity near the current gate, without qualifying PeleF or
relaxing acceptance. The canonical fixed-dt16/32/64 spatial analysis now
has decreasing L1/Linf differences in all128 positive-time state/increment
rows for both solvers at25/50ns. Momentum-increment maximum orders remain
as low as0.059, so this bounded trend does not establish asymptotic convergence.
Original reference gates fail at n16/25ns and n64/50ns; n32 passes both times.
The n16/n32 MPI3 comparison reproduces all12 serial state/budget CSVs
exactly and retains the same reference failures.
At the n16 maximum, paired observations trace the first25ns discrepancy mainly
to normal hydro flux and H2O2-controlled ghost limiting. Native-input reactor
replay on actual candidate inputs attributes the leading local H2O2 difference
to the integrator. A canonical-mechanism full-reactor component-tolerance
control passes n8/n16 but fails four n64 gates, with x-momentum2.214e-6.
It is not adopted. A separate native-pressure-order control reduces the
n64 maximum to1.467e-6 but still fails momentum and both energy gates; it
also remains isolated. Small reactor replays show that tighter absolute
tolerances reduce differences between candidate integrations but leave a
canonical/native trace-species discrepancy. Reaction-rate convention differences
must be separated from integrator error; C2 remains open.
The generic restart checkers were separately updated for timestep-history
headers and uninterrupted/restarted step-count equality.
Convergence qualification and C2 remain open. See the
[dated evidence](validation/tgreact-amr-reference-20261019.md).

The 2026-10-29 regular-AMR wall comparison now passes original n8/n16/n32
trajectory criteria after correcting transport primitive recovery and no-slip
tangential ghost reflection. Serial/MPI outputs agree exactly in the bounded
wall cases. Wall-transport MOL restart also passes, including MPI2-to1, with
explicit validation of covered-cell internal-energy reset. Spatial convergence
remains open: both solvers retain nonmonotone z-momentum increment maxima. Stage/flux and
native-donor diagnostics identify shared limiting of ghost normal velocity. An
isolated momentum-limiter variant has decreasing three-grid differences through300ns but
fails original reference limits and is not adopted.
The extrapolation/outflow n4/n8/n16 comparison also passes original trajectory
and candidate open-ledger gates through300ns, while density/energy/temperature
increment maxima remain nonmonotone at physical boundary cells.
Halving the time step preserves reference agreement but does not remove that
spatial nonmonotonicity. The n16/dt50ns MPI3 comparison passes original gates
and exactly matches all14 serial state/budget CSVs. Serial n32/dt50ns also
passes reference gates; 8/16/32 maximum norms decrease at250/300ns but retain
early-time nonmonotonicity. The n32 MPI3 comparison also passes original
gates and exactly matches all14 serial state/budget CSVs.
See the [dated evidence](validation/tgreact-amr-reference-20261019.md); these
results do not qualify EB, physical validation or a full release matrix.

The 2026-10-28 ordinary regular-AMR MOL correction preserves synchronized old
states through nonmutating thermodynamic validation and uses total-energy
pressure for nonrotating faces. Current n4/n8/n16 reacting-transport reference
comparisons pass through 300 ns; focused serial/MPI tests and rank-changing
restart pass. Direct MPI2/MPI3 reference comparisons also pass, including a clean OpenMPI build, common clocks, n16 and n8/dt25 ns. Stationary adiabatic slip/no-slip walls pass field and boundary-ledger checks, and 11 focused particle/LES coupling regressions pass. These bounded checks do not constitute physical validation. The n8 100/50/25-ns study shows decreasing L1/Linf temporal
self-differences for all 16 observables in both solvers. Spatial state norms
also decrease over n4/n8/n16, but x-momentum evolution-increment Linf is
nonmonotone at 100/200 ns in both solvers. This does not close C2 spatial
convergence or the separate EB gates. See the
[regular-AMR evidence](validation/tgreact-amr-reference-20261019.md).

The 2026-10-27 explicit native LiDryer activation-basis compatibility run
now passes all 17 fixed-PeleC field times through 200 ns on the 32-cell mesh
with 12.5-ns steps and unchanged canonical solver tolerances, plus the original
inventory gates. Prefix/restart and checker-input hashes were independently
verified. The frozen MPI2 executable independently reproduces all 17 serial CSVs and passes the same field, clock and inventory gates through 200 ns; this does not certify a newly rebuilt release matrix. This uses the pinned native calorie gas constant for activation
conversion; the canonical SI mechanism remains unchanged. The corresponding
25-ns trajectory passes through150ns but fails total/internal energy at175/200ns,
while inventory passes. Focused nonuniform serial/MPI2-to-MPI3 restart
and trace compatibility checks passed. The6.25-ns-step comparison through100ns
fails field gates from50ns although inventory passes. Three-spacing L1
self-differences decrease for all16 gas/species observables; some Linf
species differences are nonmonotone in both candidate and pinned native.
The native repeat matches all17 exports; a one-ULP initial-temperature
perturbation remains within the original gates. Matching the candidate pressure
arithmetic to native CGS order increases the 50-ns momentum/energy errors by
about22% and is not adopted. Same-input inverse-EOS replay of1848 captured
states differs by at most1.82e-12 K. A pressure-order-only native control also
fails from50ns, with matching worst-cell locations and strongly aligned error
fields. This demonstrates shared arithmetic sensitivity, but does not excuse
candidate failures or close the unchanged fixed-reference gates. First-step
corner-halo tracing identifies a trace-H2O-controlled shared limiter: even
same-input reactor differences below9e-22 kg m^-3 affect slopes of order
1e-20–1e-19. A reusable stencil-conditioning audit now quantifies component
substitutions; it changes neither interpolation nor acceptance thresholds. C2 mesh/time
convergence and C3–C5 remain open; this is reference agreement, not physical
validation. See the [dated evidence](validation/eb-full-cvode-application-20261024.md).

The 2026-10-23 ordinary nonrotating EB-tree MOL correction now uses
consistent total-energy-derived pressure for face and wall fluxes. The rebuilt
serial probe passes the original pinned TGReact 300 ns field and budget gates;
seven focused MPI/restart checks and kernel/tree tests pass. Rotating and
experimental halo-policy paths retain their held-energy convention. The 32-cell, 100 ns original-reference pilot fails six field gates, with
maxima at full-volume coarse/fine boundary cells. Mesh/time
convergence and broader C2 qualification remain open; see the
[dated pressure/flux evidence](validation/eb-default-localization-20261023.md).

The 2026-10-24 ordinary selected-mechanism EB application now passes the
pinned 32-cell/100-ns initial/final field and inventory gates for 32256
composite fluid cells, with imported Simple transport and the documented
chemistry settings. This is one-step reference agreement, not long-trajectory
or mesh-convergence qualification. Nonuniform 8-cell MPI2-to-MPI3 checkpoint
continuation now matches uninterrupted serial CSVs through 200 ns; this is
restart equivalence, not a 200-ns fixed-reference comparison. A separate
ordinary 8-cell trajectory passes all fixed-reference outputs through 300 ns
with 100/50/25-ns steps and final inventory gates. Temporal self-comparison
shows roughly first-order hydro/thermal behavior shared with pinned PeleC;
H2/N2 maximum differences are nonmonotone at the coarse/fine interface in
both solvers, so all-field time/mesh convergence remains open. See the
[ordinary application record](validation/eb-full-cvode-application-20261024.md).

The 2026-10-26 ordinary nonrotating MOL transport path now recovers normal-face
temperature/pressure from total energy in serial and MPI. Saved-field replays
identify a held-energy mismatch after restriction upstream of trace-species
limiter amplification. The rebuilt serial n32/dt50ns/100ns application passes original field and
inventory gates, with all five ordinary outputs byte-identical to the isolated
pilot. All n8/n16/n32 cases at dt100ns/100ns pass original field and inventory gates.
Fixed-solver maximum differences remain nonmonotone; conservative inter-grid
self-comparisons reproduce the pinned solver trend, including a growing
z-momentum maximum despite declining volume-weighted norms. The n32 dt25ns
run also passes all five field times through100ns and inventory gates. Across
100/50/25ns steps, all15 fields have decreasing temporal L1/L2 differences,
but z-momentum Linf remains nonmonotone in both solvers. Removing initial-profile
sampling from the mesh comparison leaves a y-momentum evolution Linf ratio
of0.9483, so maximum-norm spatial convergence remains weak. The n32 dt25ns
checkpoint extension passes through150ns but fails energy gates at175ns and
density/z-momentum/energy gates at200ns; cumulative inventory and restart
identity pass. Stage diagnostics reproduce the upper-boundary momentum error growth from
actual forcing. At fine(26,14,32), the dominant137.5ns-onward discrepancy is
the non-hydro transport contribution, opposing the hydro contribution; a
local replay now reproduces both solvers' transport fluxes and identifies
upper-halo tangential velocity/gradient differences. Timed donor replay identifies a z-momentum-controlled shared x-limiter
(candidate0.491725, native0.495741). The donor discrepancy grows across125ns
coarse/fine synchronization. Both synchronization increments close from their
recorded face transfers; the discrepancy is dominated by incoming fine
transport flux differences. Activation-basis-only sensitivity does not resolve the discrepancy: at150ns
y-momentum Linf grows by1.284, although z Linf falls. A native density-atol
comparison changes y Linf by only0.14%. First-predictor donor replay
identifies a separate H2-controlled y-limiter; same-input reactor replay
matches native H2 closely when both chemistry conventions are aligned. A
combined first25ns decision test reduces y-momentum Linf by32%, but z Linf by only0.65%. Three-cell source capture reproduces both ordinary outputs; coarse(1,5,17) H2 forcing differences are dominated by hydro, initially in the z direction. Native coarse-face capture reproduces both outputs and locates a roughly7e-7Pa reconstructed-pressure difference at the first predictor; pressure-only face replay removes the first-predictor H2 z-flux discrepancy. Converged EOS replay attributes the pressure difference to species-state differences. The combined-chemistry trace removes the coarse first-predictor discrepancy but retains the corrected-stage neighboring-cell pressure difference and fine-halo tangential velocity error. Direct native reactor-input replay reproduces the matched chemistry endpoints; remaining discrepancies originate in forcing. At coarse(1,5,16), pressure substitution removes the first-predictor x-face mismatch; the adjacent(2,5,16) has a1.82e-5Pa pressure difference despite near-identical density/composition. Conservative-state capture confirms roundoff agreement; common old/native temperature guesses recover matching pressure. Actual-temperature replay reproduces the pressure error. Valid-state admission was found to rewrite the stored inverse-EOS seed. Both preserve-seed pilots pass original25ns gates. Canonical z Linf falls83.5%, while y Linf rises0.8%. The scoped serial/MPI admission fix is implemented; the new serial/MPI admission-and-rollback regression and standard held-reset test pass. Rebuilt application diagnostic checks pass; fresh rank-changing restart matches serial through200ns. The new canonical200ns trajectory still fails175/200ns field gates despite passing inventory; z Linf at200ns is10% higher, so late amplification remains unresolved. The matched-chemistry sensitivity also fails175/200ns gates; at200ns z Linf is4.987e-6,38.2% above the canonical preserve-seed result. Its100..200ns fine-cell momentum growth closes against native face-divergence differences within2e-15; transport dominates the z-error growth. A canonical donor capture reproduces seven ordinary outputs and shows that coarse mz differences change the shared limiter; replacing only mz donors removes the transverse halo velocity discrepancy. The largest donor change occurs across125ns synchronization; captured incoming fluxes reproduce that jump within2.8e-16, with two fine-face transport contributions dominant. Local transport replay reproduces both fluxes and removes99.2% of the dominant face discrepancy by replacing transverse velocity gradients. Saved100ns donors exactly reproduce the offending halo velocities; the upstream coarse mz discrepancy already appears at25ns and is dominated by synchronized hydro fluxes. Halo pressure substitution removes their main discrepancy. Captured-state EOS replay attributes most of that pressure difference to the temperature seed and stopping path; the25ns halo-seed pilot passes original gates and reduces the target coarse mz error91%, while global maxima barely change. The scoped serial/MPI implementation passes3 halo tests,8 application diagnostics and n8 rank-changing restart; related MOL/transport tests pass13/14 with one pre-existing retained-energy failure. Production and pilot outputs plus all four diagnostic streams match exactly through25ns. The isolated200ns pilot passes all gates through175ns and all momentum/density gates at200ns; energy remains outside the unchanged limit (scaled1.198e-8 versus1e-8). The production25ns checkpoint continuation reproduces all eight pilot outputs and its budget CSV exactly through200ns. At the energy-maximum cell,100..200ns energy-error growth closes against actual forcing differences within1.1e-12; the non-hydro transport contribution dominates. Transport-state replay through150ns reproduces both fluxes and locates the remaining energy discrepancy in upper-halo temperature/composition, dominated by conduction and species transport. Timed donor replay traces these differences to the mz-controlled shared limiter: replacing only coarse mz donors reduces the150ns halo temperature discrepancy from1.18e-3K to5.13e-8K. The next check isolates the remaining125ns incoming fine-face flux discrepancy. Chemistry changes remain unadopted. Existing pinned-EOS ledger and two Godunov reactor failures were independently reproduced before the change. Chemistry conventions remain unchanged. First-iteration face-average energy forcing was superseded by direct input capture. No limiter, mechanism
or synchronization correction has been adopted. Focused tests report42/43 passes, with the remaining hydro bitwise
failure reproduced on pre-change source. Eight diagnostic checks and a fresh
nonuniform MPI2-to-MPI3 restart comparison against serial through200ns pass.
These are distinct evidence classes, not C2 closure; see the
[dated record](validation/eb-full-cvode-application-20261024.md).

The full-coordinate CVODE RHS now recovers temperature from affine energy at each evaluation, matching the pinned reactor policy; its temperature coordinate is the EOS seed. An independent off-manifold analytic oracle and full/reduced reactor tests pass2/2; application diagnostics and rank-changing restart pass. However, the numerical-Jacobian200ns trial fails original field gates from25ns onward, despite closed inventory. Restoring the pinned independent-temperature approximate Jacobian policy passes original25ns field/inventory gates, full/reduced unit tests, application diagnostics and rank-changing restart. The verified25ns checkpoint continuation passes all field gates through175ns and density/momentum gates at200ns; energy alone remains outside its unchanged limit (scaled1.2165e-8 versus1e-8). Relative to the pre-RHS halo-seed result, energy error rises1.6% and z-momentum maximum14.5%. The current-RHS native-chemistry sensitivity also fails200ns energy, increasing its error13.2%; it is not adopted. Production chemistry remains unchanged. Current-source8/16/32 comparisons at dt100ns/100ns pass all original field/inventory gates. Conservative adjacent-grid z-momentum Linf grows by1.90965 while L1/L2 shrink; pinned PeleC reproduces the same trend. Initial-sampling subtraction leaves y-Linf ratio0.94832. At dt25ns the z-Linf ratio remains1.72881 in both solvers, with declining L1/L2 and y-evolution Linf ratio0.94897; smaller time steps do not remove the spatial maximum issue. A n32 dt12.5ns200ns temporal comparison is pending. Long-trajectory reference agreement and C2 remain open.

The 2026-10-27 fixed-step final-clock repair permits only machine-scale
roundoff at completion without changing the stored time or numerical state.
Five serial/MPI/restart acceptance/rejection cases passed. The saved 32-cell,
12.5-ns-step trajectory was recovered without another numerical step; its
state and cumulative ledger match the original step-16 snapshots exactly.
All 17 reference field checks and the inventory check are now available:
fields pass through 112.5 ns, z momentum fails from 125 ns, and at 200 ns
density and total/internal energy also fail the original limits. Halving
the timestep therefore does not resolve C2; see the
[dated recovery and comparison](validation/eb-full-cvode-application-20261024.md).

## Supported boundary

The gas EB-tree app exposes an experimental `separate_momentum` coarse-halo
policy for held MOL with redistribution `none`. Fixed/selected common and
recursive lifecycle, serial/MPI retry/replay, changed-rank restart and policy
mismatch refusal checks pass; the default remains false. This changes the
reference numerical policy and does not establish original-PeleC parity or
physical convergence. See the [application evidence](validation/eb-momentum-application-20261021.md)
and [design](design_decisions/eb-halo-momentum-limiter.md).
The additional `consistent_velocity` option now connects density-coupled
momentum interpolation to ordinary serial/MPI held MOL, including common and
recursive stages, saved-policy identity and checkpoint admission. It requires
`separate_momentum`; both default to false. Fixed/selected common and recursive lifecycle,
selected reacting/regrid and changed-rank restart checks pass; see the
[velocity application record](validation/eb-velocity-application-20261022.md).
`paired_pressure` now connects fluid-volume-restricted EOS pressure to both
covered coarse spatial states and fine halos in ordinary serial/MPI
common and recursive held MOL. Recursive stages use synchronized covered
pressure preparation and endpoint EOS pressure interpolation for ancestor
donors. It requires both momentum options and a nonrotating frame. Input valid leaves remain unchanged
by reconstruction, and saved stages/checkpoints bind the policy. See the
[pressure application record](validation/eb-pressure-application-20261022.md)
and [recursive/native MOL evidence](validation/eb-pressure-recursive-20261022.md).
The 200 ns nonisothermal native control now has completed 8/16/32 spatial
refinement evidence: paired-pressure error falls from 0.001271 to 0.0001390
to 0.00001456 Pa, with correction-accounted normalized budgets below 3.11e-13.
This is bounded pressure convergence, not formal third-order accuracy or
reacting fixed-reference qualification; see the
[mesh evidence](validation/eb-pressure-mesh-20261022.md).
The completed 300 ns reacting fixed-reference comparison fails the original
field gates: paired temperature difference is 0.239 K versus 0.000110 K
with the default policy, despite closed energy budgets and close H2 consumption.
The four-policy ablation places most of this difference at the prerequisite
separate-momentum limiter (0.241 K); pressure reconstruction slightly reduces
it. The option remains experimental; see the
[reacting comparison](validation/eb-pressure-reacting-20261022.md).
A recorded-stencil audit independently reproduces both limiter outputs to
1.12e-16 and finds 48/384 scalar factors changed by momentum separation;
this is a broad interpolation-policy change in resolved variable-velocity
flow, not just a roundoff fix. See the
[attribution](validation/eb-momentum-attribution-20261022.md).
The independent analytic-shear 8/16/32 comparison is complete for all four
policies: at 32, paired x-momentum and pressure L2 errors are 31.4% and
27.9% of default, but z-momentum L2 error is 151.1% and transverse Linf
errors increase under refinement. All 12 correction-accounted budgets pass
(maximum 2.03e-13). This supports bounded integrated-error improvement,
not all-field convergence; see the
[analytic shear results](validation/eb-shear-analytic-20261022.md).
The completed 16-grid localization/time-step/regular comparison places the
transverse maxima beside coarse/fine interfaces. Halving the cut-case step
changes paired maxima by less than 0.3%; regular geometry reduces paired
y/z-momentum and pressure maxima by about 47%/79%/79%, but leaves interface
errors and nearly unchanged default maxima. The initial-only operator has
zero transverse rates in all eight cases, while streamwise/density/energy
rates are already nonzero. The first-step trace locates resolved pressure
error in the first fine predictor before closure/reflux, followed by transverse
rates in its next spatial evaluation. The sampled predictor balance closes
to roundoff: pressure drops before recover, predominantly through the density
update, while the subsequent E=H+KE reset leaves rho/H/pressure unchanged.
An isolated connected-direction halo candidate improves every nonzero cut
field norm with momentum separation at 50 ns and preserves the regular
control, but worsens default y-momentum L2 by 5.47%. The default remains unchanged.
The 8/16/32 study now reaches 200 ns: separated-policy L1/L2 errors decrease
on both refinements, but y-momentum Linf is nonmonotone and z-momentum Linf
orders remain only 0.02–0.23. The explicit `connected_directions` option now
reaches serial/MPI common and recursive stepping and checkpoint/history identity;
all nine focused fixed/selected application and unit checks pass, including
changed-rank restart and policy mismatch refusal. See the
[application record](validation/eb-connected-application-20261022.md).
The production serial recursive path now reproduces all 32 recorded diagnostic
CSV tables at 16 cells/50 ns for cut/regular shear: enabled policies match the
frozen candidate and default matches the original. Budgets stay below 6.06e-14;
see the [production-path comparison](validation/eb-connected-shear-integration-20261023.md).
A saved-data audit locates all 144 transverse/pressure maxima beside
coarse/fine interfaces. Paired y maxima lie in cut cells, while z maxima lie
in full-volume cells at the x/z coarse/fine intersection; shrinking weighted
error concentration explains why L2 improvement does not establish Linf
convergence. See the [local-error audit](validation/eb-shear-local-convergence-20261023.md).
A 16-grid/50-ns stage and closure audit now reconstructs both fine-probe
momentum trajectories from face rates and applied closure increments to 2.51e-16.
Their initial transverse changes arise before closure; saved initial halos
identify a 0.190915 m/s velocity mismatch at the limited smooth shear extremum.
The numerical policy is unchanged; see the
[operator attribution](validation/eb-sync-increment-20261023.md).
An isolated screened-quartic halo candidate reduces transverse Linf errors
at 16 cells/50 ns to about half the connected baseline, but increases
z-momentum L1 by 38.4% and changes admission after the first stage. It remains
outside production; see the [candidate evidence](validation/eb-quartic-evolution-20261023.md).
A donor audit attributes that switch to near-zero transverse momentum setting
its own roundoff scale. Using the momentum-vector magnitude in the isolated
screen retains admission and reduces y/z Linf to 6.74%/14.0% of the connected
baseline at 16 cells/50 ns; z L1 falls to 25.7%. This remains a single-grid
candidate result; see the [roundoff-scale evidence](validation/eb-quartic-roundoff-20261023.md).
A recovered 32-cell/50-ns comparison now covers all four policies. Candidate
y/z Linf decreases by about two from 16 to 32 cells, while pressure Linf
decreases by about four; maxima remain coarse/fine adjacent. No candidate
norm worsens against the connected baseline. This remains isolated serial
evidence, not C2 closure; see the [mesh comparison](validation/eb-quartic-mesh-20261023.md).
General geometry and local-error convergence qualification remain open. Geometry and wall position change
together, so this does not isolate a single operator or close convergence.

The selected H2/O2/N2 EB-AMR application now has common/recursive reacting
regrid and serial/MPI changed-rank restart evidence, with independent species,
element and UEINT checks. A higher-temperature CVODE case also completes
regridding and 2-to-3-rank restart, consuming 0.023% of the initial hydrogen
over 300 ns. This is short-time reacting integration evidence, not ignition
validation; see the [application record](validation/reacting-eb-amr-20261019.md).
Gas EB-tree applications also support opt-in intermediate composite fields,
budgets and physical-time metadata, including absolute-step continuation
after restart. The [snapshot record](validation/eb-tree-snapshots-20261019.md)
checks intermediate serial/MPI agreement and unchanged final results;
nonuniform EB-AMR trajectory agreement with PeleC remains open.

A matched LiDryer case now compares four physical times against the pinned
PeleC EB-C3 application. The uniform 2000 K reacting planar-EB case agrees
within `7.78e-8 K` in temperature and `4.82e-10` in species mass fraction;
hydrogen consumption agrees within 0.0023%. This closes only the uniform,
single-level reacting comparison, not spatial-gradient or AMR parity. See
the [reference record](validation/eb-lidryer-reference-20261019.md).

A spatial-gradient, two-level regular AMR hydro control now matches pinned
PeleC through three recursive one-iteration MOL steps (`1.14e-15` maximum
scaled field difference). Two-iteration nonreacting MOL fails the unchanged
reference gate: PeleF conserves the accepted two-half transfer while the
pinned branch registers only the old half. This remains an open numerical
policy difference, not qualified parity; see the
[AMR comparison record](validation/amr-mol-reference-20261019.md).

The nonuniform LiDryer TGReact control also passes three recursive two-level
MOL2 steps against pinned PeleC, with `1.21e-5 K` maximum temperature difference
and 0.00225% hydrogen-consumption difference. Initially this imported the
complete reaction history ([operator record](validation/tgreact-amr-reference-20261019.md)).
Ordinary serial/MPI tree applications now compute fresh MOL reaction histories
independently with level-specific intervals, transactional initialization and
no reseeding on restart. The computed-history reference comparison and selected
LiDryer two-to-three-rank/serial restart pass; see the
[initialization record](validation/mol-initial-reactions-20261019.md).
Gas EB tree startup now also initializes fluid-cell reaction history with
level clocks, preserves public gas/T/UEINT and restored history, and rolls
back failed initialization. Recursive hot CVODE regrid/rebalance and
2-to-3-rank continuation pass; see the
[EB initialization record](validation/eb-initial-reactions-20261019.md).
A nonuniform planar EB-AMR comparison with pinned TGReact now runs reaction,
hydro and matched Simple molecular transport through three recursive MOL2
steps. Initial fields, volume fractions and mass/element budgets match their
gates, but density, momentum, energy and temperature trajectories fail the
unchanged field limits. A transport-off control retains the main discrepancy,
largest initially in coarse cut cells; see the
[open comparison record](validation/tgreact-eb-amr-reference-20261019.md).
Ordinary held closure now applies the missing flux re-reflux redistribution
in serial/MPI. The 300 ns density and x-momentum errors decrease by about
72 and 10 times respectively, while the field gate remains open. Remaining
maxima are on the fine level; see the
[correction and qualification record](validation/eb-no-redist-rereflux-20261019.md).
Held MOL now applies limited coarse-to-fine interpolation in connected EB
neighborhoods, in serial/MPI. Initial fine-halo conserved quantities match the
native diagnostic to roundoff, but the 300 ns temperature difference remains
`2.78e-3 K`; endpoint field gates still fail. Inlet-adjacent and missing
same-level stencils retain PCM; see the
[halo interpolation record](validation/eb-mol-linear-halo-20261019.md).
Ordinary held MOL transport now uses the same prepared, spatially interpolated
neighborhood as hydro in serial/MPI
([halo record](validation/eb-mol-transport-halos-20261020.md)). Prepared transport
also reconstructs open-face centroid fluxes before conservative species/enthalpy
limiting and uses both transverse stencils at patch interfaces. Analytic
transport checks and hot reaction/changed-rank restart pass. At 300 ns,
x-momentum error falls to `3.20e-5 kg/m2/s` and temperature error to
`1.97e-3 K`, with 19 reference field limits failing at that stage. See the
[transport correction record](validation/eb-transport-centroid-20261020.md).
Ordinary reacting NoRedist MOL now preserves fine-level energy after parent
re-reflux, validates old-clock states without redundant energy resets, and
resets only the closing level in serial/MPI. Rotating held-energy protection,
MPI partial recovery and hot changed-rank restart pass. The 300 ns temperature
error is `1.10e-4 K` and total-energy error `0.0202 J/m3`; all first-step field
limits pass. Twelve density/momentum/energy limits still fail at later steps,
so nonuniform EB-AMR parity remains open. See the
[energy synchronization record](validation/eb-energy-sync-20261020.md).
Covered coarse interpolation donors now also receive their closing energy
reset after restriction, without a composite inventory charge. Their energy
gap falls to roundoff and focused serial/MPI/restart checks pass, but the
twelve reference field failures are nearly unchanged; this was not the main
remaining discrepancy. See the
[covered-parent record](validation/eb-covered-energy-20261020.md).
A controlled native spatial-EOS intervention now identifies the dominant
remaining gap: using actual UEINT for both hydro pressure and transport
temperature brings every field within the original limits through 300 ns.
This diagnoses the energy-policy difference; the unmodified pinned reference
still fails twelve limits. Production retains held-energy cancellation safety.
See the [policy isolation record](validation/eb-held-policy-20261020.md).
The matched-policy diagnostic also completes 50/25 ns time refinement and
100 ns stepping to 1.2 microseconds with closed conservation budgets.
Observed time orders are about 1.1--1.2 in both solvers; second-order temporal
accuracy is not established. O2/H2O differences exceed the unchanged field
limit from 1.0 microsecond, so longer trajectory agreement remains open.
See the [temporal study](validation/eb-temporal-study-20261020.md).
A separately generated native activation-basis mechanism now removes those
six longer-trajectory O2/H2O exceedances against the held-EOS diagnostic:
all thirteen times through 1.2 microseconds satisfy the field and conservation
limits. The canonical Cantera mechanism is unchanged; unchanged-PeleC parity
and convergence remain open. See the
[activation-basis record](validation/eb-native-activation-20261020.md).
The variable-grid probe now also passes a 16-cubed root/fine-patch comparison
through 100 ns, with all 8192 initial states matched and closed budgets.
V1/V2 8-cubed results retain byte identity. This is a finer-grid integration
check; a common-time-step three-grid study is still needed for convergence.
See the [mesh pilot](validation/eb-mesh-pilot-20261021.md).
Conservative fluid-volume comparison now separates initial, endpoint and
evolution differences across grids. The common-50-ns 8/16 pair is complete;
the 32-cubed run finished but fails twelve field comparisons at 50/100 ns,
despite closed budgets. Stage captures locate amplification of small
water-profile differences by the shared coarse-to-fine slope limiter.
The original mass-fraction tolerance is scaled by local density; explicit
V3 density controls now bypass that conversion while retaining V1/V2 defaults.
The strict coupled pilot timed out; exact-input initialization caching removes
redundant solves. A subsequent coarse-corrector stall is now reproduced and
fixed: negative frozen species sources no longer impose nonnegative CVODE
internal substeps, while signed-trial bounds and endpoint rollback remain.
The analytic unit, 48 forced-cell solves and 8-cubed coupled check pass.
The 32-cubed five-evaluation prefix now completes in 734 s. Its water-limiter
gap shrinks, but momentum-controlled halo amplification remains; saved-state
caloric checks identify temperature-inversion residuals as the next target.
An explicit caloric convergence mode now preserves the default pinned inverse
while allowing tighter physical and signed-trial temperature recovery, with a
distinct restart policy identity. A private five-evaluation intervention reduces
the fine-halo energy gap from 8.92 to 0.00975 J/m3. The formal opt-in path passes
analytic/rollback checks and the 8-cubed coupled comparison; default-mode CSVs
retain byte identity. Matching 16/32-cubed trajectories now also pass through
100 ns. Evolution-increment orders are about 1.09--1.52 for the main gas
quantities in both solvers; second-order accuracy is not established.
The frozen 33-point table explains the initial temperature order of 1.012.
A separate 129-point cosine sampling study also passes all three grids.
It gives initial order 1.997 and endpoint temperature/density/H2O orders
near two, but main gas evolution-increment orders remain about 1.1--1.5
in both solvers. A regional saved-field audit does not confine the error to
the wall. The ordinary EB-tree namelist now exposes the opt-in caloric mode.
Its selected-mechanism gas application check passes serial/MPI field identity,
2-to-3-rank restart, and changed-precision/mixed-rank policy refusal; see the
[application record](validation/eb-eos-application-20261021.md). See the
[three-grid record](validation/eb-three-grid-20261021.md),
[caloric convergence record](validation/eb-eos-convergence-20261021.md),
[forced-CVODE record](validation/eb-forced-cvode-20261021.md) and
[restriction record](validation/eb-mesh-restriction-20261021.md).

Ordinary held MOL supports source-frozen replay of complete root and
sub-parent windows through the serial/MPI recursive driver. Explicit spatial
weights reconstruct gas and actual UEINT faces/rates together. Root and
sub-parent fixtures pass weights 0, 0.5 and 1 with exact serial/MPI outputs,
independent inventory checks and rollback. Native histories now retain actual
ancestor clocks as well as gas/T/UEINT/source endpoints. Sub-parent replay
restores pending stock and open parent face families; unit weight matches the
next original native clock entry. See the
[sub-parent record](validation/native-mol-subparent-replay-20261012.md),
[weighted-root record](validation/native-mol-weighted-root-replay-20261011.md)
and preceding [unit-root record](validation/native-mol-root-replay-20261011.md).
Native gas/UEINT source admission now spends current physical stock rather
than an unobserved future coarse/fine reflux credit. The previously refusing
future-fine producer passes its full-mechanism one/two-correction fixture in
both owner layouts, with exact pending gas/H payment and independent inventories;
see the [current-stock funding record](validation/native-mol-future-source-funding-20261012.md).
The [trace-debt record](validation/native-mol-trace-debt-admission-20261012.md)
retains the preceding donor-guard evidence. Automatic source-frozen weighted
recovery of completely observed ordinary held MOL parent windows now has
[natural closure and rollback evidence](validation/native-mol-automatic-parent-20261014.md).
Automatic recovery also carries resolved gas/UEINT receipts through two
fully observed sub-parent windows in a controlled three-level fixture; see the
[material automatic recovery record](validation/native-mol-automatic-material-20261014.md).
Selected ordinary native MOL common steps can also retain their observed
external and gas/UEINT auxiliary packets after a native failure, while declined
caller packets and native-admitted reactor refusals publish no recovery prefix.
See the [source-prefix capture record](validation/native-mol-source-prefix-20261014.md).
Explicit common-clock reuse now consumes observed packets, rebuilds native
spatial operators and evaluates newly reached source stages once; see the
[prefix consumer record](validation/native-mol-prefix-consumer-20261015.md).
Recursive ordinary native MOL common clocks now retry failed spatial families
while retaining observed source packets, with bounded retries and transactional
rollback. See the [recursive retry record](validation/native-mol-recursive-prefix-retry-20261015.md).
A newly recovered partial single-root clock now passes exact serial/MPI and
complete-frozen-replay comparison, independent gas/H inventories and refusal
rollback for two owner layouts; see the
[partial root recovery record](validation/native-mol-partial-root-recovery-20261016.md).
A recovered partial root now also completes both child clocks with exact
serial/MPI states and faces, independent composite gas/H inventories and
late-child whole-transaction rollback; see the
[descendant record](validation/native-mol-partial-descendants-20261016.md).
Partial fine-clock recovery under an open root now also passes mixed-leaf
serial/MPI inventories, exact explicit-prefix reproduction and late-child
rollback; see the [fine-clock record](validation/native-mol-partial-fine-20261016.md).
Partial sub-parent recovery with descendants and opt-in continuation from a
complete-clock prefix now have exact serial/MPI, inventory and rollback
[evidence](validation/native-mol-window-prefix-20261018.md). Nonzero-material
partial recovery, broader window replacement and source/physical qualification
remain open.

Explicit full-mechanism root and first-fine replay now pass material-receipt
cases at weights 1, 0.5 and 0 in both owner layouts. Native adjacent donor
queries now enter collective admission on all ranks, with owner-local
publication and idle-rank refusal rollback. See the
[material replay record](validation/native-mol-material-replay-20261013.md).

The pinned 53-species n-dodecane bundle now imports Simple transport fits.
All 6,148 coefficients match compiled upstream output; four transport states
agree within `1.73e-15` relative. This prepares the transport dependency for
Spray-A-Wbreakup
([import evidence](validation/dodecane-simple-import-20260927.md)).
A 53-species Spray-A adaptation now completes its first two injection events
with variable liquid properties and KHRT. Fuel accounting and a breakup-off
intervention pass; full-duration/reference comparison remains open
([case and evidence](validation/spray-a-connection-20260927.md)).
The KHRT injection-reference normalization now includes the pinned 0.1
multiplicity. A corrected eight-event Release run creates 1,232 new particle
IDs, finishes with 1,264 parcels and conserves fuel within `3.13e-16` relative
([correction and child release](validation/spray-a-khrt-threshold-20260927.md)).

An extracted Jet failure revealed internal-energy inversion stagnation at
adjacent representable temperatures near zero energy. The inversion now
terminates on that resolved bracket; focused thermodynamic checks and a bounded
Jet replay through 0.430 ms pass. Full prepared 1.5 ms continuation is running
([evidence](validation/temperature-inversion-resolution-20260927.md)).

Rotating-frame momentum/energy kernels and optional MOL energy-flux correction
now have analytic checks and a 60-case compiled-reference comparison
(maximum scaled error `4.37e-15`; [kernel evidence](validation/rotating-frame-kernels-20260927.md)).
EOS recovery, wall pressure and single-patch MOL now accept rotating energy.
A bounded 400 rad/s patch check passes; a 1000 rad/s LF condition exceeds the
existing absolute conservation gate ([patch evidence](validation/rotating-mol-patch-20260927.md)).
Serial tree initialization, conserved transfer, halos and MOL spatial rates now
carry the frame ([tree evidence](validation/rotating-eb-tree-transfer-20260927.md)).
Serial inviscid, nonreacting MOL now integrates rotation on common or recursive
time steps with bounded composite conservation checks
([SDC evidence](validation/rotating-sdc-20260927.md)). Forced chemistry now carries
the rotating energy convention, including density forcing
([reaction evidence](validation/rotating-reaction-20260927.md)). Viscosity, heat/species
transport and embedded-wall sources now connect to serial rotating SDC
([transport evidence](validation/rotating-transport-20260927.md)). MPI gas
transfer, MOL stepping, transport and SDC now carry the frame. Serial/MPI gas
applications accept rotating MOL input, including inlet, automatic timestep,
scheduled regrid and checkpoint continuation across rank counts. MPI EB spray
now converts particle and film exchange, regular-cell LES states, and the
application's particle/film inventory and injection ledger to rotating energy.
Film SDC stages now reject rank-disagreeing or mid-interval changes to the
rotating-frame settings before using cached source rates.
Particle outflow energy now uses the potential at the EB flight exit; direct
drift and collective impact checks cover this boundary ledger conversion.
Two recursive two-rank fixtures and bounded two-step spray application runs
with and without LES pass. Particle-trajectory and fixed-reference boundary-ledger parity,
broader conservation qualification, and fixed-reference application comparison
remain open.

Regular serial/MPI AMR MOL now accepts rotating slip, no-slip and isothermal
domain walls and prescribed-pressure outflow on common and recursive intervals.
The pressure-boundary EOS receives ordinary energy; modified ghosts then recover
rotating total energy while retaining ordinary internal energy. Untouched cells
retain their exact payload. Mach-inlet internal energy is now computed before
the rotating total-energy conversion. These connections do not establish fixed
PeleC boundary-trajectory parity. The focused wall/pressure matrix passes
common/recursive 2-to-3-rank restart in eight conditions; the three focused
boundary and internal-energy unit tests also pass.

Regular MPI AMR spray now selects MOL through both fixed and selected-mechanism
frontends, forwarding LF flux selection, molecular transport, LES, external
forcing and common/recursive particle/film sources. Its particle stages follow
the pinned old/new endpoint source construction, including reacting corrections
that do not rebuild the new particle family. Common MOL supplies parent endpoint
halos with one interval per level. The library also forwards the rotating frame
to particle/film sampling and source conversion, recovering gas temperature from
an ordinary-energy copy after physical exchange correction. The remap library
preserves rotating conserved energy before recovering temperature; its focused
coarsening/owner-migration and frame-mismatch rollback test passes, including
particle and reaction/heat history migration. The regular spray application now
accepts `reactive_rotation` input, forwarding the frame through MOL/LES,
remapping and automatic timesteps. Native rotating energy is retained in
checkpoints; EOS validation and gas CSV output use an ordinary-energy copy.
Frame parameters and cell coordinates bind checkpoint identity. Fixed and
selected LF/LES two-to-four-rank continuation fixtures pass, including cumulative
coordinate energy, scheduled regrid/injection, changed-frame rejection and
nonfinite-history rejection. The library now exports
`frame_energy_change` separately from physical particle-force work. It accounts
for the gas source's cell-centered potential, liquid mass and position changes,
and the actual exterior exit potential, counting only physical particles and
films. The strict rotating moving-particle exchange test passes without changing
its `1e-19` tolerance or adding rotating forces to the particle update. Analytic
checks cover signed deposited sources, covered-cell/halo exclusion, pure drift,
exterior exit energy and transactional rejection. Recursive EB particle/film
wrappers expose the same ledger. The EB application reports a separate
`coordinate_energy` budget column and retains its cumulative value in the coupled
checkpoint. Its fixed-mechanism two-step changed-rank continuation and pre-ledger
checkpoint rejection fixture previously passed. Native forced Newton trials now
resolve coupled chemical intermediates below the absolute tolerance, using the
largest component for mass closure and a stiffness-aware correction criterion
within the existing bulk residual gate. The water/N2 trace balance rejects the
preceding solver and passes the current one. The actual MPI rotating reaction
fixture and selected regular LF/LES changed-rank restart also pass; the reaction
fixture's mode buffer now retains its full 17-character name.
Forced chemistry now also accepts the held old internal-energy scalar as its
initial energy, matching the pinned `React.cpp` packing of `UEINT`. External
energy forcing and the published endpoint still use the total/momentum
projection. Normal MOL forwards this scalar; cell/EB native trajectory and
rollback checks, and exact serial/MPI owner/idle reaction checks pass. Common EB
MOL now optionally returns a separate composite numerical-energy correction for
old, half-weighted predictor, final and post-sync UEINT resets. A nonreacting,
inviscid three-level forced-kinetic fixture checks its independent inventory,
exact owner/idle gas/T/scalar faces and atomic reset/policy failure. Recursive
MOL now carries this correction privately through canonical clock anchors and
leaf synchronization, publishing it only with the outer transaction. A focused
three-level no-redistribution fixture checks updates `[1,2,4]`, a compensated
independent inventory, exact owner/idle gas/T/held and registered/physical/scalar
faces, and atomic geometry/policy/prefix/presence rejection. Gas and held MOL
faces share one native stage evaluation; the same-level conservative correction
carries a signed held face/rate with the gas donor coefficient while preserving
independent volume terms. Common no-redistribution MOL now connects the held old
energy, admitted held spatial mean, native molecular transport and forced
chemistry. Its detailed three-rank fixture now passes the independent inventory,
exact distributed/owner-zero gas/T/held faces and correction checks, and native
reaction-refusal rollback (956.60 seconds). Selected common intervals retain
owner/level-sparse old and accepted held energy, spatial rates/faces, reset policy
and numerical correction. Manufactured paired face changes and native pure-N2
fixed-mechanism cases pass on two/four ranks, including capture/reconstruction,
atomic rejection and removal of stale held payload on gas-only replacement.
Ordinary recursive MPI MOL now archives the full gas/T/reaction/held clock and
numeric correction prefix. No-redistribution held interval replay reconstructs
both face families with the same interval weight, restores the clock prefix and
recomputes endpoint/leaf resets without charging the previous descendant tail.
Two/four-rank native pure-N2 cases with one/two corrections pass full clock
capture, independent inventory, idle/owner policy rollback and a native-admitted
child provider refusal after the parent completed. The latter retains all input
fields, the numeric prefix and prior accepted clock/face payload without retry.
The recursive driver now uses a shared atomic native held-interval replay
operator. Two/four-rank tests replay independently signed gas/held face impulses,
checking projected gas energy forcing, the held old UEINT anchor, distinct
chemical/numeric ledgers, complete clock/heat history and native-provider refusal
without publication. Ordinary interval archives now retain initial reaction
history too, and replacement clocks are captured after their canonical start
repair. Common prepared MOL now retains each independent
held old/raw/spatial source, native/low face family and scalar register stage.
A three-rank two-correction test passes paired native reconstruction, weights
1/0.5/0, owner-zero/idle layouts, remote NaN and stage-control refusal, and
native-admitted provider rollback (72.22 seconds). Recursive prepared held
chemistry now restores complete scalar clock/prefix/open parent registers and
uses a private paired gas/scalar closure. A one-correction three-level diagnostic
passes `[1,2,4]` clocks, exact owner-layout gas/T/held endpoints, provider rollback
and independent inventory `1.97e-16`. A two-correction nonuniform diagnostic
initially failed its unchanged inventory gate (`8.12e-6`); isolated provisional
coarse-stage registers explain the energy residual `-0.84317331168`. Held
reacting prepared intervals now publish the actual old/final mean register
family while preserving every intermediate native stage separately in the
archive. A collective `accepted_coarse_stages` contract governs reconstruction
and replay; mixed owner/idle contracts reject without publication. The same
two-correction fixture now passes its unchanged `3e-14` gate at `2.56e-16`, exact
owner-layout gas/T/held endpoints and native-provider rollback (108.68 seconds).
Common weighted replay and legacy gas-only interval history also pass (73.95
and 15.71 seconds). This accepted coarse-stage policy differs from the pinned
direct provisional coarse accumulation, which remains in the legacy gas-only
path. It is a conservation fix, not established pinned trajectory parity;
that numerical comparison remains required. Serial registration reconstruction
and the MPI `register_provisional_coarse` selector now assemble the source-policy
registration independently: all admitted coarse evaluations in chronological
order, but only old/final fine dm, with the stored half-scale applied once.
The nonuniform native fixture has a nonzero provisional scalar coarse delta.
Its focused three-rank tests pass exact serial/MPI registered gas/UEINT deltas,
unchanged physical RHS, owner/idle policy refusal, late malformed-stage and
partial-allocation rejection with output rollback (39.64 and 83.20 seconds).
The five local source files defining this schedule match byte-for-byte the
pinned PeleC and AMReX revisions. This is registration/source-policy evidence;
the current coupled state synchronization still uses accepted old/final deltas.
Integrating the complete pinned coarse registration with the state-update and
conservation requirements, and pinned numerical trajectory comparison, remain
open. Source/binary identity, preceding passes and commands are retained
in `build/eb-held-mol-registration-20261006`. Naturally failed prepared held
clock-group replay remains unqualified. Serial recursive no-redistribution MOL
now connects held forced chemistry and native transport through `[1,2,4]`
clocks, forwarding embedded-wall controls. A three-rank two-correction pure-N2
gradient fixture passes independent serial/MPI inventory (`1.90e-15` and
`3.64e-16`), exact gas/T/held/IR/heat and physical/registered/scalar face families,
and native provider-refusal rollback in serial, distributed and owner-zero/idle
layouts (77.81 seconds). Serial ordinary no-redistribution held MOL now retains
complete gas/T/UEINT/IR/heat clock starts and numeric prefixes. Native paired
interval/group replay rebuilds gas forcing and the independent scalar family,
using one spatial weight and the first old face anchor on each level. The
recursive driver now invokes this private full-window replay for failed gas,
scalar or energy-reset closure, and can defer a completed final child window to
its enclosing parent. Supplied-provider refusal on a native-admissible input
aborts without weight retries. A three-rank pure-N2 test passes weights 1/0.5/0
with independent inventory residuals `1.88e-15` to `1.92e-15`, `[1,2,4]` physical
clock metadata, replacement of fabricated field/correction tails, complete
provider/clock-gap rollback and exact heat against independent native reaction
calls (78.80 seconds). It also retains the native serial/MPI endpoint/face
comparison and serial archive refusal checks. Its built source/test/binary
identity and logs are in `build/eb-held-serial-parent-replay-20261006`; the
initial exact-zero inert-heat assumption failed and is not counted as a pass.
Ordinary reacting held MOL now accepts frozen auxiliary gas sources with zero
UEINT volume source in serial and MPI, preserving the actual density/species
volume source during replay. A controlled covered-coarse H2O2 injection with
the pinned full mechanism causes a real native chemistry/reflux
closure failure. Its three-rank diagnostic passes automatic replay of all seven
physical intervals at spatial weight zero, exact serial/MPI gas/T/UEINT/IR/heat
and physical/registered/scalar faces, and exactly zero composite reactive
inventory. Independent inventory residuals are `1.20e-14` (serial) and `5.20e-15`
(MPI), within the unchanged `3e-14` gate. Source refusal, including an idle rank,
and native-admitted reactor refusal inside actual automatic replay preserve
caller fields, ledgers and accepted clock/scalar-face payloads without retry
(181.00 seconds). Built identity and logs are retained in
`build/eb-held-aux-natural-replay-20261006`. This controlled fixture qualifies
ordinary no-redistribution replay, not physical validation or pinned trajectory
parity. Serial prepared held
MOL now connects native forced chemistry, transport and gas/UEINT StateRedist.
Its old/final source, face and register families share one conservative spatial
weight; local reaction/EOS admission and supplied-provider refusal are private
until accepted. Serial common/recursive and MPI common/recursive prepared
closures now assemble matching gas and scalar rereflux before leaf energy reset
and average-down. Two focused three-rank tests pass (93.18 and 190.22 seconds),
including exact serial/MPI gas/T/UEINT/IR/heat and paired face families, provider
refusal rollback and independent common/recursive inventory. Serial inventory
residuals are `1.87e-15` (common) and `1.74e-15` (recursive); MPI residuals are
`1.76e-15` and `2.56e-16`. Source/binary identity and completed logs are retained
in `build/eb-held-serial-prepared-20261006`. Serial prepared intervals now retain
each native gas/scalar stage and replay with frozen external volume rates,
rebuilding StateRedist with one prescribed spatial weight. The three-rank
weighted comparison now passes with nonzero force at weights 1/0.5/0, including
exact gas/T/UEINT/IR, paired raw halos/faces/registers, malformed scalar-stage
rejection and native-admitted provider-refusal rollback (102.13 seconds).
Recursive serial archives now retain all seven complete clock starts and open
ancestor gas/scalar registers; exact serial/MPI clock/register and existing
common/recursive inventory checks pass (199.21 seconds). Built identity and
completed logs are in `build/eb-held-serial-prepared-replay-20261006` under
`sparse-fixed-identity.json` and `sparse-fixed-tests.log`. These fixed elementary pure-N2 diagnostics are not pinned trajectory
parity or reacting-mixture validation. Serial prepared full-window replay now
uses the same recursive physics driver, restoring complete clocks, ancestor
trajectories and open gas/scalar parent registers. Its focused three-rank test
passes weights 1/0.5/0 (`1.74e-15` to `2.17e-15` inventory residual), exact native
unit-weight endpoints/faces, an internal parent window and malformed-time/
register rejection before chemistry. A native-admitted child-provider refusal
after completed root physics preserves every field, ledger and prior clock
without retry (96.84 seconds). The internal-window comparison applies an
independent native next-clock reset; direct comparison across that boundary
initially failed and is not counted as a pass. Built identity and completed logs
are in `build/eb-held-serial-prepared-parent-20261006`. Prepared held intervals
now archive and replay actual auxiliary gas rates with zero UEINT volume source,
without repeating source callbacks. The frozen callback-dispatch diagnostic
passed its no-source recursive test (45.64 seconds). Its reacting natural replay
accepted serial and both MPI layouts, but timed out at 1800 seconds before the
refusal checks finished; that test is not passed. This executable predates the
subsequent direct held-energy EOS and transfer corrections; its identity and
terminal log are retained in `build/eb-held-prepared-aux-20261006`.
Explicit paired auxiliary sources now use an opt-in `Nv+1` packet: gas rows
retain their layout and the last row carries actual held UEINT, its independent
volume rate and physical source budget. Serial/MPI common and recursive MOL
forward this family, retain frozen old/predictor packets and route off-level
paired impulses through private gas/scalar closure and numerical reset accounting.
A focused three-rank fixture passes exact serial/MPI nonreacting and reacting
routed endpoints, idle/NaN/bounds refusal rollback, signed formation-energy input,
zero-weight common prepared replay, serial full-window prepared routed replay
without repeating source callbacks, and no-ledger routed temperature recovery with
independent held ancestor/time-fraction input (23.54 seconds). Its pure-N2 analytic forced
reactor isolates source/replay structure; it does not qualify chemical integration.
The existing gas-only routed-history fixture also passes (20.01 seconds).
Source/binary identity and logs are in `build/eb-held-auxiliary-family-20261006`.
A selected real-reactor recursive StateRedist application check also passes exact
serial/MPI/restart fields and independent inventory (`2.165e-15`); it uses no
auxiliary source and predates the final optional-ledger/halo-entry guards. Its
actual binary/command hashes and preceding source snapshot are retained separately.
MPI failed-closure prepared routed replay, physical particle/film source consumers
and pinned source/trajectory comparison remain unqualified or unconnected.
Native root initialization now forms UEINT from supplied temperature and
composition before recovering gas primitives, and carries it through paired
PCM prolongation. Coupled volume-weighted restriction, serial/MPI tree averaging
and remap retain UEINT alongside gas/T and reaction/heat fields through owner
changes. Tree halos now use the held scalar in donor and peer EOS recovery,
including width-one faces, wider target metrics, wall reflection, prescribed
and Mach inlets, and ancestor time interpolation. Fixed elementary-mechanism
two/three-rank diagnostics pass high-kinetic ordinary/rotating cases, exact
serial/MPI fields, independent remap inventories within the unchanged `3e-14`
gate, and atomic EOS/presence refusal. Build/source identities and logs are in
`build/eb-held-tree-{average,remap,seed,halo}-20261006`. Root timestep estimation
now forwards UEINT through acoustic, molecular, halo and isothermal-wall bounds.
Two/three-rank tests pass independent NASA/geometry bounds for common and
recursive coupled/split clock policies, exact serial/MPI timesteps and atomic
scalar/presence refusal (`build/eb-held-tree-timestep-20261006`). Patch checkpoint
schemas now retain independent UEINT, signed cumulative numerical-energy
correction and all reset controls, binding those controls on reads. Fixed/selected
split/SDC records pass high-kinetic ordinary/rotating and native negative-formation
energy roundtrips with policy, scalar, prefix, trailer and presence rollback
(`build/eb-held-patch-checkpoint-20261006`). Serial V7 and indexed MPI V8 tree
envelopes retain the same fields and account for the composite correction once;
each embedded patch prefix must match the envelope. Two/three-rank diagnostics
and an actual two-rank write followed by three-rank read pass exact gas/T/UEINT,
reaction/heat data, serial/MPI interchange, idle policy/presence refusal, late
scalar/prefix rollback and preservation of an accepted file after writer refusal
(`build/eb-held-tree-checkpoint-20261006`). These are library serialization and
transfer checks, not full application continuation or pinned trajectory parity.
The opt-in gas application now binds UEINT/reset controls to initialization,
common/recursive MOL steps, timestep selection, scheduled regrid, cost-based owner
migration, serial/MPI checkpoint continuation and CSV/numerical-energy budgets.
Fixed gas/common and selected full-mechanism chemistry-enabled pure-N2 recursive application
checks pass three steps with collapse/recreation of a child, exact gas/T/UEINT CSVs,
serial/MPI interchange, two-to-three-rank restart, exact numerical-prefix restore
without advancing, control mismatch refusal and protected outputs. Independent
composite inventory residuals are `1.09e-15` and `1.64e-15`; MPI summation can differ
from serial by one final bit, so continued prefixes are compared with the native
destination run. Owner migration and CSV tests also pass two/three ranks with
independent signed UEINT and late scalar/idle presence refusal
(`build/eb-held-app-binding-20261006`). The selected chemistry-enabled pure-N2
common StateRedist lifecycle also passes exact serial/MPI/restart fields,
numerical-prefix/control refusal and regrid (`2.43e-15` composite residual).
Native MOL cell recovery, cut-wall pressure
sources, common held reflux/reset and average-down now carry actual UEINT;
MOL/wall unit checks pass independent NASA/pressure bounds, unresolved rotating
thermal energy, signed formation energy and late scalar rollback
(`build/eb-held-pressure-wall-20261006`). An extreme `omega=1e10` hydro-only
application diagnostic passes the same lifecycle/rank-changing checks with
`2.80e-16` inventory residual. Molecular cell, physical-donor and isothermal-wall
EOS recovery now accepts actual UEINT, including serial/MPI narrow and prepared
transport halos and quadratic embedded-wall library donors. Independent Fourier
flux and late finite-scalar rollback checks pass. The same extreme-rotation
common-step diagnostic with all molecular transport passes exact gas/T/UEINT
output, serial/MPI and two-to-three-rank restart, numerical-prefix restore and
protected-output refusal (`2.80e-16` composite residual; chemistry disabled).
Recursive entry/exit restriction now averages gas/T/UEINT together, and ordinary
no-redistribution synchronization closes both flux families before reset/EOS.
The extreme chemistry-enabled recursive diagnostic now passes the same exact
serial/MPI fields, regrid, changed-rank restart, prefix and refusal checks
(`2.80e-16` composite residual). The chemistry-enabled recursive StateRedist
lifecycle also passes at this rotation (`4.04e-16`). Nonreacting recursive
StateRedist now passes one/two SDC corrections at this rotation, retaining the
old-half registered gas family separately from the accepted physical transfer
(`5.27e-16`, exact native/rank-changing restart fields). Source identity, failed
attempts and completed logs are in `build/eb-held-molecular-eos-20261006`;
`paired-chemistry-identity.json` binds its preceding phase, and `identity.json`
binds the completed nonreacting extension. Shared gas applications now connect
actual UEINT to the embedded-wall timestep, retaining native thermal/viscous
sources and wall/reset budgets. Independent quadratic heat/normal-stress and
late finite donor rollback checks pass. Common and reacting recursive
StateRedist two-correction wall lifecycles pass exact gas/T/UEINT output,
serial/MPI interchange and two-to-three-rank restart at 400 and `1e10` rad/s;
composite residuals are `2.01e-15` and `1.16e-15`. The MPI heat-only API also
accepts the held/ancestor family. Its two/three-rank owner/idle tests pass
independent composite heat/traction, exact fields and remote scalar refusal
without publishing rates or ledgers. Source/build identity and completed logs
are in `build/eb-held-embedded-wall-20261006`.
Held prepared StateRedist MOL now also connects spatial LES callbacks in common
and recursive serial/MPI gas steps. Wide filter donors carry the same held/ancestor clock;
LES itself retains the total-energy EOS and UEDEN-only diffusion source seen in
the pinned code. A classified MOL source is evaluated at old/predictor times;
later reacting corrections retain its old half rather than rebuilding new LES.
This policy is retained in interval records and checked during replay. The
actual seven-clock reacting Smagorinsky diagnostic records exactly seven old
and seven predictor callbacks, with no second-correction callback. Common
Smagorinsky and reacting recursive Dynamic two-correction lifecycles pass exact
serial/MPI gas/T/UEINT, owner migration, scheduled regrid, changed-rank and
serial/MPI restart, independent held/T/inventory checks and protected-output
refusal (`1.60e-15` and `1.73e-15`). Smagorinsky changes the forced reflecting-wall
state relative to disabled LES. Changed LES enable/coefficient controls reject
without altering CSV, budget or checkpoint. Serial and MPI now assemble prepared
molecular, MOL and LES contributions in the same order; an earlier momentum-only
`6.94e-18` difference was corrected without loosening the exact output gate.
All-positive nodal cells/faces now retain exact regular volume/aperture/center
identities rather than rounded tetrahedral sums. Independent identity and
legacy cut/coarsening/connectivity checks pass. Source/build identities and
completed logs are in `build/eb-held-les-binding-20261006` and
`build/eb-regular-nodal-identity-20261006`. These short native diagnostics are
not pinned solver trajectory parity or physical LES/rotating-flow/ignition
validation. LES temporal/CF registration and broader model/reference cases still
require the final qualification.
Ordinary no-redistribution MOL now also connects held spatial diffusion
callbacks on common and recursive serial/MPI steps, using actual paired
six-layer donors and the same molecular/MOL/additional-source order. MPI MOL
LES uses the same spatial evaluator as serial; its valid-only callback API also
accepts UEINT and forwards the callback and source-stage policy to child clocks.
Common Smagorinsky, reacting recursive Dynamic, and reacting recursive
Smagorinsky with molecular transport disabled pass the short lifecycle checks:
exact serial/MPI/restart gas/T/UEINT, regrid/migration, independent inventory and
protected-output refusal (`1.60e-15`, `1.73e-15`, `1.58e-15`). Logs and identities
are in `build/eb-held-diffusion-none-20261006`; this is native integration
evidence. The three-rank valid-only callback unit passes (106.12 seconds) with
three-level common/recursive clocks, generic/classified stage counts, molecular
transport on/off and two owner layouts including an idle rank. It checks exact
UEINT fields/faces against a native no-callback baseline, nonzero additional gas
faces, an actual old-energy reset and late callback refusal preserving every
field, face and ledger. Its pure-N2 reaction provider uses the exact analytic
forced inert solution to isolate this interface; the preceding full native
reaction variant timed out at 180 seconds and is not a pass.
The Godunov raw patch operator now captures actual UEINT Riemann fluxes through
normal and transverse CTU stages, applies viscosity to the independent held
scalar and uses the shared centroid stencil. It returns advection minus pinned
regular/cut-cell pressure work, integrated and low-order scalar faces, and
paired forward-flux correction while preserving the pressure-work residual.
Gas/scalar outputs publish together; incomplete inputs and late invalid scalar
data retain prior rates, faces and ledgers. Physical boundary wrappers forward
actual held donors and accept optional independent inlet UEINT. The primitive
forcing continues to use the pinned `pc_srctoprim` total-energy projection.
Focused unit checks pass regular/three-axis cut flux inventories, independent
held viscosity, nonzero gas forcing, exact gas-family identity, six supported
physical boundary kinds and active paired face limiting (11.88 seconds); the
legacy PCM/PLM viscosity unit also passes (8.31 seconds). The local `Hydro.H`
matches the official frozen commit byte-for-byte. Its actual 3-D `pc_eb_pdivu`
function matches the Fortran helper exactly on 192 regular, inclined-plane,
spherical and covered fixture cells. Fixed/selected serial/MPI gas frontends
compile; source/binary identities, commands and failed probes are retained in
`build/eb-held-godunov-spatial-20261006`. This is lower-operator/kernel evidence:
Serial/MPI Godunov tree source APIs now carry actual held halos, ancestor
endpoints, raw RHS/old populations, integrated/low-order scalar faces and
pressure work on common or selected levels. Peer recipients use exchanged
valid scalar rates. The prepared StateRedist API applies the pinned factor-one
registration to the independent UEINT component; it retains raw populations
separately and publishes gas/scalar fields, faces and registers together.
A three-rank focused check passes regular/cut three-level peer trees, two owner
layouts with an idle rank, selected-level gas/source/UEINT ancestor clocks,
exact serial/MPI fields/halos/faces/registers, independent raw inventory
(`6.94e-18`), and collective idle/partial-family/late-NaN/zero-based scalar
refusal with every output retained (184.52 seconds). The existing gas-only
CTU tree unit also passes (31.22 seconds). Fixed/selected serial/MPI gas
frontends compile. Identities, commands and the diagnosed ancestor-clock fixture
failure are in `build/eb-held-godunov-tree-20261006`.
Common and selected-level SDC now connect held Godunov rates and faces, actual
UEINT diffusion/EOS halos, independent scalar StateRedist stages and the paired
forward correction. A successful step can return the signed composite pressure
work before redistribution. Its scalar registrations follow the gas coarse
stage chronology and final fine weights; all outputs publish with the accepted
state. A nonreacting inviscid three-level cut-tree diagnostic passes exact
serial/MPI state and face comparisons, two owner layouts with an idle rank,
collective query refusal and independent UEINT inventory (`5.91e-16`). Composite
reductions differ only by a few ULP; field and face equality remains exact.
The common, selected-level and native-diffusion three-rank checks pass
(60.63, 26.40 and 57.15 seconds), including exact scalar/gas registers and
retention of all prior fields, faces, registers and diagnostic queries on
collective presence mismatch. Exact identities and results are retained in
`build/eb-held-godunov-step-20261006`.
Recursive serial/MPI stepping now forwards actual held ancestors, uses the
Godunov factor-one scalar registration and accumulates signed pressure work
privately over accepted level clocks. The native recursion had returned only
immediate-child faces; it now retains the complete descendant gas, physical and
scalar interval family, with matching prepared MPI replay return. A three-level
cut-tree nonredistributed diagnostic passes `[1,2,4]` clocks, two owner layouts
with an idle rank, exact gas/T/UEINT and face comparisons, query rollback and
independent scalar inventory (`1.61e-15`). The complete recursive no/StateRedist
and native-diffusion checks pass (165.86 and 157.24 seconds), as does the existing
MOL reset-policy check (19.44 seconds). Results and identities are in
`build/eb-held-godunov-recursive-20261006`.
Fixed/selected serial/MPI gas applications now enable held Godunov SDC and
report cumulative UEINT pressure work separately from the UEDEN residual.
Serial V9 and indexed MPI V10 checkpoint records retain this history; readers
keep legacy V7/V8 MOL continuation while Godunov requires the new history.
A selected pure-N2 native-reactor/reset common StateRedist lifecycle passes
regrid, owner migration, two/three-rank and serial/MPI continuation, exact
physical CSVs and scalar prefix restore, protected policy-refusal outputs and
nonfinite pressure-history rejection. Its independent gas inventory residual
is `1.64e-15`. A three-step serial recursive StateRedist native-reactor/reset
smoke also passes. Two checkpoint units pass; actual previously saved V7/V8
MOL read-only restarts retain exact gas/T/UEINT CSVs and zero Godunov pressure
work. Evidence is in `build/eb-held-godunov-app-20261006`. These are bounded
pure-N2 integration checks, not active-fuel chemistry or solver trajectory
parity. Godunov reacting/LES trajectories and paired reactive parent replay
remain unfinished. Full pinned solver trajectory parity remains open.
The native Godunov raw and prepared StateRedist source APIs accept one optional
interval face weight for gas and UEINT together. They blend integrated faces
before divergence and redistribution, keep native pressure work and embedded
wall volume sources independent, and retain exact PCM/native endpoints at
weights zero/one. Owner and idle ranks must agree on the weight and its presence;
nonfinite or out-of-range inputs retain all prior output families. This is the
lower-operator foundation for native reacting parent replay. Common and
selected-level Godunov SDC now forward this hydro face weight in every native
evaluation and paired forward-correction retry, while retaining the old/new
diffusion and external volume families. Each actual iterate recomputes native
pressure work; a frozen MOL mean cannot supply it. MPI checks weight presence,
bounds and owner/idle agreement before native providers or energy reset.
Focused three-rank common no/StateRedist SDC checks pass absent/one/zero/0.375
weights, exact native weight-one identity, two owner layouts with an idle rank,
independent UEINT inventory (`6.45e-16`) and full collective weight rollback
(243.10 seconds). A selected pure-N2 StateRedist native-reactor/reset check
also passes exact gas/T/UEINT, faces and registers, bounded reset/pressure
ledgers and retention of all prior publications on invalid weight/presence
(83.51 seconds). The initial reactor fixture supplied an empty reaction table;
it was corrected to load the real mechanism without changing solver guards.
These checks exercise the native call/reset path, not active-fuel chemistry.
Results and exact identities are in `build/eb-held-godunov-sdc-weight-20261006`.
Common/selected Godunov SDC now offers a success-only clock archive with actual
pre-reset/initial/accepted gas/T/UEINT, independent source stages and hydro
passes, predictor RHS, pressure work, gas/scalar faces and register chronology.
Serial and sparse MPI histories retain ancestor/source endpoints and classify
provider presence. A selected cut-tree native-diffusion/pure-N2 reactor check
passes exact stage data and preservation of the archive on idle weight/query
disagreement (24.80 seconds), with identities in
`build/eb-held-godunov-history-20261006`. Recursive Godunov stepping now connects
preorder clock collection, parent/subtree indices, pre-reset physical state,
pressure/reset prefixes, open parent registers and routed pending sources,
plus synchronized endpoint states and complete descendant faces. Publication
moves the complete family only after the outer transaction succeeds; an absent
common-step scalar output retains true optional forwarding. The focused recursive
native-diffusion check passes seven `[1,2,4]` clocks, exact owner-local histories,
parent/subtree/time metadata and synchronized descendant face payloads, with an
independent UEINT residual of `1.23e-15` (154.74 seconds;
`build/eb-held-godunov-clocks-20261006`). Its initial nonreacting archive query
exposed an unallocated serial reaction source; a canonical zero source fixes
that query without changing the non-query native path.
Common/selected actual-held Godunov replay now uses accepted auxiliary packets
and explicitly declared external provider packets while reevaluating native
hydro, diffusion, pressure work and StateRedist from fresh state. Physical
providers are reevaluated rather than frozen. A changed-state pure-N2/native
transport check passes exact gas/T/UEINT, faces and registers against an
independent live-provider calculation, two MPI owner layouts with an idle
rank, late cached-NaN rollback, physical-provider refusal and native-admitted
opaque reactor refusal without repeated calls (93.03 seconds;
`build/eb-held-godunov-source-replay-20261006`). Full-clock replay now connects
the fresh recursive native driver with whole-family cache/topology/owner and
parent/time/subtree preflight. Its first MPI runtime check exposed rank lookup
after cache validation; rank lookup now precedes that guard, and the test checks
acceptance before inspecting output. Actual off-level auxiliary impulses are
now archived after fluid/leaf filtering and replayed through the existing
pending buckets. Eligible reactive native steps privately retain their complete
clock family; after a refused closure with all native clocks accepted, a root
trial recomputes native physics at paired PCM/native face weights and uses the
unchanged acceptance gates to choose a conservative interval weight. State,
faces, records, counts and pressure/numerical ledgers publish together. These
connections build for all four fixed/selected serial/MPI gas applications.
Full-clock and routed-source replay checks pass (331.88 and 334.57 seconds;
`build/eb-held-godunov-root-repair-20261006`); they prove fresh native replay,
not natural automatic repair selection. Nonroot replay now restores the
complete pre-reset gas/T/UEINT/IR/heat clock, numerical prefix, unchanged
external ancestor/source endpoints, per-level clocks, open gas/scalar parent
registers and routed sibling buckets. It avoids repeating the outer global
reset, preserves global preorder metadata and stages parent return allocations
before publication. Its four application builds pass. The focused second-parent
group check passes `[0,1,2]` clocks, exact physical/stage/IR/face payloads,
independent open-parent registration and routed sibling buckets, and rollback
for native-scale/mask/IR/pending-NaN/idle-time failures (145.21 seconds).
Explicit numerical and absolute-pressure prefix comparisons also pass
(147.49 seconds; `build/eb-held-godunov-subtree-replay-20261006`). Automatic
parent repair now connects completed native clock groups to a fresh recursive
trial and replaces gas/T/UEINT/IR/heat, actual old/new source-predictor endpoints,
open parent registers, routed pending sources, faces and numerical/pressure
ledgers together. A classified provider refusal stops further parent trials.
All four gas frontends build, and the subtree check including actual predictor
endpoints passes (144.47 seconds;
`build/eb-held-godunov-parent-repair-20261006`). These checks establish native
replay; natural automatic root/nonroot selection remains unqualified.
Accepted-prefix replay now validates the complete real clock schedule, freezes
available accepted source tables and evaluates unvisited clocks for the first
time. Its second-parent group check passes native state/stage/face and ledger
comparisons, subsequent frozen-tail replay and owner/idle rollback (192.83
seconds; `build/eb-held-godunov-prefix-20261006`). Automatic parent trials seed
later weights from that first accepted complete table and include new clocks
in counts and all registration/source/pressure ledgers. A failure-only recovery
channel now returns complete accepted common source tables without publishing
failed physical state. Parent escalation retains original accepted packets and
the original incoming clock, adding only newly accepted tables from the failed
quarantined trial. Provider refusal or incomplete source stages cannot produce
a recovery table. The focused failed-closure/recovered-source check passes
first evaluation, frozen reuse and physical/history/ledger rollback (196.07
seconds; `build/eb-held-godunov-recovery-20261006`); automatic cache extension
and natural accepted-blend selection remain unqualified. Common native steps
now count complete accepted source stages, return a failure-only source prefix
and reuse just those packets on fresh physical caller state. Later stages run
natively for the first time. A failed new source/provider or capture cannot
publish a prefix that would authorize repeating that attempt. Full tables remain
required by existing full-clock replay. Four gas frontend/unit builds pass.
The focused source-stage check passes source0-only recovery after a finite
rho/N2 sink refuses the native predictor, new-stage evaluation on a fresh
admissible state, complete-table frozen reuse, exact serial/MPI fields/stages/
faces and idle count/NaN/provider-refusal rollback (37.55 seconds;
`build/eb-held-godunov-stages-20261006`). This is direct common-step recovery;
recursive stepping now accepts a last-clock source-stage prefix, restores
actual off-level old/last source rates before integration and connects common
predictor failures to native parent replay. Failed trial recovery appends only
new source stages while preserving original accepted packets and the incoming
clock. Missing/refused new provider attempts stop escalation. Four gas frontend
and focused unit builds pass. The directed child-prefix/parent check passes
four accepted provider stages without repeated calls during automatic parent/
root attempts, fresh seven-clock replay with only seventeen new source stages,
actual old off-level packet and independent integrated impulse, exact serial/MPI
fields/stages/faces, and malformed routed/count rollback (193.87 seconds;
`build/eb-held-godunov-stage-parent-20261006`). The deliberately inadmissible
sink refuses its original step; the fresh replay does not establish natural
automatic accepted-blend selection. Source cache metadata now records auxiliary
and provider acceptance independently from completed native source snapshots.
A last-clock frontier may retain either accepted component at stage zero;
packet flags agree collectively before selecting MPI callbacks. Native rates,
StateRedist and hydro still reevaluate on the fresh trial. Provider refusal,
malformed data or unsuccessful packet saving blocks failure-cache publication.
Quarantined recovery extends components without replacing original accepted
packets or incoming clocks, including trials stopping before an older accepted
provider is reached. Four gas frontend and focused unit builds pass. The
per-packet check passes count0 source recovery after finite native species-density
summation refuses, frozen packets on repeated native refusal, independent
auxiliary-only first-provider evaluation, child/parent source frontier and routed
old-rate retention, exact serial/MPI native completion, and malformed/idle
packet-flag rollback (110.24 seconds; `build/eb-held-godunov-packets-20261006`).
An initial negative raw sink completed Sources and refused its predictor; that
failed fixture assumption was replaced by the actual source-stage refusal.
The finite summation probe deliberately exceeds the native gate; it does not
establish physical source admissibility or natural automatic blend selection.
The EB particle/film sampling library now accepts independent UEINT valid and
ancestor snapshots, constructs the matching scalar halos, and recovers particle
cell primitives or film samples from that scalar. It retains background Eden and
physical gas deposition separately. Owner-local shape/bounds/finiteness and
optional-presence checks precede collective sampling; cached intermediate stages
reject changed UEINT mode. Four gas frontend and sampling runtime builds pass.
The four-rank directed 900K-Eden/700K-UEINT check passes exact 700K particle/film
exchange and deposition controls, malformed scalar rollback and staged mode
rejection (0.52 seconds); the existing no-UEINT sampling/stage/film check also
passes under runtime/overflow trapping (1.57 seconds;
`build/eb-held-spray-sampling-20261006`). Regular AMR particle/film stage and
halo wrappers also forward caller-prepared UEINT, retaining its sampling mode
across cached stages. Owner/idle sampling and rollback checks pass
(`build/regular-held-spray-wrappers-20261006`). The regular upper MOL coordinator
now receives native actual-UEINT endpoint packets and full ancestor snapshots,
keeps physical Nv particle/film sources separate, and publishes gas, populations,
UEINT and numerical reset together. A missing actual-scalar input in the regular
MOL halo adapter is corrected for serial/MPI AMR hydro. A three-rank directed
root/refined/wall-film test passes independent 700K UEINT versus 1000K UEDEN
sampling and endpoint agreement, plus idle-prefix refusal. Coupled conservative
remap, owner rebalance, checkpoint and CFL routes retain the actual scalar;
a rotating topology check passes 700K EOS/CFL, both stage histories and actual
owner migration. Shared fixed/selected MOL application builds pass. A bounded
nonreacting fixed application with explicit retained-scalar policy, scheduled
coarsening and 2-to-3-rank restart matches all 64 gas/UEINT checkpoint rows and
12 parcels exactly (`build/regular-held-spray-upper-20261006`). This is focused
connection evidence. Native serial/MPI regular MOL now passes the actual old/new
UEINT into typed diffusion callbacks and retains it in transport/LES ancestor
snapshots. The all-level common-step scalar updates no longer skip the loop when
no target depth is selected. Molecular transport and Smagorinsky, WALE, Vreman
and Dynamic LES recover primitives from actual UEINT; Dynamic filtering carries
the same independently filtered scalar. Stored halos use the existing linear
admissibility contract. Common/subcycle hydro and molecular bounds, LES estimates
and shared gas/spray application routes now carry UEINT too. A three-rank
nonuniform N2 directed check verifies identical physical transport/LES rates,
faces and bounds for actual 700K-family UEINT versus distinct 1000K UEDEN,
native common/recursive callback chronology, scalar endpoints and late refusal
rollback (`build/regular-held-diffusion-20261006`). Six fixed/selected gas/spray
frontends build; the original no-UEINT population control also passes. These are
focused implementation checks. Standalone MOL now supplies the actual old/new
UEINT to hydro and molecular-transport halos. Local and MPI periodic, moving-wall,
Mach-inlet and pressure-outlet constructors carry that scalar; Mach/pressure donor
EOS and application Mach timestep bounds use it too. A three-rank directed check
passes independently retained-scalar boundary states, native MOL callback stages
and late rollback (`build/regular-held-boundaries-20261006`); eight fixed/selected
standalone/tree gas/spray frontends build. This does not qualify whole-time Mach
applications or selected-mechanism runtime behavior. Cartesian MPI Dynamic LES
also supplies the same SI coefficient floors as serial; its strict boundary
RHS/SSPRK2 equality and rollback control passes. The rotating Cartesian coupled
reacting MOL connection now freezes evaluated sources on interval failure and
conservatively blends the new spatial family toward an admissible old-time family.
A shared MPI coefficient applies to gas, actual UEINT faces/rates and boundary
ledgers. Directed checks pass dry-tracer recovery, fractional trace-water
conservation and coefficient-failure rollback. The rotating coupled aggregate
previously reached its 180-second limit. Exact native-reactor input/timing
observation now identifies two failed LF/Euler integrations per rank, each
consuming the 100000-substep limit near the signed-trial bound. The native forced
interval now rejects a proven negative conserved molar pool before adaptive
stepping. The certificate uses connected integral-stoichiometry components whose
net mole count is conserved; inconclusive cancellation and nonconserved components
retain ordinary integration, and affine-negative species can still be replenished
by chemistry. The aggregate passes in 2.71 seconds without changing its timeout,
field/budget gates or spatial-source replay. A 3156-step replenishing-reactor
case has byte-identical baseline/new state and temperature
(`build/mol-forced-profile-20261007`). Full AMR recovery and physical/pinned
trajectory qualification remain separate.
Callback-presence
disagreement now refuses collectively before recovery reductions. This is focused
connection evidence (`build/regular-mol-window-recovery-20261007`), not qualification
of a complete AMR recovery envelope or pinned application histories.
Regular PPM/CTU now supplies actual UEINT to cell EOS and exports the final
Riemann scalar flux and `pdivu`, including artificial viscosity and the inflow
PLM face replacement. The shared AMR Mach-inlet donor now supplies actual UEINT
to its pressure EOS. Serial/MPI stage-only tree checks pass nonuniform and signed
formation energy, prescribed boundaries, parent endpoint interpolation and
collective rollback (`build/regular-godunov-held-spatial-20261007`). Native
regular AMR chemistry now forwards a separate old UEINT input to the reactor;
unpacking retains the pinned total/momentum energy projection. Common/recursive
serial/MPI inert-reactor checks pass gas/T/UEINT/IR/heat equality and late refusal
rollback (`build/regular-held-old-reaction-20261007`). Native nonreacting Godunov
common/recursive intervals now carry the full scalar flux family and final reset
correction; the MPI synchronization carries gas and UEINT together. Focused
root/refined/idle and signed-energy checks pass
(`build/regular-godunov-held-native-20261007`). Native regular intervals now export
signed pressure work over accepted composite leaf clocks. The independent UEINT
inventory, serial/MPI equality, presence refusal and late reaction rollback checks
pass; MOL exports zero. Optional V2 retained-energy checkpoints carry this signed
history beside the numerical correction. Serial/MPI, changed-rank, indexed and
late-corruption checks pass (`build/regular-godunov-pressure-work-20261007`). Fixed
and selected serial/MPI gas tree frontends now accept Godunov retained-energy
controls, preserve the scalar through remapping and bind its cumulative pressure
work to restart. Four bounded nonreacting N2 conditions pass common/recursive
three-interval execution, exact serial/MPI final CSV and 2-to-3-rank continuation
of complete state/history. Nonfinite pressure history and changed-policy restarts
refuse (`build/regular-godunov-gas-app-20261007`). These checks do not establish
reacting application or pinned trajectory parity. Regular MPI AMR spray now
retains Godunov/MOL UEINT, uses parent old/new scalar samples in common and
recursive particle/film clocks, and commits pressure work with outer success.
Coupled retained-energy V2 checkpoints restore this history with rank-neutral
state and reject late film-capacity failure without publishing it. Fixed/selected
three-interval nonreacting common/recursive Godunov and fixed common MOL
2-to-3-rank restarts match gas/T/UEINT, 12 parcels/component masses and signed
histories exactly; nonfinite history and changed-policy restarts refuse
(`build/regular-godunov-spray-history-20261007`). These are bounded connection
checks with evaporation disabled. Inspection of the byte-verified pinned
`Source/Setup.cpp` now establishes that rotation requires MOL; the previously
listed rotating Godunov assembly would be a separate extension. Its completion
scope awaits the user's decision. Broader reacting/transport/LES qualification
and whole-time/pinned numerical and particle comparisons remain open.
The EB spray coordinator now unpacks
actual UEINT and ancestor snapshots, declares zero direct UEINT particle/film
source, retains the physical Nv gas correction, and commits the scalar,
numerical reset and pressure-work histories with outer gas/particle/film success.
Recursive native gas stepping dispatches the particle predictor while keeping
UEINT prediction separate; private post-level populations follow paired
synchronization. A directed 700K-UEINT/1000K-Eden four-rank fixture passes common
and recursive Godunov/MOL, early prefix disagreement and late split-capacity
rollback (31.69 seconds; maximum composite residual `2.15e-14`). It exposed
reset forcing charged from a provisional SDC iterate. Serial/MPI ledgers now
record only the reset actually injected into the accepted iterate; a three-rank
cut-tree check passes 1/2 iterations of Godunov and MOL, matching actual UEINT,
state and histories (12.72 seconds; maximum residual `1.05e-14`). Four gas
frontends, fixed MPI EB spray and fixed/selected spray runtimes build. The spray
application and coupled checkpoint carry UEINT plus separate numerical and
pressure-work histories; a bounded nonreacting common-Godunov 2-to-4-rank restart
matches 216 gas/UEINT rows, 2 particles and 12 budget components exactly
(`build/eb-held-spray-coordinator-20261006`). That application check precedes the
accepted-iterate ledger repair. Paired MPI interface repair now carries actual
UEINT through the same complete gas donor and face packets; its EOS queries use
the independent scalar. A two-owner/one-idle directed check passes signed scalar
conservation and atomic EOS refusal. The previously refused nonreacting recursive
application continues from its accepted step-3 checkpoint to the unchanged
`2e-8 s` end time, with reported inventory error `1.95e-18`
(`build/paired-held-interface-budget-20261006` and
`build/eb-held-spray-recursive-app-20261006/paired-budget-fixed-restart`). This is
checkpoint continuation evidence, not a fresh whole-time or reference comparison.
The existing gas-only donor control now passes in 21.85 seconds within its
unchanged 60-second limit. Serial/MPI alternate-route search stops when the
complete packet cannot replenish any actual density/species deficit, avoiding
exhaustive face exclusions in an unfundable covered endpoint. The original
covered/corridor success, refusal and rollback fixtures remain unchanged.
The lower reflux and signed held-packet checks also pass (0.62 and 0.47 seconds);
source/binary identity and logs are in `build/eb-covered-route-pruning-20261007`.
An earlier test started before build completion was interrupted and is not
counted. Bounds-check diagnostics also exposed a four-node test guard referencing
`ownership(5)`, now corrected with explicit nested branches. Serial paired repair now matches MPI admission, gas,
temperature, actual UEINT and both face families exactly. A directed registered
impulse check also verifies that gas and signed scalar debt are applied once.
Prepared serial/MPI StateRedist closure forwards both impulse families and their
ancestor spill range to paired repair. A nonreacting pure-N2 cut-tree control
passes Godunov/MOL with 1/2 iterations (39.75 seconds), matching state, actual
UEINT, reset and pressure-work histories. Paired repair now funds exposed
negative ancestor endpoints outside the adjacent-level face loop, using the
same complete gas/UEINT donor packets and the donor level's actual route extent.
A three-level directed register fixture matches serial/MPI admission, gas,
temperature, actual UEINT and both face families exactly for distributed and
owner-zero/idle layouts. Post-debt actual-UEINT EOS, signed composite inventory,
read-only debt, insufficient-inventory refusal and EOS-refusal rollback pass.
The existing lower gas, three-rank covered-route and expanded held tests pass
(0.67, 21.24 and 0.84 seconds); identities and the pre-fix rejection are in
`build/eb-held-ancestor-budget-20261007`. This qualifies directed exposed-ancestor
funding. Gas-only registration closure now also forwards the ancestor spill
range into serial/MPI repair, with collective range/presence agreement. A cut-cell
register fixture creates negative water debt on an exposed ancestor through the
actual re-reflux/spill assembler. Its two/three-rank distributed, reversed-owner
and owner-zero cases match serial gas/T/faces exactly, conserve every composite
component and preserve the unrelated older clock ledger. The same fixture rejects
with independently rebuilt pre-change modules. Directed gas debt/EOS refusals
and mixed/invalid MPI range controls also preserve all inputs.
The older positive register-close control's physical-boundary loss is now
resolved. A compiled probe using the pinned AMReX EB2 geometry and re-reflux
kernel retains `0.83367781155015186` of a source next to one nonperiodic domain
face and `0.72279635258358643` at a two-face corner; the interior retains exactly
one. The kernel distributes through geometric connections into exterior cells
that are absent from the accepted state. Production serial/MPI register closure
now supplies the topology's physical domain and periodicity. Nonperiodic exterior
recipients are excluded before the unchanged kernel normalizes its weights;
periodic images and peer/ancestor halo recipients remain available. The original
positive and ancestor register closures pass on two/three ranks with exact
serial/MPI gas/T/faces and retained older ledgers. The signed three-component
low operator checks full boundary/corner inventory, exact interior and periodic
image output, and invalid-domain output rollback. The native pure-N2 held
Godunov/MOL register/reset control also passes (39.48 seconds). The boundary state
update is a conservation correction relative to pinned AMReX, not pinned boundary
trajectory parity. Source byte checks, independent compiled-reference measurements and
before/after failures are retained in `build/eb-rereflux-physical-boundary-20261007`
and `build/gas-ancestor-register-range-20261007`.
Older exposed ancestors now also accept covered donors on their own patch;
adjacent-level interfaces retain their existing fine-donor selection. Donor
inventory comes from physical leaves after pending gas/UEINT debt, including
private fine debt once. Each descendant is qualified before spending its mean;
pre/post-register bounds and actual-UEINT EOS constrain the routed packet before
trial state and faces commit. Two/three-rank directed gas and held fixtures cover
children and grandchildren, complete signed packets through a covered transit,
valid descendant debt, insufficient inventory, a hidden bad leaf EOS and a
negative leaf budget, with exact serial/MPI fields/faces and unchanged live debt.
The original gas donor control also passes (22.73 seconds). Identities and the
pre-change rejection are in `build/covered-ancestor-donor-20261007`.
The donor fixtures now also cover extruded and nodal sloped cut geometry,
with physical-volume-weighted source inventories, child/grandchild providers
and cut covered transit. Gas register closures pass on two/three ranks and
paired held budgets on three ranks; exact serial/MPI states/faces, once-only
leaf debt and rejection rollback remain checked at their existing tolerances.
Evidence is in `build/cut-covered-ancestor-20261007`. A registered endpoint
with uniform pre-state can now borrow a complete physical donor packet when
mean contrast cannot improve its deficit. Gas and signed UEINT travel together;
all gas deficits must close with donor inventory retained. The borrowed density
covers the existing EOS density floor plus represented roundoff, without
clipping state or weakening EOS. Generic flux-only contrast pooling stays
unchanged. Low-level zero/negative-density and missing/insufficient-inventory
controls, two/three-rank ancestor closures, cut held budgets and the native
pure-N2 Godunov/MOL register/reset control pass. The pre-change refusal and
initial below-floor EOS failure are retained in
`build/uniform-register-inventory-20261007`. Registered repair also scans
exposed cells on the closing level away from adjacent fine interfaces, followed
by every registered older ancestor. Gas/held far-endpoint fixtures cover both
owner/idle layouts and once-only full-component source inventory; two/three-rank
closures, the original gas donor control and native pure-N2 register/reset control
pass. Flux-only closure retains its original scan range. The pre-change rejection
and exact identities are in `build/all-exposed-register-endpoints-20261007`.
Exposed registered actual-UEINT endpoints now fund thermal-only deficits using
signed NASA caloric bounds at the mechanism's common fitted temperature limits.
One fraction moves every gas component and UEINT; healthy pre-register targets
can pool several physical donors, while complete post-debt EOS is required before
publication. Donor pre/post EOS and gas bounds remain strict. Nearly exhausted
stock is scaled as a whole material vector, preserving rho=sum(rhoY); target and
donor debts are reconstructed in the same order as the callers. Native prepared
closure retries paired repair after a private reset/EOS refusal and charges only
the accepted reset. Existing covered-leaf EOS/inventory refusals remain intact.
Cold, hot, signed water/negative-UEINT, multi-donor and impossible thermal fixtures
pass with exact serial/MPI fields/faces, original conservation tolerances and
once-only live debt. Six focused controls pass, followed by the same five thermal
cases using actual `h2o2_full` bounds (300--3500 K). Evidence and all intermediate
failures are in `build/held-eos-register-funding-20261007`. Its read-only native
observer measures a coarse signed UEINT register correction of `-13220.3972`,
but the post-register state is healthy and no repair call occurs. This is a
nonzero-register control, not a native forced-deficit qualification. Native forced
SDC deficits, fine-leaf EOS-only repair, broader descendant thermal
stock distribution, whole-time recovery and pinned trajectories remain open.
Exposed registered gas-only endpoints now also borrow whole physical material
packets to fund thermal or kinetic EOS deficits. Signed NASA caloric bounds
include momentum and intersect every pre/post target and donor constraint;
disconnected admissible fraction intervals are retained for counterstreaming.
Several donors may pool a thermal deficit while their pre/post EOS remains
healthy; complete post-debt EOS is required before publication. Two/three-rank
`h2o2_full` fixtures cover cold/hot bounds, water with negative formation energy,
multiple providers, counterstreaming and insufficient stock with atomic refusal.
Serial/MPI fields/faces match exactly, full-component source inventory and live
debt remain checked at their original tolerances. Both new tests and seven
existing reflux/ancestor/held/native-register controls pass. The initial failing
gas fixture, corrected declaration failure and exact tested identities are in
`build/gas-eos-register-funding-20261007`. These gas funding fixtures are
nonrotating; combined gas-negative/thermal partial pooling and native
forced-deficit qualification remain open.
Registered exposed fine leaves now fund thermal-only gas/actual-UEINT deficits
before final EOS publication. Their private state already includes the register;
funding reconstructs its healthy pre-state and qualifies donor pre/post stock,
then retains post-state privately until the original read-only debt is removed
for strict closure. Complete gas/UEINT packets travel along straight same-level
paths, including covered transit and terminal providers. Covered stock comes
from qualified physical post-register leaves; routed allocation checks every
changed descendant's actual pre/post EOS and rolls back a refused route before
trying another donor. Exposed transit cells retain zero net divergence.
Cold/hot, signed water, multi-provider and insufficient-stock fixtures pass for
regular and cut cells, including healthy donor debt, exact serial/MPI fields/
faces, once-only inventory and atomic refusal at the existing tolerances.
Covered fine fixtures also force remote descendant stock, two distinct covered
providers through a covered transit, and hidden bad-leaf EOS refusal. Gas two/
three-rank and held three-rank checks pass. Evidence is in
`build/fine-register-eos-funding-20261007` and
`build/fine-covered-eos-funding-20261007`. Cross-owner same-level peer/bent-path
qualification, broader nonuniform descendant thermal allocation, coupled coarse/
fine deficient-stock ordering and native forced SDC deficits remain open.
Exposed raw fitted-range EOS deficits now also borrow complete material packets
with absent or explicitly zero registers. Density/species and composition must
already be physical; zero-debt opt-in permits private thermal pooling while
donor pre/post EOS and complete final target EOS remain strict. Registered
targets retain the existing pre/post admission. Closing, fine and covered fine
regular/cut fixtures pass cold/hot, signed water, multi-provider and atomic
refusal checks with exact serial/MPI state/faces and unchanged inventory
tolerances (`build/unregistered-eos-funding-20261007`). A native Godunov SDC
one-cell UEINT sink of `-3e13` at `dt=1e-8` remains inside the fitted caloric
range. Increasing it to `-3.1e13` produces a predictor EOS refusal before
prepared closure, on both the old diagnostic binary under controlled forcing
and the current implementation; prefix/work outputs remain uncommitted.
The pinned Fuego/LiDryer source audit then identifies a different reference
contract: `GET_T_GIVEN_EY` linearly extrapolates outside 300--3500 K, and Fuego
accepts the returned temperature despite `ierr=1`. The exported headers match
the pinned upstream bytes. A compiled generated-header probe returns 290 K,
3525 K and 1000 K for below-range, above-range and interior energies respectively
(`build/unregistered-eos-funding-20261007/reference-eos`). This is an isolated
EOS branch check, not full PeleC/CFD parity. Explicit
`pinned_eos_extrapolation` now binds per-run EOS tables to regular 3D,
serial/MPI two-level and regular AMR-tree, native EB-tree and regular single-EB
paths, including selected variants. The inverse uses the common fitted bounds
by default; `pinned_eos_temperature_min/max` accepts a positive ordered pair
for mechanism.H exceptions, without changing fit metadata or coefficients.
Outside those inverse boundaries it uses linear boundary-Cv extrapolation;
positive finite polynomial property evaluation shares the same table policy.
Its in-range inverse uses the pinned default 200 Newton iterations, +/-100 K
updates and a 1e-6 K step criterion. The strict library inverse retains its
existing residual criterion. Mechanism sources and SHA remain unchanged.
LES callbacks share the per-run table. Restart identities include effective
EOS mode and bounds; strict checkpoint bytes remain unchanged in seven
completed regular/tree baseline comparisons. Unsupported frontends reject
the input instead of ignoring it. Existing MPI tree transport/full/timestep
parameter packets now reject differing EOS policies before physics callbacks
without adding collectives; the two-level application also checks agreement
before initialization. These checks do not yet audit every lower-level MPI API.

Nine focused EOS/checkpoint/native-ledger tests pass after the wall fix.
The native cooling fixture accepts positive sub-300 K UEINT and its independent
source/face/pressure-work balance has residual `3.24e-16`. Its complete gas
ledger exposed a wall-pressure mismatch: Godunov flux used total-energy pressure
while its wall ledger used UEINT pressure. Serial/MPI wall accounting now follows
the gas provider, consistent with pinned `pc_ctoprim`; exact serial/MPI states
and ledgers pass. Three MPI controls rerun after policy-packet changes pass,
including differing mode/bounds rejection with untouched state, ledgers and
zero source evaluations. Eight regular/tree CLI cases and two invalid-input
controls pass. The two-level connection builds in all four fixed/selected,
serial/MPI targets: 18 short successful CLI runs and eight checkpoint-policy
refusals pass, including 3-to-2-rank and serial-to-MPI restart. Eight corresponding
serial/MPI CSV comparisons pass at the existing `2e-12` gate. These are bounded
application checks, not full-window trajectory qualification.

Isolated compiled LiDryer/N2 inversion agrees at 290, 3525 and 1000 K with
identical header coefficients. The modern heptane_lu_qss header uses 300--4000 K
rather than its YAML table's 300--3000 K common fit interval: an explicit override
agrees at 290, 4025 and 1200 K with maximum relative difference `1.13e-15`.
The complete exported HP_DME and heptane_lu_qss headers match pinned upstream
bytes. The legacy HP_DME header instead rounds polynomial coefficients after
integration; its 4025 K comparison still differs by `7.90e-6` K (relative
`1.96e-9`) and fails the unchanged `2e-12` gate. That legacy-kernel parity gap
remains open. Probe evidence does not establish mechanism import, QSS chemistry
or whole CFD parity. Logs, source identities, earlier failures and restart
instructions are in `build/pinned-eos-extrapolation-20261007`.
Spray regular serial, regular MPI, MPI regular AMR and single-EB owning
libraries now configure a private EOS table before property evaluation;
fixed/selected frontends, film/LES/reaction callbacks and restart use that table.
MPI input/normalized-policy agreement covers the new controls. Forty bounded
runs/restarts and 16 isolated mode/bounds refusals pass across eight frontends,
including cold/hot gas, particle heating and 2-to-4-rank AMR continuation.
Strict regular serial/MPI checkpoints match saved binaries byte for byte;
fixed/selected AMR component bytes match after excluding only each transaction's
fresh generation and its exactly verified context-length dependency. Source,
commands, initial input/checker failures and artifacts are retained in
`build/pinned-eos-spray-bindings-20261007`. Four fixed-frontend EOS application
regressions are registered. These checks do not establish reference particle
trajectories or physical validation.
The namelist-driven 1D and 2D regular/AMR/EB owning simulations and output
paths now use private EOS tables as well. Strict-mode 1D outputs match saved
prebinding binaries for ten regular/two-level/multilevel/multipatch/sparse-MPI
cases; ten cold and ten hot cases, rank-changing sparse restart and policy
refusals pass. The 2D two-level, multipatch, static/dynamic three-level and tree
checkpoint families retain strict header bytes and bind active EOS identity.
The shared MPI 1D policy check covers regular, replicated-tree and sparse-tree
public physics/regrid entries; registered two/three-rank fixtures reject differing
mode, bounds and mixed per-species policies before modifying states. Selected
1D/2D checkpoint baseline energy may be negative under the NASA reference:
finite values, density/species conservation and the actual cell EOS remain the
admission checks. Strict fitted-temperature water with negative energy is accepted,
while a finite EOS-invalid cell is refused with existing outputs retained.
Artifacts and exact bounded source/binary identities are in
`build/pinned-eos-1d-bindings-20261007` and `build/pinned-eos-2d-bindings-20261007`.
These are connection/restart checks, not pinned full-trajectory or physical
qualification; low-level 2D MPI refusal qualification remains separate from its
completed application-level checks.
The gas, actual-UEINT and biased-reactor donor search now uses the active EOS
inverse limits rather than treating NASA fit limits as admissibility bounds.
Pinned extrapolation searches above its zero-temperature energy intercept
without an artificial fitted ceiling; every proposal still passes finite-property,
positive-temperature and Cp/Cv EOS checks. Eighteen cold/hot and explicit-bound
material-packet fixtures plus six strict/invalid-Cv refusal fixtures pass in the
existing reflux unit. General narrow Cp/Cv admissibility pockets still need
further search qualification; selected/prepared/deeper clock funding must also
retain offlevel gas and independent-UEINT pending impulses before its gates can
be extended. Bounded identities and logs are retained in
`build/native-extrapolated-donor-20261007`.
Common two-level nonreacting Godunov/MOL native predictors and correctors now
fund represented material-source shortfalls with complete gas/actual-UEINT donor
packets and matching conservative faces. Corrector admission uses the actual mean
source and spatial face family, including diffusion; a new-time source can exhaust
stock even when its old-time predictor was admissible. Virtual covered coarse operator values are excluded from physical source
debt and their stages are reconstructed from funded fine states before subsequent
source evaluation. Nonreacting three-rank fixtures (two owners and one idle rank)
retain the original complete material source, close gas and independent UEINT
balances, and match serial states exactly. Godunov frozen-source replay also
reproduces gas/actual-H/temperature, gas and scalar faces, energy prefixes and
source ledgers exactly without callbacks in the material-source regressions.
Corrector funding refreshes accepted hydro/face records while retaining the
original native Godunov predictor input. Three staged-source regressions pass,
including covered coarse forcing and diffusion. Globally unfundable predictor
sources retain public states, ledgers and energy prefixes. Evidence and
rejected/invalid diagnostic attempts are retained in
`build/native-physical-recovery-20261007` and
`build/native-corrector-recovery-20261007`. This is bounded native-step evidence;
common reacting held-source recovery now also checks actual old-reactor energy
with a read-only caloric offset, while transporting gas and actual UEINT together.
Inert chemistry-path regressions pass for physical and covered corrector sources,
including serial/MPI ownership and exact Godunov replay. Covered reacting forcing
is reconstructed from admitted physical descendants before reactor packing;
physical leaf sources and the actual old-UEINT clock remain unchanged. A supplied
reactor receives the ordinary native attempt without a forced-composition caloric
preview preempting it. Four bounded Godunov/MOL refusal probes retain public
states and ledgers. These inert and manufactured-provider checks do not establish
pinned-trajectory parity. The manufactured caloric fixture passes, including a
reaction that replenishes signed forced water and accepts a rejected caloric
preview without material loans. Actual H2/O2/N2 kinetics now passes the bounded
physical, covered and diffusive corrector-source fixtures and the built-in
no-callback admission path: nonzero chemical change, mass/element closure,
serial/MPI states and exact Godunov replay are checked. This does not establish
physical validation or pinned numerical trajectories. Artifacts are in
`build/native-reacting-source-recovery-20261007`. Reacting diffusion, selected
clocks, prepared registers and deeper hierarchies still need source-recovery
integration/qualification. Remaining
other-dimensional frontend EOS binding, application qualification of newly
EOS-bound reference preparations, legacy HP_DME kernel parity and
full-window/pinned trajectories remain unfinished.
The forward limiter now applies its existing safety factor to a saturated positive
outgoing budget too, covering a smaller face lost in the summed allowance's
rounding; the extracted regression and shared-face rollback check pass. Native
and CVODE forced reactors now
evaluate signed Newton trial concentrations with integer reaction orders and
their analytic concentration derivatives. Accepted native interval endpoints
retain strict physical admission; Richardson extrapolation keeps its physical
check and falls back to the half-step result when it reverses a trace product's
sign. The actual negative-forcing cell, the existing N2/water trace tests and
the no-UEINT two-rank reacting SDC fixture pass, as do two clean CVODE controls
(`build/forced-signed-trials-20261006`). These checks do not establish negative
third-body/falloff or broader mechanism qualification. No source clipping or
tolerance change was used;
reacting/StateRedist/rotating/selected application coverage and fixed-reference
particle trajectories remain unfinished.
Focused patch checks pass
regular and three-axis cuts at weights zero, 0.375 and one, retaining exact
native pressure/wall volume work (18.78 seconds). A three-rank selected cut-tree
check passes raw/StateRedist, two owner layouts with an idle rank, exact
serial/MPI fields, halos, faces and registers, independent UEINT inventory
(`1.73e-18`), and complete rollback for idle weight/presence disagreement, NaN
and bounds refusal (67.56 seconds). Evidence is retained in
`build/eb-held-godunov-replay-flux-20261006`.
Valid-only held
StateRedist, routed/nonzero-UEINT auxiliary replay,
remaining Godunov and particle/film coupling, natural prepared refusal
qualification and pinned trajectory comparison remain required.
The selected-mechanism rotating recursive EB fixture previously rejected its
first step at `dt=1e-8`: chemistry consumed provisional coarse H2O2 inflow before
reflux removed it, leaving a negative composite inventory that conservative
redistribution could not fund. Coarse-only replays were rejected and removed.
MOL exports its accepted local forcing, density/species volume source and
physical old spatial flux. The recursive driver retains sparse owner/level
records and can reconstruct reaction forcing from changed face fluxes.
Accepted records identify replayed solves explicitly; their initial state,
forcing and endpoint describe the replacement solve, while the captured old
spatial flux remains the original interval-group anchor.
For reacting recursive MPI EB MOL with no redistribution, an exhausted
interface repair now privately replays all
coarse/fine subintervals in the failed parent interval. Their time-consistent
old spatial flux is the low-order anchor; one blend for all conserved face
components retains as much of the provisional MOL flux as the unchanged
reaction/EOS gates admit. Frozen wall/volume/particle source families are kept,
while gas, temperature, reaction/heat history, chemical inventory and registered
flux are committed together. Routed auxiliary impulses are retained separately:
accepted descendant closures are replayed in order, while the failed parent's
unapplied bucket remains for its original post-sync closure exactly once.
The serial and MPI routed source weights now match the selected-level MOL
forcing: the first reacting pass retains both endpoints; subsequent reacting
passes retain only the old half. A focused test rejects the preceding weights.
The reacting MOL spatial-flux budget safeguard now also covers the first pass:
reconstructed new-time transport can otherwise withdraw an absent species before
any level closure. Its first-pass budget includes both accepted volume-source
endpoints; later passes use the surviving old family. If the preceding chemical
prediction would overdraw the old-state low-order anchor, that anchor uses only
the accepted external family and chemistry is recomputed after flux limiting.
Admissible passes keep their original fluxes. This additional conservative
safeguard differs from upstream negative-species cleanup; pinned trajectories
remain unverified.
The selected rotating evaporating two-step fixture completed at `dt=1e-8`
before the first-pass budget extension. That run's scheduled regrid/injection and
two-to-four-rank restart passed with exact gas,
temperature and reaction/heat history; the final two-step composite inventory error is
`1.38e-17`. Normal local and recursive history checks retain the original results
and source callback counts, including idle-rank and late failure rollback.
Prepared StateRedist retains its prior synchronization. Routed impulse archive and
normal three-level source accounting have focused checks. A synthetic two-pass
predictor case also completes without replay or extra native reacting solves.
A synthetic
three-level one-pass case now exhausts parent repair and replays all seven
native reacting intervals, with frozen routed sources, the middle-closure
impulse checked before the next middle interval, exact owner-layout gas/T/IR/heat
and final-callback rollback. Redistribution interval archives now retain old
and accepted spatial/volume source families, full owner-local transfer impulses,
and prepared StateRedist register masks, scales and deltas. Prepared MOL
archives also retain every old/new source stage separately, including six-layer
state/raw-RHS halos, unweighted faces and registers before reacting-stage
coarse-register accumulation. Focused stage re-evaluation reproduces the
native RHS and registers exactly. A reconstruction API now rebuilds the accepted
mean forcing and all intermediate coarse-register contributions from supplied
stage inputs, retaining the accepted forcing exactly for an unchanged spatial
anchor. A focused changed-stage case agrees with native re-evaluation;
an idle-rank control refusal preserves output forcing and registers.
Prepared intervals also retain the initial reaction predictor and each stage's
external/auxiliary volume sources. The native `replay_interval` path rebuilds
fresh MOL halos, spatial operators, registers and chemistry, without calling
frozen auxiliary providers again. Focused identity and owner-local reactor
rejection preserve gas/T/IR, forcing and register/ledger rollback. A changed
momentum state is retained in the rebuilt initial halo and produces a new
native raw RHS and reacting endpoint on distributed and idle-rank layouts.
Native prepared replay now accepts one rank-neutral spatial weight for every
MOL stage. It mixes first-order old-state LF divergence and native non-face
sources with the freshly evaluated StateRedist RHS, and applies the same weight
to all face components and coarse/fine register deltas. Saved stages retain
both families; focused weights 0, 0.5 and 1 reproduce their full forcing and
registers, while weight 1 retains the original gas/T/IR exactly. Interval
archives also retain the actual wall ledger for replacement after replay.
The native prepared StateRedist parent-replay implementation now restores all
owner-level clock states and the still-open parent register before rebuilding
fresh old/new stages, chemistry and descendant closures. Wall/chemical ledgers,
gas/T/IR/heat, records and register families publish together. The first
three-level fixture rejected a local withdrawal of absent peroxide before
reaching synchronization. A local complete-family spatial budget now addresses
that path, retaining accepted volume forcing and limiting all native stages,
faces and register deltas with one weight. The zero-peroxide fixture now passes
all seven native intervals with no extra reacting solves, independent middle
register/face closure, exact gas/T/IR/heat across owner and idle-rank layouts,
source/post schedules and root/middle callback rollback. Its composite inventory
error is `4.15e-16`. This fixture exercises the local budget and normal closures.
Local family selection now also checks native chemistry/EOS admission when the
conservative spatial budget alone leaves insufficient reactant. Rejected trials
remain private; the accepted gas/T/IR/heat and every spatial stage, face and
register use the same selected weight. A supplied reactor refusal on a
native-admissible cell aborts the transaction instead of changing that weight.
The same supplied-reactor refusal rule now covers ordinary MOL and FluxRedist
parent trials: a refusal on an independently native-admissible input aborts
without another weight trial, preserving the caller's gas/T/IR/heat, archives
and ledgers. Focused distributed and owner0/idle fixtures exercise that gate.
Selected local spatial means now remain authoritative through chemistry,
face accounting, interval capture and off-level transfer publication; encoding
an endpoint and averaging it back had exhausted a remaining trace-water budget
through cancellation. FluxRedist checks the actual composed forcing before
accepting a single weight for the complete conservative family.
The longer-clock trace fixture accepts all seven local intervals, but does not
trigger parent replay and therefore fails its replay-coverage assertion;
whole-clock replay qualification remains pending. Unprepared StateRedist replay
is still open. A focused archive
check exercises nonzero cross-level corrections for both methods, with exact
gas/T/IR and exported source/transfer/register families against calls without
capture, including late source-failure rollback on distributed and idle-rank
layouts. FluxRedist joint replay now applies its complete spatial operator to
blended canonical faces and frozen accepted non-face sources, retaining closure
spills for replacement inside open ancestor clock windows. Its first/later local budget
path also blends these families together rather than retaining an unfunded
face-independent redistribution residual. Its low family uses first-order
old inviscid Lax–Friedrichs MOL hydro with frozen accepted non-face diffusion/wall sources
with no redistribution in the low family. The accepted high family retains
FluxRedist; one weight limits both active-cell withdrawal and off-level
withdrawal against actual old species inventory and accepted volume sources.
Reaction prediction is not treated as available transported mass. The new
three-level native fixture now passes with root-clock replay after a descendant
closure cannot be funded locally, exact gas/T/IR/heat across distributed and
idle-rank layouts, composite conservation and final-callback rollback.
Only completed end-of-clock physics closures may escalate; provider and
post-callback rejection still aborts the outer transaction. Deferred descendant
post callbacks run once after accepted ancestor replay.
StateRedist replay,
strongly coupled three-level qualification, temporal accuracy, full-duration
runs and pinned trajectories
remain open. These bounded runs do not establish full physical qualification.
The focused nonrotating coupled library test passes, including films, reacting
callbacks and collective rollback. Fixed common and selected recursive LF/LES
application restart pass for the inert split-particle fixture, with 2-to-4-rank
continuation, regrid and MOL-settings rejection. These runs disable drag, heat
and evaporation and do not qualify evaporating spray. The fixed common and
selected recursive LF/LES evaporating restart fixtures now pass after a
conservative coarse/fine budget correction. An inadmissible coarse correction
returns one conserved vector to a connected finer leaf region; all conserved
components use the same weights. The boundary supply and internal transfers
also update every corresponding face copy, including periodic seams and
covered ancestor faces. Separate internal energy can travel in that same
vector. The two-/three-rank mixed-level, sibling-patch and periodic fixtures
check independent face divergence, composite inventory, exhausted boundary
donors and atomic insufficient-budget rejection. This conservation policy
extends the upstream species cleanup; pinned-reference equivalence remains
unverified. Strict accepted-species checks, tolerances and golden hashes are
preserved. Broader EB positivity and rotating-particle conservation qualification
remain open.

MOL-specific characteristic slopes, primitive face reconstruction and the
mixture Riemann interface flux are now implemented with focused analytic and
failure checks ([reconstruction](validation/mol-reconstruction-20260927.md),
[interface flux](validation/mol-riemann-20260927.md)).
The interface also implements the pinned Lax–Friedrichs diffusion speed and
flux formulas, with analytic checks. `mol_lax_friedrichs` now selects this
branch throughout the MPI tree application and binds it to restart identity.
Gas and single-level particle checks pass; a nested particle case produces
negative fuel density with the original centroid flux
([LF kernel](validation/mol-lax-friedrichs-20260927.md),
[application evidence](validation/mol-lf-application-20260927.md)).
Prepared EB patches now assemble centroid-corrected MOL spatial rates and
boundary/wall ledgers with stationary balance and conservation checks.
The application connections are described below
([spatial rate](validation/eb-mol-spatial-rate-20260927.md)).
MOL now uses its pinned centroid-stencil rules separately from CTU; focused
analytic checks and acoustic-MOL nested restart pass. The LF negative-fuel
case is traced to neighboring diffusive fluxes entering centroid interpolation
([stencil and diagnosis](validation/mol-centroid-stencil-20260927.md)).
An additional conservative limiter now bounds LF centroid corrections on
interior faces by available species mass. MPI MOL also exchanges the donor
allowances at same-level peer faces before publishing boundary fluxes and
rates; a cut-peer fixture checks shared flux, nonnegative forward fuel and
species conservation. Serial tree MOL now uses the same allowance at peer
faces; the active cut-peer fixture matches MPI rates and fluxes in the
registered two-rank MOL test. Fine patch boundary corrections also consult
an available coarse donor allowance in common-level spatial evaluation.
A periodic cut-cell fixture matches serial/MPI rates and fluxes and retains
nonnegative fuel after peer and seam synchronization. In the stressed full
coarse/fine synchronization, level-1 node 3 cell `(1,3,1)` starts with
`5.15e-5` H2 density and receives a `-6.18e-5` reflux correction. The
strict transaction rejects that uncorrected state. The existing MPI interface
budget repair now passes a test that requires actual repair activation and
accepted nonnegative species in this stressed nested periodic fixture. Both
two- and three-rank runs use one interface repair, retain the unrestricted
signed composite vector to `6.94e-18` relative error, and match the independent
open-face divergence to `1.64e-16` relative error. The test also checks finer
leaf inventories through the covered coarse region. These focused results
qualify this conservative repair path; broader cut-cell geometries, complete
MOL trajectories and fixed-reference equivalence remain unverified.
The previously failing two-step nested particle case
completes with a `5.42e-20` maximum budget residual; gas CSV and particle rows
match across 1-to-2-rank restart. Broader periodic and coarse/fine
positivity, recursive coarse allowance exchange, broader multilevel
qualification and fixed-reference numerical
comparison remain open.
Accepted MPI SDC states now reject even tiny negative species before
publication; a two-rank red/green rollback test and nested restart pass.
The upstream clipping cleanup changes species inventory and is not applied
([acceptance repair and conservation boundary](validation/sdc-negative-species-20260927.md)).
Serial SDC now applies the same strict accepted-species contract before output
publication; a negative-source rollback test and the existing serial tree test
pass ([serial acceptance repair](validation/serial-sdc-negative-species-20260927.md)).
Owner-local MPI MOL rates and fluxes now match serial evaluation in a nested
wall fixture, including empty-rank ownership and collective rollback
([MPI spatial evaluation](validation/mpi-eb-mol-rate-20260927.md)).
Nonreacting common MPI steps now combine old/new MOL rates, transport and wall
exchange with trapezoidal weights. Focused heat, viscosity and nonuniform
inviscid cases conserve inventory ([time update](validation/mpi-mol-time-20260927.md)).
Reacting common MOL now re-evaluates the preceding reacted state between
iterations; a homogeneous reactor comparison and late-reactor rollback pass.
([reacting update](validation/mpi-mol-reaction-20260927.md)).
Recursive MOL now forwards the same update through level subcycling; focused
nonreacting wall and state-redistribution budgets pass on two ranks.
Recursive gas MOL is connected
([recursive update](validation/mpi-mol-recursive-20260927.md)).
Recursive particle MOL now passes a short two-level case and exact rank-changing
restart after preserving exact face symmetry in invariant nodal geometry
([repair](validation/eb-face-symmetry-20260927.md)).
Common particle MOL now passes its short evaporating fixture and exact
rank-changing restart after adding a positive-face slope bound. This safeguard
extends the pinned reconstruction; full reference parity remains open
([repair](validation/mol-positive-faces-20260927.md)).
Serial tree common and recursive MOL evolution now also connects reaction,
transport, redistribution and restart. Four three-level cases match MPI state
and temperature exactly; serial-to-MPI continuation passes
([serial connection](validation/serial-mol-tree-20260927.md)).
The MPI gas application now exposes `mol_mode`, binds it to checkpoint identity,
and passes a short heated-wall case with exact rank-changing continuation
([application checks](validation/mpi-mol-application-20260927.md)).

The pinned Spray-EB plane and five initial physical droplets now have an
auditable SI input preparer preserving the default 64³ grid and final time.
An 8³ short case now completes with five parcels after fixing reflected-wall
source donor positions; full-duration/reference parity remains open
([case](../cases/spray_eb_reference_3d/README.md),
[repair and evidence](validation/spray-eb-source-halo-20260927.md)).
The preparer now selects MOL as in the pinned input. Its 8³, 200 ns two-rank
case passes with wall heat, no-slip traction and five particles
([MOL connection](validation/spray-eb-mol-connection-20260927.md)).

A 3D embedded-wall quadratic gradient and scalar heat-flux kernel now follows
the pinned type-0 stencil. Its coefficients agree with a compiled reference in
24 orientations ([kernel evidence](validation/eb-boundary-gradient-20260927.md)).
An optional SDC transport source now returns gas energy rates and integrated
wall heat, with analytic heat/ledger and rollback checks
([source integration](validation/eb-wall-heat-source-20260927.md)).
MPI two-layer donor exchange and composite wall-heat integration now pass a
cross-rank stencil case and exclude covered coarse cells
([MPI integration](validation/mpi-eb-wall-heat-20260927.md)).
Common and recursive SDC evolution now include wall heat in their energy ledger;
a two-rank cooling case closes total inventory to `1.8e-15` relative
([time integration](validation/eb-wall-sdc-integrator-20260927.md)).
MPI gas/particle SDC applications now expose embedded wall temperature, cap
thermal intervals and bind restart settings. A heated particle case, exact
rank-changing continuation and changed-temperature rejection pass; full
Spray-EB parity remains open
([application evidence](validation/eb-wall-application-20260927.md)).
A separate embedded viscous-traction kernel now includes normal bulk viscosity
and tangential shear, matching the compiled pinned function to `1.44e-14` in
24 cases
([traction kernel](validation/eb-wall-viscous-kernel-20260927.md)).
Stationary viscous traction now supplies gas momentum sources and a composite
MPI force ledger; cross-owner and covered-coarse exclusion checks pass
([force source](validation/eb-wall-viscous-source-20260927.md)).
Common and recursive MPI SDC now integrate this force with wall impulse and
zero stationary-wall mechanical energy. A two-rank adiabatic drag fixture closes
inventory to `3.1e-17` relative
([traction time integration](validation/eb-wall-viscous-integrator-20260927.md)).
MPI gas/particle applications now expose optional embedded no-slip, with a
viscous interval cap and restart identity. Combined and adiabatic particle
steps, exact rank-changing continuation and changed-setting rejection pass
([no-slip application](validation/eb-noslip-application-20260927.md)).
Serial CTU/MOL common and recursive steps now also include embedded heat and
viscous traction. Eight wall cases match MPI state/temperature exactly;
serial-to-MPI restart, wall-setting rejection and the wall timestep cap pass
([serial connection](validation/serial-eb-wall-20260927.md)).

A pinned four-port Spray-Multijet adaptation completes two SDC steps and eight
events with 32 parcels and closed fuel inventory. Original liquid-property
values are preserved; full-duration reference parity remains open
([case](../cases/spray_multijet_3d/README.md),
[execution evidence](validation/multijet-connection-20260927.md)).

The pinned reaction-free NC10H22/O2/N2 mechanism for Spray-Jet now imports
with Simple transport fits. Six-state compiled coefficient checks pass; a
decane jet reference comparison remains open
([import evidence](validation/decane-spray-import-20260927.md)).
A short source-rate-driven decane application now completes two injection
events and three SDC steps; eight parcels and gas fuel close the injected
mass balance ([connection](validation/decane-jet-connection-20260927.md)).
External MPI EB tree rate files now accept all 7,502 reference knots, with
exact numerical restart identity and a two-rank continuation check
([long profiles](validation/long-spray-rate-profiles-20260927.md)).
External event schedules now support 65,536 entries; a 1,100-entry bounded
application check and old/new schedule restart checks pass. Full positive-rate
jet execution remains open ([schedules](validation/spray-event-files-20260927.md)).

The pinned reaction-free NC7H16/O2/N2 mechanism imports with source/license
checks; selected serial spray and MPI EB tree SDC paths require chemistry OFF
and temperatures inside the source range
([import](validation/nonreacting-spray-import-20260927.md),
[application boundary](validation/nonreacting-spray-app-20260927.md)).
Selected EB tree gas/particle transport can use pinned Simple fits, with
coefficient provenance and restart identity checks
([integration](validation/simple-transport-integration-20260927.md)).

The periodic 4³ two-step Spray-Evaporation adaptation now uses explicit Euler
particle updates and the pinned single-cell source deposition. Particle
temperature/velocity differences are `4.0e-13 K` / `3.9e-14 m/s`; fuel-density
difference is `1.4e-11` relative to its reference maximum. Rank-changing restart
matches continuous output. Gas momentum and stored-temperature differences
remain, so full spray parity and physical validation are still open
([current comparison](validation/spray-cell-source-20260927.md)).
Fitted Simple transport also supplies bulk viscosity to 3D stress, viscous work
and timestep limits. Six-state compiled-reference checks pass; the remaining
axial-momentum difference is `2.7e-6` relative to the reference maximum
([bulk-viscosity evidence](validation/simple-bulk-viscosity-20260927.md)).
Fitted face coefficients now follow the pinned harmonic mean; three nonuniform
faces agree to `6.2e-16`, without explaining the remaining evaporation residual
([face checks](validation/simple-face-transport-20260927.md)).
The reference temperature inversion retains corrections below `1e-6 K`;
a compiled 64-state audit reproduces its `8.87e-7 K` stored-temperature/EOS
discrepancy. This bounds the temperature comparison without changing conserved
energy; the momentum residual remains unresolved
([inversion audit](validation/spray-temperature-inversion-20260927.md)).
A full reference build with tighter temperature inversion reduces the maximum
temperature difference to `3.75e-9 K`; momentum remains at `2.66e-6` relative.
This separates the temperature stopping rule from the remaining momentum
residual ([intervention](validation/spray-tight-temperature-reference-20260927.md)).
EB tree inputs now expose individual gas transport operators with unchanged
defaults; a two-rank continuation preserves the gas output and rejects changed
transport settings ([controls](validation/eb-tree-transport-controls-20260927.md)).

The fixed regular application now exposes a nonreacting `gamma_law` pseudo-air
model for HIT/SDC/LES, with independent ideal-gas state and transport-scale
checks for four LES closures ([evidence](validation/gamma-law-app-20260927.md)).
Its prescribed transport needs no invented molecular database; gamma-law
MPI/AMR support is not established. A full pinned PeleC HIT build now agrees
on a 4³ one-step comparison for hydro, constant transport and WALE (maximum
normalized difference `3.92e-15`; [evidence](validation/hit-full-reference-20260927.md)).
Long-time turbulence parity remains unestablished.
Dynamic Smagorinsky reaches `3.92e-15` with explicit pinned heat-component
compatibility. Its default face-normal heat flux differs from the pinned
reference's all-direction x-component selection; this explains the preceding
`7.67e-11` residual ([intervention](validation/hit-dynamic-heat-component-20260927.md)).
Its SI denominator regularizers also have an independent unit-conversion check
([details](validation/hit-dynamic-reference-20260927.md)).
An optional regular timestep cap now supports matched multistep reference
runs. Dynamic, Smagorinsky and Vreman agree over three steps on 4³ cells
(maximum normalized difference `7.72e-14`;
[checks](validation/hit-multistep-reference-20260927.md)).
The pinned state-dependent linear momentum forcing now reaches regular SDC
source stages and inventory accounting. Three-step hydro and transport/WALE
comparisons agree with the full reference (`3.92e-15` maximum normalized
difference; [checks](validation/hit-linear-forcing-20260927.md)).
HIT input also accepts an explicit coordinate-unit factor; the gamma-law
template converts centimetre reference lengths to metres, with exact
manufactured-input equivalence ([checks](validation/hit-coordinate-20260927.md)).

Periodic 3D transport face/flux and timestep APIs now accept explicit constant
viscosity and conductivity, with analytical flux and coefficient-scaling checks
([evidence](validation/constant-transport-20260927.md)). Periodic regular SDC
and fixed/selected applications now accept `reactive_constant_transport`, with
active coefficient, WALE coexistence, conservation and exact fixed/selected
output checks ([integration](validation/constant-transport-app-20260927.md)).
`mode='hit'` derives the coefficients from initialized HIT scales and supplied
Taylor length/Reynolds/Prandtl parameters, with independent coefficient and
explicit-input equivalence checks
([scaling evidence](validation/hit-transport-scales-20260927.md)).
With species diffusion disabled, prescribed transport no longer evaluates
unused molecular coefficients; thermodynamic validation remains active
([dependency checks](validation/constant-transport-independent-20260927.md)).

Primitive conversion now preserves exact forward-EOS endpoint states when
inverse rounding lands outside the mechanism interval. Strict ranges and
out-of-range rejection remain unchanged; endpoint unit checks and 300 K
fixed/selected uniform applications pass
([repair evidence](validation/eos-endpoint-20260927.md)).

Regular fixed/selected gas applications now accept pinned-format HIT velocity
cubes through `reactive_hit`, with ASCII/binary input, periodic interpolation
and Mach scaling. Manufactured interpolation and pinned 32³ input checks reach
WALE SDC with exact fixed/selected output
([evidence](validation/hit-input-20260927.md), [case](../cases/reactive_hit_3d/README.md)).
This is initialization/integration evidence, not evolved HIT reference parity.
The pinned C++ `pc_initdata`/`locate` functions also agree with the SI gamma-law
initial field on 64 cells (maximum normalized difference `3.77e-15`;
[direct kernel comparison](validation/hit-initialization-reference-20260927.md)).

The forced reactor now has a transactional CVODE BDF backend for changing
species densities and internal-energy sources. Inert analytical, early-time
native-reactor comparison and failure/reuse checks pass
([backend evidence](validation/forced-cvode-20260927.md)). Cell and regular-grid
SDC now accept an explicit reactor callback, with CVODE source/history and
rollback checks ([dispatch evidence](validation/sdc-cvode-bridge-20260927.md)).
Fixed/selected regular gas applications expose `sdc_reactor='cvode'` when built
with SUNDIALS, with exact matched-composition output checks
([application evidence](validation/sdc-cvode-app-20260927.md)).
Serial/MPI AMR common and recursive drivers now forward reactor callbacks;
focused native-adapter dispatch and rollback checks pass
([AMR dispatch evidence](validation/amr-sdc-reactor-20260927.md)).
Fixed/selected serial/MPI AMR applications now expose CVODE selection, with
exact common/recursive serial-MPI outputs and changed-rank restart checks
([AMR application evidence](validation/amr-sdc-cvode-app-20260927.md)).
Single-grid and serial/MPI EB common/recursive SDC now forward reactor
callbacks, with cut-cell native-adapter dispatch and rollback checks
([EB dispatch evidence](validation/eb-sdc-reactor-20260927.md)).
Fixed/selected single-grid EB applications now expose CVODE selection with
exact restart and backend-mismatch rejection. Integration also repaired forced
CVODE density closure and Newton trial-state handling
([EB application evidence](validation/eb-sdc-cvode-app-20260927.md)).
Fixed/selected serial/MPI EB tree and MPI particle applications now expose
CVODE selection with reactor-bound checkpoints. Focused gas common/recursive
and reacting particle recursive LES checks establish exact changed-rank
continuation ([tree/particle evidence](validation/eb-tree-particle-cvode-20260927.md)).

KHRT now connects reference-count/azimuth children, deterministic samples and
conservative gas momentum/energy compensation to application input through
`khrt_parcel_factor`. The default -1 retains four symmetric children; [0,1]
enables the reference count policy. Fixed-mechanism MPI EB common/recursive
SDC has conservation and changed-rank continuation checks
([application evidence](validation/khrt-reference-app-20260927.md)).

MPI AMR SDC now accepts owner-local diffusion callbacks for common and recursive
steps, with collective validation and combined-flux synchronization
([callback checks](validation/mpi-amr-sdc-diffusion-20260927.md)). The MPI LES
RHS and common/recursive SDC wrapper now connect all four closures, with
owner-local stage evaluation and temporal parent donors
([MPI LES checks](validation/les-mpi-amr-sdc-20260927.md)). The opt-in MPI gas
application now accepts `amr_tree_les`, limits the timestep and preserves LES
settings in checkpoint identity ([application checks](validation/les-mpi-amr-app-20260927.md)).
The serial gas application now uses the same input and checkpoint identity,
with exact serial/MPI output and checkpoint interchange in focused common and
recursive cases ([serial application checks](validation/les-serial-amr-app-20260927.md)).
MPI EB SDC now accepts additional diffusion rates and integrated face fluxes
for common and recursive steps ([callback checks](validation/eb-sdc-diffusion-callback-20260927.md)).
A regular-geometry LES adapter now reaches common/recursive MPI particle SDC;
all four closures have focused mixed-parcel conservation/rollback checks
([particle LES checks](validation/les-particle-sdc-20260927.md)). The particle
application now accepts `eb_spray_les`, applies a diffusion timestep bound and
binds LES settings to restart identity ([application checks](validation/les-particle-app-20260927.md));
non-regular geometry is rejected by the adapter.
The selected H2O2 gas and particle applications also have focused common and
recursive LES checks, including serial/MPI gas interchange and changed-rank
particle restart ([selected integration](validation/les-selected-integration-20260927.md)).

Regular fixed/selected gas input now selects periodic SDC with
external forcing through `time_integrator='sdc'`
([checks](validation/sdc-app-20260926.md)). Molecular transport now participates
as an old/new source family ([checks](validation/sdc-transport-20260926.md)).
Regular SDC also accepts additional source callbacks together with molecular
transport, preserving transport boundary inventories and exported fluxes
([focused checks](validation/sdc-callback-transport-20260927.md)).
The opt-in LES runtime now supplies a regular SDC wrapper for Smagorinsky,
WALE, Vreman and Dynamic sources alongside molecular transport and forcing
([library checks](validation/les-sdc-20260927.md)). Fixed/selected regular gas
applications now accept `reactive_les`, with diffusion timestep limits and
separate LES boundary inventories. Four-closure CVODE application checks give
exact fixed/selected CSVs and reject unavailable or incompatible selections
([application checks](validation/les-regular-app-20260927.md)).
The regular LES SDC wrapper now also forwards selectable forced reactors, with
four-closure CVODE dispatch, conservation and rollback checks
([reactor integration](validation/les-sdc-cvode-20260927.md)).
Regular LES SDC now exports combined hydro/transport/LES full-face fluxes for
subsequent AMR synchronization ([reconstruction checks](validation/les-sdc-flux-20260927.md)).
It also accepts caller-supplied stage halos and externally retained source
estimates for coordinated patch iterations
([halo checks](validation/les-sdc-halos-20260927.md)).
Serial AMR common-step SDC now connects all four LES closures to stage sources
and combined-flux reflux ([AMR checks](validation/les-amr-sdc-20260927.md)).
Serial recursive LES SDC now forwards parent temporal donors and combines
substep fluxes at synchronization
([recursive checks](validation/les-amr-sdc-recursive-20260927.md)).
The library couples chemistry into periodic SDC
([checks](validation/sdc-reactions-20260926.md)), and fixed/selected application
input now enables it with persistent reaction sources
([checks](validation/sdc-reacting-app-20260926.md)). Broader boundary/AMR/MPI
integration remains open.
The SDC library additionally supports zero-gradient outflow with separate
boundary inventories ([checks](validation/sdc-outflow-20260926.md)). Fixed/selected
application routing now exposes this path, including a nonuniform reacting
check ([evidence](validation/sdc-outflow-app-20260926.md)). Outflow molecular
transport now participates with separate boundary accounting
([checks](validation/sdc-outflow-transport-20260926.md)); other physical
boundary types remain open.
Prescribed outlet pressure now reaches both SDC hydro and transport
([checks](validation/sdc-pressure-outlet-20260926.md)); existing ghost-state
admissibility limits still apply.
Specified-state and Mach inflow now reach SDC hydro and molecular transport
([checks](validation/sdc-inflow-20260926.md)). Slip, no-slip and isothermal
walls now connect reflection, viscous stresses and heat transport
([checks](validation/sdc-walls-20260926.md)).

The serial AMR library now runs common-step SDC across the whole patch tree,
including hydro, molecular transport, forcing, reactions, reflux and final
restriction ([checks](validation/sdc-tree-loop-20260926.md)). All patches finish
one stage before the next begins. Face-flux export and caller-owned halo APIs
remain available ([fluxes](validation/sdc-flux-export-20260926.md),
[halos](validation/sdc-halo-callbacks-20260926.md)). Static serial fixed/selected AMR applications now expose this path with separate
external-source inventories ([checks](validation/sdc-tree-app-20260926.md)).
Serial SDC checkpoint/restart now preserves reaction estimates and forcing
inventories in plain/indexed records ([checks](validation/sdc-tree-restart-20260926.md)).
Scheduled/automatic serial SDC regridding now transfers signed reaction estimates
with the gas ([checks](validation/sdc-tree-regrid-20260926.md)).
SDC hydro stages now support level selection and separate coarse temporal
donors for state and predictor sources ([checks](validation/sdc-level-hydro-20260926.md)).
Transport stages now select levels and interpolate parent endpoints; reaction
stages update one level while retaining other levels' state and source estimates
([checks](validation/sdc-level-sources-20260926.md)). Coupled level SDC now advances
one interval using parent state/source endpoints and exports integrated fluxes
([checks](validation/sdc-level-loop-20260926.md)). Serial recursive SDC now connects refinement-ratio substeps, parent temporal
donors, integrated-flux synchronization and application restart
([checks](validation/sdc-subcycle-20260926.md)). MPI now has an owner-local SDC
reaction stage with collective rollback and optional level selection
([checks](validation/sdc-mpi-reaction-20260926.md)). MPI SDC hydro now communicates
signed predictor sources and returns owner-local rates/integrated fluxes
([checks](validation/sdc-mpi-hydro-20260926.md)). MPI transport also returns
owner-local SDC rates/fluxes with collective limiter reduction
([checks](validation/sdc-mpi-transport-20260926.md)). Common-step MPI SDC now
couples hydro, transport, forcing, reactions and reflux while retaining owner-local
state ([checks](validation/sdc-mpi-loop-20260926.md)). The MPI checkpoint library
now reads/writes SDC V3/V4 records with owner-local source payloads and changed-rank
interchange ([checks](validation/sdc-mpi-checkpoint-20260926.md)). MPI fixed/selected applications now expose common-step SDC, source-preserving
regridding/rebalancing and changed-rank restart
([checks](validation/sdc-mpi-app-20260926.md)). MPI SDC hydro now selects a level
and samples separate parent state/source temporal endpoints
([checks](validation/sdc-mpi-level-hydro-20260926.md)). MPI transport now also
selects levels and interpolates parent endpoints
([checks](validation/sdc-mpi-level-transport-20260926.md)). The coupled MPI level
SDC driver now uses these stages and exports interval fluxes and predictors
([checks](validation/sdc-mpi-level-loop-20260926.md)). Recursive MPI SDC now connects child substeps, accumulated-flux synchronization
and application dispatch ([checks](validation/sdc-mpi-subcycle-20260927.md)). EB now has a covered-cell-aware forced SDC reaction stage
([checks](validation/sdc-eb-reaction-20260927.md)); EB transport also exposes a transactional source/flux stage
([checks](validation/sdc-eb-transport-20260927.md)). EB CTU prediction now accepts conserved SDC source rates
([checks](validation/sdc-eb-hydro-source-20260927.md)). EB hydro now exposes conservative rates, integrated fluxes and boundary/wall
ledgers ([checks](validation/sdc-eb-hydro-stage-20260927.md)). Physical boundary filling now feeds the EB SDC hydro stage
([checks](validation/sdc-eb-boundary-stage-20260927.md)). Small-cell redistribution now has a transactional source-rate interface
([checks](validation/sdc-eb-redist-20260927.md)). The serial single-grid EB SDC loop now couples hydro, transport, forcing,
reactions and redistribution ([checks](validation/sdc-eb-loop-20260927.md)).
Fixed/selected single-grid applications now accept SDC input and report
separate forcing inventories ([checks](validation/sdc-eb-app-20260927.md)).
Single-grid EB SDC checkpoint/restart now preserves reaction estimates and
forcing inventories ([checks](validation/sdc-eb-restart-20260927.md)).
The EB tree now exposes provisional SDC reaction updates across all levels or
a selected level, preserving refined-parent estimates until synchronization
([checks](validation/sdc-eb-tree-reaction-20260927.md)).
Signed EB tree SDC source halos now support peer/ancestor donors, ancestor
time interpolation and scalar physical extrapolation
([checks](validation/sdc-eb-tree-source-halo-20260927.md)).
The common-step EB tree now exposes source-aware provisional hydro rates and
time-integrated fluxes without averaging, redistribution or reflux
([checks](validation/sdc-eb-tree-hydro-20260927.md)).
Provisional EB tree transport now returns rates, integrated fluxes and a
limiter without averaging or synchronization
([checks](validation/sdc-eb-tree-transport-20260927.md)).
The serial EB tree now couples those stages in a common-step SDC loop, with
final peer/coarse-fine synchronization and separate source inventories
([checks](validation/sdc-eb-tree-loop-20260927.md)).
The serial EB tree checkpoint envelope now preserves SDC forcing inventories
and reaction fields, with transactional reads
([checks](validation/sdc-eb-tree-checkpoint-20260927.md)).
The serial tree application now routes static common-step SDC input and
checkpoint/restart, with an independent forcing budget column
([checks](validation/sdc-eb-tree-app-20260927.md)).
EB tree SDC regridding now transfers signed reaction estimates together with
gas state; scheduled coarsen/refine restart is checked
([checks](validation/sdc-eb-tree-regrid-20260927.md)).
EB SDC hydro and transport stages now support a selected level with ancestor
time interpolation; hydro separately interpolates source predictors
([checks](validation/sdc-eb-tree-level-stages-20260927.md)).
The coupled EB tree SDC loop now also advances a selected level, retains
inactive levels and exports interval fluxes plus old/new source predictors
([checks](validation/sdc-eb-tree-level-loop-20260927.md)).
Recursive EB tree SDC now accumulates interval fluxes, synchronizes adjacent
levels and routes through the serial application with checkpoint/restart
([checks](validation/sdc-eb-tree-subcycle-20260927.md)).
MPI EB trees now expose an owner-local SDC reaction stage, including selected
levels and collective failure handling
([checks](validation/sdc-mpi-eb-reaction-20260927.md)).
MPI EB transport now returns owner-local SDC rates and integrated fluxes while
preserving provisional parent states, with collective rollback
([checks](validation/sdc-mpi-eb-transport-20260927.md)).
MPI EB SDC hydro now exchanges signed predictor sources and returns owner-local
CTU rates and integrated fluxes, retaining provisional parent states
([checks](validation/sdc-mpi-eb-hydro-20260927.md)).
Common-step MPI EB SDC now couples reactions, old/new transport and external
sources, CTU hydro, final reflux and composite inventories
([checks](validation/sdc-mpi-eb-loop-20260927.md)). MPI EB SDC checkpoint V4 now stores owner-local reaction estimates and forcing
inventories with direct patch offsets, reads serial V3 records and supports
serial readers and changed-rank reads
([checks](validation/sdc-mpi-eb-checkpoint-20260927.md)). MPI EB application routing now supports fixed-topology common-step SDC with
rank-changing checkpoint continuation
([checks](validation/sdc-mpi-eb-app-20260927.md)). MPI EB SDC load rebalancing now migrates reaction estimates atomically with gas
and temperature ([checks](validation/sdc-mpi-eb-rebalance-20260927.md)). MPI EB SDC regridding now transfers signed reaction estimates with fluid-volume
restriction and PCM prolongation; scheduled changed-rank restart is checked
([evidence](validation/sdc-mpi-eb-regrid-20260927.md)). MPI EB SDC source halos now interpolate ancestor time endpoints with collective
control agreement and rollback
([checks](validation/sdc-mpi-eb-source-time-20260927.md)). MPI EB SDC hydro/transport stages now select a level and receive ancestor gas
endpoints; hydro additionally accepts separate source endpoints and fractions
([checks](validation/sdc-mpi-eb-level-stages-20260927.md)). The MPI EB SDC loop now advances one selected level, forwards distinct gas/source
times and exports fluxes plus old/new predictor fields
([checks](validation/sdc-mpi-eb-level-loop-20260927.md)). MPI EB SDC now recursively advances child levels, accumulates their fluxes and
synchronizes each parent interval. MPI application subcycling, coupled transport
timestep bounds and rank-changing continuation are connected
([checks](validation/sdc-mpi-eb-subcycle-20260927.md)). A serial composite flux-redistribution
operator transfers signed rates across peer, coarse/fine and paired periodic
leaf faces. Common-step serial SDC now applies it to hydro, molecular transport
and the final reflux increment
([operator](validation/sdc-eb-tree-redistribution-20260927.md),
[stage checks](validation/sdc-eb-tree-redistribution-loop-20260927.md)).
Serial application input and checkpoint identity now expose `flux_redist`
for common-step SDC ([checks](validation/sdc-eb-tree-redistribution-app-20260927.md)).
MPI common-step SDC and application dispatch now use owner-local composite
flux redistribution through sparse neighbor-rate and excess exchanges;
rank-changing continuation retains the selector
([checks](validation/sdc-mpi-eb-redistribution-20260927.md)). MPI particle intervals now export independent owner-local integrated gas sources
with collective rollback ([checks](validation/spray-sdc-source-export-20260927.md));
the serial/MPI EB SDC gas driver now consumes old/new auxiliary source callbacks
with a separate exchange inventory
([checks](validation/sdc-eb-auxiliary-source-20260927.md)). The MPI common-step particle wrapper now uses one old half-kick/full-drift
and reevaluates each new half-kick from that saved endpoint, publishing particle
state only with accepted gas SDC
([checks](validation/spray-sdc-kick-drift-20260927.md)). Physical gas boundary
data and final computational particle splitting now reach this wrapper
([checks](validation/spray-sdc-boundaries-split-20260927.md)). Independent uniform
particle acceleration now contributes accepted impulse/work to the external
inventory ([checks](validation/spray-sdc-force-20260927.md)). The MPI spray
application now dispatches common-step SDC and preserves reaction history and
external inventories through rank-changing coupled restart
([checks](validation/spray-sdc-app-20260927.md)). Film integration and recursive
particle scheduling remain open. MPI AMR/EB film intervals now export independent
integrated gas exchange for that connection
([checks](validation/film-sdc-source-export-20260927.md)); the common-step wrapper
now evolves existing stationary films from saved midpoint states
([checks](validation/spray-sdc-film-20260927.md)). Optional EB impact deposition/splash now uses the shared private particle-ID
transaction ([checks](validation/spray-sdc-impact-20260927.md)). The MPI spray
application now dispatches this film/impact path and preserves it through
regridding and rank-changing restart
([checks](validation/spray-sdc-impact-app-20260927.md)). Serial/MPI recursive
gas SDC now exposes timed per-level auxiliary sources and a transactional
exchange inventory ([checks](validation/sdc-recursive-auxiliary-source-20260927.md));
optional routing now retains off-level exchanges until coarse-interval
synchronization ([checks](validation/sdc-cross-level-source-20260927.md)).
Recursive MPI EB spray application input now connects particle impacts, splash
and newborn-film thermal updates with focused two/four-rank conservation checks
([evidence](validation/spray-recursive-impact-app-20260927.md)). A film-bearing
two-to-four-rank restart exactly matches continuous gas/particle/film output
([focused check](validation/spray-recursive-film-restart-20260927.md)). A serial
composite weighted-state operator now supports outward receivers and connected
neighborhood fallback across coarse/fine leaf faces, with conservation and
single-grid comparisons
([checks](validation/eb-tree-state-redistribution-20260927.md)); serial common-step SDC and application input now expose `state_redist` with
checkpoint-bound target fraction
([checks](validation/sdc-eb-tree-state-redistribution-20260927.md)). MPI now matches
serial common-step StateRedist through sparse value/mean exchanges and supports
rank-changing application restart
([checks](validation/sdc-mpi-eb-state-redistribution-20260927.md)). Conservative
neighborhood displacements now supply selected-level StateRedist transfers;
serial/MPI recursive SDC and rank-changing application restart use them
([checks](validation/sdc-state-redistribution-recursive-20260927.md)). Serial/MPI
source-level redistribution now separates immediate rates from time-integrated
cross-level transfers with a checked composite reconstruction
([checks](validation/sdc-eb-redistribution-level-transfer-20260927.md)); the
selected-level serial/MPI SDC loop now emits the accepted final-iteration
transfer alongside integrated fluxes
([checks](validation/sdc-eb-redistribution-level-loop-20260927.md)). The serial recursive
driver now accumulates transfers until their coarser source/recipient interval
closes, redistributes reflux increments and exposes application restart
([checks](validation/sdc-eb-redistribution-recursive-20260927.md)). MPI recursive
interval registers and application dispatch now follow the same owner-local
schedule, with exact serial/MPI state comparison and rank-changing continuation
([checks](validation/sdc-mpi-eb-redistribution-recursive-20260927.md)).

The variable-density forced chemistry RHS is now available as a library
building block ([checks](validation/forced-reactor-rhs-20260926.md)). An adaptive
implicit substep now integrates these rates
([checks](validation/forced-reactor-integrator-20260926.md)). Transactional
full-interval subcycling is also implemented
([checks](validation/forced-reactor-interval-20260926.md)) and feeds this library coupling.

Regular fixed/selected gas applications now accept force-per-volume input
`external_forcing(3)` with PPM/PeleC, time-centered work and separate source
inventories ([checks](validation/external-force-app-20260926.md)).

The full twelve-solve PPM CTU assembly now agrees with compiled pinned PeleC
bodies in eighteen small-grid cases, including independent mixture variation,
active shock flattening, artificial viscosity and inflow/outflow/slip/no-slip
boundaries ([comparison](validation/ppm-pinned-boundary-20260926.md)). This
includes normal tracing, Riemann solves, transverse correction, inflow PLM
replacement and final species normalization. Complete solver and physical
validation remain open.

Optional primitive sources now reach the periodic and supplied-halo PPM CTU
flux APIs, with 28 compiled-reference assembly cases
([evidence](validation/ppm-ctu-source-20260926.md)). Cell-local conserved-source
conversion now follows the pinned policy
([checks](validation/ppm-source-conversion-20260926.md)). Conserved source grids
now feed periodic/halo CTU directly, with 38 compiled-reference assembly cases
([checks](validation/ppm-conserved-source-20260926.md)); application source
generation and integration remain open.

The pressure path preserves independent and signed face species weights and
uses contact-speed passive upwinding
([projection](validation/ppm-mixture-projection-20260926.md),
[transverse](validation/ppm-signed-transverse-20260926.md),
[face flux](validation/ppm-passive-flux-20260926.md)). Physical wall/outflow
Riemann corrections reach regular serial/MPI AMR and MPI spray dispatch, while
explicit pressure outlets retain their prescribed donors
([checks](validation/ppm-physical-riemann-20260926.md)).

EB PLM normal prediction also preserves independently limited species weights
at open faces, with final CTU species-flux normalization retained
([regression](validation/eb-plm-species-20260926.md)).

Serial single-level EB KH eligibility now previews the complete boundary flight,
including periodic wraps and reflecting domain walls. Fixed-mechanism periodic
and wall cases preserve children and pending inventories through restart
([checks](validation/serial-eb-boundary-khrt-20260926.md)).

The serial EB impact interval now supports periodic crossings and specular
domain reflection before EB deposition, using the shared periodic film-image
counter ([checks](validation/serial-eb-periodic-impact-20260926.md)). Single-level
coupling and fixed/selected application restart now connect periodic film
exchange and reflecting-domain EB impacts
([evidence](validation/serial-eb-periodic-film-application-20260926.md)).

MPI EB deposition/splash and stationary films can now coexist with specular
outer-domain walls. Remaining flight and separate exterior/EB budgets are
preserved through regrid and changed-rank continuation
([checks](validation/mpi-eb-film-domain-reflection-20260926.md)). Periodic film thickness queries now sum all overlapping images across ranks
([checks](validation/mpi-periodic-film-query-20260926.md)). Periodic EB impact
flight and film/gas exchange now connect through regrid and changed-rank restart
([evidence](validation/mpi-eb-periodic-film-20260926.md)). Film creation on outer
domain faces remains open.

MPI EB KH child generation now accepts explicit nonperiodic domain walls,
with wall-before-release ordering and CTU regrid/changed-rank continuation
([checks](validation/mpi-eb-wall-khrt-20260926.md)). Explicit-boundary PCM
tree checkpoints now use a scoped embedded-record context, with changed-rank
continuation checked ([repair](validation/eb-tree-pcm-checkpoint-20260926.md)).

Patch Euler transport now accepts isothermal wall temperatures with transactional
heat-flux checks ([evidence](validation/eb-isothermal-patch-20260926.md)).
Serial common/recursive tree temperature forwarding now passes heated-wall
conservation checks ([evidence](validation/eb-tree-isothermal-transport-20260926.md)).
MPI common/recursive transport now agrees on wall temperatures and matches
serial ([checks](validation/mpi-eb-isothermal-transport-20260926.md)). Thermal
timestep bounds now connect in serial with an active-bound check
([evidence](validation/eb-tree-isothermal-timestep-20260926.md)); MPI thermal
selection now matches serial
([checks](validation/mpi-eb-isothermal-timestep-20260926.md)). Full reacting
common/recursive splits now forward wall temperatures with serial/MPI state
parity and rollback checks
([evidence](validation/eb-tree-isothermal-full-20260926.md)). Fixed/selected
serial and MPI gas applications now accept thermal wall inputs and bind wall
temperatures to changed-rank restart
([checks](validation/eb-tree-isothermal-application-20260926.md)). MPI EB
spray also forwards thermal gas walls, with common/recursive conservation and
fixed-mechanism changed-rank continuation checks
([evidence](validation/mpi-eb-isothermal-spray-20260926.md)).

MPI EB particle drift now supports specular domain-wall reflection with exterior
impulse accounting ([checks](validation/mpi-eb-wall-particle-20260926.md)).
Coupled common/recursive CTU wall steps and the fixed MPI spray application
now include wall reflection and changed-rank restart checks
([evidence](validation/mpi-eb-wall-spray-application-20260926.md)).

Fixed/selected tree gas applications now expose slip/no-slip adiabatic walls,
with scheduled regridding and changed-rank restart checked on an enclosed
planar fixture ([evidence](validation/eb-tree-wall-application-20260926.md)).

Serial EB tree halo construction now accepts adiabatic slip/no-slip walls,
including reflected six-layer corner donors with nonzero momentum
([checks](validation/eb-tree-wall-halo-20260926.md)). MPI wall halo exchange
now matches serial with distributed and empty-rank ownership
([checks](validation/mpi-eb-wall-halo-20260926.md)). Serial common PCM/CTU and
recursive PCM hydro now accept tangential adiabatic walls with composite
conservation checks ([evidence](validation/eb-tree-wall-hydro-20260926.md));
MPI common/recursive CTU wall hydro now matches serial for a tangential planar
fixture ([checks](validation/mpi-eb-wall-hydro-20260926.md)). Full reacting
wall splits now match serial for common/recursive CTU with explicit/implicit
chemistry ([checks](validation/mpi-eb-wall-full-20260926.md)). Serial common/recursive molecular transport now forwards
domain-wall labels and enforces impermeable adiabatic fluxes
([checks](validation/eb-tree-wall-transport-20260926.md)); MPI wall transport runtime is covered by those full-split checks.

Serial EB tree halos also accept paired periodic faces, wrapping donor indices
through peers and ancestors, including halos wider than the root domain
([checks](validation/eb-tree-periodic-halo-20260926.md)). MPI periodic tree
halo exchange now matches serial face-only, full one-layer and six-layer
construction, including temporal ancestor donors and idle ranks
([checks](validation/mpi-eb-periodic-halo-20260926.md)). Periodic tree application
dispatch is not yet connected.
The shared EB peer-flux correction kernel now accepts a periodic domain
length on one axis, checks matching face centroids and applies conservative
state correction across the seam
([checks](validation/eb-periodic-peer-flux-20260926.md)). Serial tree traversal
now reconciles periodic peers and mirrors restricted fine transfers into the
opposite coarse face. Common-step PCM and PLM/CTU hydro connect that path
([checks](validation/eb-tree-periodic-flux-20260926.md)). Serial common and
recursive molecular transport and recursive hydro now forward periodic
controls, including provisional endpoint correction before descendant time
interpolation ([checks](validation/eb-tree-periodic-recursive-20260926.md)).
Periodic common/recursive reacting PLM/CTU splits now pass focused explicit
and implicit chemistry checks. MPI flux synchronization also accepts periodic
faces with exact serial agreement, owner-local transfers and collective
rollback ([checks](validation/eb-tree-periodic-full-mpi-flux-20260926.md)).
Common and recursive MPI hydro now forward periodic labels through provisional
endpoint and final flux synchronization, with focused PLM/CTU serial agreement
and collective rollback checks
([checks](validation/mpi-eb-periodic-hydro-20260926.md)). MPI common/recursive
transport and reacting splits also connect periodic synchronization, with
serial agreement and collective rollback checks
([checks](validation/mpi-eb-periodic-full-20260926.md)). Fixed/selected serial
and MPI gas applications now accept paired periodic input, wrap geometry
padding and retain boundary identity across scheduled regridding and changed-rank
restart ([checks](validation/eb-tree-periodic-application-20260926.md)).
MPI EB spray also accepts periodic input, with fixed-mechanism 2/4-rank
comparison and 2-to-4-rank restart checks
([checks](validation/mpi-eb-periodic-application-20260926.md)). Periodic
EB impact/film combinations are also connected, as described above.
The composite EB particle stencil now accepts periodic axes and checks
visibility along segmented, wrapped paths to donor centroids
([checks](validation/eb-tree-periodic-stencil-20260926.md)). MPI particle field
sampling and gas deposition now route those periodic donors to their owners,
with conservation and rollback checks
([checks](validation/mpi-eb-periodic-field-20260926.md)). Mixed-particle drift
and reacting MPI EB spray coupling now forward periodic faces through field
exchange, drift and migration, with common/recursive CTU conservation checks
([checks](validation/mpi-eb-periodic-step-20260926.md)). Periodic KH child
generation now follows wrapped parent flight and passes changed-rank continuation
([evidence](validation/mpi-eb-periodic-khrt-20260926.md)).

Serial EB tree halo construction now accepts prescribed inflow states for
one-cell and six-cell halos while retaining interior peer/ancestor donors
([checks](validation/eb-tree-prescribed-halo-20260926.md)). Common and recursive
hydro/transport/reacting library updates now forward those conditions
([checks](validation/eb-tree-inflow-step-20260926.md)). MPI halo exchange also
accepts collectively agreed prescribed states, with serial agreement and idle-rank
rollback checks
([checks](validation/mpi-eb-prescribed-halo-20260926.md)). Common and recursive
MPI hydro now forward these controls with serial agreement
([checks](validation/mpi-eb-inflow-hydro-20260926.md)). MPI molecular transport,
full reacting splits and donor timestep bounds are also connected
([checks](validation/mpi-eb-inflow-full-20260926.md)). Fixed/selected serial and
MPI gas application inputs now connect prescribed inflow and boundary-bound
changed-rank restart across scheduled regridding
([checks](validation/eb-tree-inflow-application-20260926.md)). Tree application
inlets also accept per-face named mole fractions independently of initial gas,
with resolved composition bound to restart
([checks](validation/eb-tree-inlet-composition-20260926.md)). MPI EB spray now
also forwards prescribed inflow through the coupled split and application,
with fixed-mechanism changed-rank continuation checks
([checks](validation/mpi-eb-spray-inflow-20260926.md)).

Single-level fixed/selected EB gas applications now accept named mole fractions
per inlet face, with exact input-composition restart binding
([checks](validation/eb-inlet-composition-20260926.md)). Single-level spray
also constructs these inlet states and binds composition to its coupled restart
([checks](validation/eb-spray-inlet-composition-20260926.md)).

Single-level fixed/selected EB applications now offer `hydro_method='plm_ctu'`
with the `pelec` Riemann solver, outflow halos, reaction/transport composition
and method-bound checkpoint continuation. PCM remains the default. Fixed and
selected H2O2 coupled continuation passed focused checks
([evidence](validation/eb-ctu-application-20260925.md),
[input example](../cases/reactive_eb_3d/plm_ctu.nml)). This does not qualify a physical wall-flow solution. Explicit inflow/outflow input,
donor timestep bounds and boundary-bound restart are also connected for fixed and
selected single-level gas applications
([checks](validation/eb-inflow-application-20260925.md)). Fixed domain-wall inputs also
connect reflected geometry, thermal timestep limits and temperature-bound restart
([checks](validation/eb-wall-application-20260925.md)). Paired periodic faces also connect
wrapped geometry, coupled transport/CTU and changed-boundary restart rejection
([checks](validation/eb-periodic-application-20260925.md)).
MPI EB spray also forwards CTU through the gas/particle split, with common and
recursive library checks and fixed-mechanism changed-rank application continuation
([evidence](validation/mpi-eb-spray-ctu-20260925.md)). Serial EB spray fixed/selected
applications also connect CTU with method-bound continuation
([evidence](validation/serial-eb-spray-ctu-20260925.md)). The allocated serial
EB spray library also forwards explicit inflow/outflow gas boundaries and
escaped particle inventories
([checks](validation/eb-spray-open-coupling-20260925.md)). Fixed/selected serial
spray input, inlet timestep bounds and boundary-bound restart are connected
([application checks](validation/eb-spray-inflow-application-20260925.md)). Periodic
serial spray library coupling now wraps field exchange and parcel drift with
matching gas faces
([checks](validation/eb-spray-periodic-field-20260925.md)). Fixed/selected serial
periodic input and restart are also connected
([checks](validation/eb-spray-periodic-application-20260925.md)); periodic
impact/film dispatch now has the focused checks linked above.
Serial domain-wall gas conditions now connect specular particle reflection,
wall thermal timestep bounds and fixed/selected temperature-bound restart
([checks](validation/eb-spray-wall-20260926.md)).
Serial/MPI fixed and selected EB tree gas applications now also accept CTU,
including scheduled regridding and method-bound changed-rank continuation
([application checks](validation/eb-tree-ctu-application-20260925.md)).
Serial EB tree library operators support common and recursive CTU steps with
reaction/transport composition ([checks](validation/eb-tree-ctu-subcycle-20260925.md)).
MPI common and recursive tree hydro and coupled reaction/transport stepping
also support CTU with owner-local state and serial agreement
([common hydro](validation/mpi-eb-ctu-common-20260925.md),
[recursive/coupled checks](validation/mpi-eb-ctu-recursive-20260925.md)).

The pinned nozzle Mach-inlet state kernel is now available for six faces, with
interior composition and extrapolated or prescribed pressure
([state-kernel checks](validation/mach-inlet-state-20260925.md)). Serial hydro
now constructs its ghosts at both SSPRK2 stages
([hydro checks](validation/mach-inlet-hydro-20260925.md)). Molecular transport
also rebuilds Mach inlet states per stage
([transport checks](validation/mach-inlet-transport-20260925.md)). The full serial
reaction/transport/hydro composition forwards the same controls
([split checks](validation/mach-inlet-split-20260925.md)). Regular fixed/selected
input and per-step inlet timestep limits are connected
([input checks](validation/mach-inlet-input-20260925.md)). Regular Smagorinsky,
WALE, Vreman and Dynamic LES operators also accept the Mach controls
([LES checks](validation/mach-inlet-les-20260925.md)). Regular fixed/selected spray
now forwards Mach input through its coupled operators and binds restart context
([spray checks](validation/mach-inlet-spray-20260925.md)). Distributed wide Mach
halos and regular MPI Dynamic LES now match serial operators
([MPI LES checks](validation/mach-inlet-mpi-les-20260925.md)). The MPI spray
two-layer ghost entry also supports inviscid and viscous Mach donors
([ghost checks](validation/mach-inlet-mpi-spray-ghosts-20260925.md)). MPI hydro,
molecular transport, static LES and mixture splitting now forward the controls
([step checks](validation/mach-inlet-mpi-stepping-20260925.md)). Regular MPI spray
input, timestep limits and rank-changing restart are connected
([application checks](validation/mach-inlet-mpi-application-20260925.md)). The
abrupt-start timestep restriction from unused Dynamic filter padding is corrected
([startup checks](validation/dynamic-active-face-timestep-20260925.md)). The AMR
ghost sampler now reconstructs Mach inlet states after donor interpolation
([ghost checks](validation/mach-inlet-amr-ghosts-20260925.md)). MPI AMR ghost
exchange forwards the same controls with collective agreement and rollback
([MPI ghost checks](validation/mach-inlet-mpi-amr-ghosts-20260925.md)). Serial AMR
hydro now rebuilds Mach inlets for common and recursive level steps
([hydro checks](validation/mach-inlet-amr-hydro-20260925.md)). Distributed common
and recursive hydro now forward collectively agreed Mach controls
([MPI hydro checks](validation/mach-inlet-mpi-amr-hydro-20260925.md)). AMR
serial molecular transport also forwards Mach controls through common and
recursive stages ([transport checks](validation/mach-inlet-amr-transport-20260925.md)).
MPI molecular transport now also forwards the controls through common and
recursive stages ([MPI transport checks](validation/mach-inlet-mpi-amr-transport-20260925.md)).
Serial AMR LES selectors and common/recursive stepping also accept Mach controls
for all four models ([LES checks](validation/mach-inlet-amr-les-20260925.md)).
MPI AMR LES selectors and common/recursive stepping now forward collectively
agreed controls ([MPI LES checks](validation/mach-inlet-mpi-amr-les-20260925.md)).
Serial/MPI reacting transport/hydro composition now forwards Mach controls through
common and recursive steps ([split checks](validation/mach-inlet-amr-split-20260925.md)).
AMR spray coupling forwards the controls through the gas composition and both
LES half steps ([spray checks](validation/mach-inlet-amr-spray-20260925.md)).
Serial fixed/selected AMR gas application input, current-donor timestep limits
and restart checks are connected
([application checks](validation/mach-inlet-amr-application-20260925.md)). MPI gas
input and collective inlet bounds also support rank-changing and serial/MPI restart
([MPI application checks](validation/mach-inlet-mpi-amr-application-20260925.md)).
AMR spray input, timestep selection and context-bound restart are also connected
([spray application checks](validation/mach-inlet-amr-spray-application-20260925.md)).
This completes Mach dispatch for these regular and patch-tree gas/spray paths;
it does not establish physical nozzle validation or complete PeleC parity.

Periodic Dynamic LES is now connected to serial/MPI regular patch-tree stepping
and the MPI AMR reacting-spray application, including filter-bound restart.
Focused fixed-mechanism MPI parity and 2-to-4-rank continuation passed on
2026-09-25 ([runtime evidence](validation/mpi-dynamic-amr-runtime-20260925.md)).
Adaptive regridding/rebalancing with Dynamic LES and selected H2O2 restart now
also have focused execution evidence ([combined checks](validation/les-dynamic-adaptive-20260925.md)).
Prescribed reservoir inflow and extrapolated outflow also run through Dynamic
AMR LES, including changed-rank spray restart
([focused checks](validation/les-dynamic-open-20260925.md)).
Dynamic wall faces are connected through the AMR selector, stepper and spray application
([wall integration](validation/les-dynamic-amr-walls-20260925.md)).
Dynamic pressure-target outlets now correct all donor layers using physical
domain coordinates ([pressure integration](validation/les-dynamic-amr-pressure-20260925.md)).
Serial LES subcycling now has a recursive driver with interval-weighted reflux
([focused checks](validation/les-amr-subcycle-20260925.md)) and a
[scaled timestep selector](validation/les-subcycle-timestep-20260925.md).
MPI LES now also schedules recursive owner-local intervals and synchronizes
accumulated fluxes ([MPI recursion checks](validation/les-mpi-amr-subcycle-20260925.md)).
The MPI spray application now connects recursive LES half steps and the scaled
diffusion selector ([coupling checks](validation/les-mpi-spray-subcycle-20260925.md)).
The regular serial spray application also connects Dynamic physical-boundary
fluxes and restart ([wall/restart check](validation/regular-dynamic-spray-20260925.md)).
Regular MPI spray also connects Dynamic physical-boundary exchange, with
changed-rank wall/restart evidence
([MPI application check](validation/mpi-dynamic-spray-walls-20260925.md)).
Detailed-fuel Dynamic LES and physical qualification remain open.
Non-regular EB LES is an extension beyond the pinned reference, which rejects
such cells in all four closures ([source audit](validation/les-reference-scope-20260925.md)).
Use the regular-grid AMR LES path for the supported LES scope.

| Area | Implemented subset | Still open |
| --- | --- | --- |
| Reactive flow and transport | 1D/2D paths; regular 3D NASA7 hydro/chemistry/molecular transport with per-face periodic, prescribed reservoir inflow, outflow, slip-wall and adiabatic/isothermal no-slip boundaries with prescribed tangential wall motion | Characteristic inflow and moving geometry, general 3D boundary/coupling qualification |
| 3D AMR | Serial/MPI two-level and multi-patch, multi-level trees; prescribed or temperature-tagged regridding, conservative remap, physical boundaries, common or recursive reaction/transport/hydro, measured-work rebalancing and rank-neutral restart ([tree design](design_decisions/amr-patch-tree-3d.md)) | Broader EB/particle qualification, large-tree scaling and convergence |
| 3D embedded boundary | Planar, extruded-section and 3D nodal level-set single-level applications with fixed/selected reaction/transport/hydro, outward-neighbor StateRedist with connected-small-cell fallback and geometry-bound restart ([redistribution evidence](validation/eb-cluster-redistribution-20260917.md), [application evidence](validation/eb-level-set-application-20260917.md)) | Broader small-cell stabilization, convergence and integrated EB-AMR-MPI qualification |
| Chemistry | Fixed H2/O2 and configure-time selected supported bundles, including Lindemann/Troe/SRI/Tsang falloff, chemically activated rates, PLOG and Chebyshev; official-Fortran CVODE for selected 0D and opt-in regular 3D spray cells; FFCM1_Red methanol, GRI-Mech and Lu n-dodecane fixtures; input/runtime capacity up to 128 species | Diesel ignition/physical validation across conditions and sparse/performance qualification |
| Parallelism/output | Rank-parity and changed-rank restart for qualified subsets; distributed two-level MPI AMR CSV output | Production scaling, broader distributed I/O qualification and other output paths |
| LES and spray | Opt-in regular 3D single-component parcels with per-port liquid identity or shared or per-port multicomponent compositions, density/pressure fits and vapor-state caps, and periodic/reservoir-inflow/outflow/specular parcel boundaries and adiabatic/isothermal tangentially moving gas walls; temperature-dependent liquid density/heat capacity and per-port breakup viscosity/surface tension, selectable Schiller-Naumann or PelePhysics drag with ETAB deformation correction, drag/heating/evaporation with Clausius or bounded Antoine saturation pressure and an optional near-boiling vapor-state cap, and ETAB and KH-RT representative-drop breakup, optional KH shed-mass children with conservative transverse velocities, and bounded computational splitting, conservative two-way coupling and uniform external acceleration with impulse/work ledger, deterministic cone/scheduled injection from point or finite circular solid/hollow outlets with swirl with up to 16 independently configured ports and optional piecewise-linear shared or independent per-port injection-rate histories with optional discharge-based injection speed, and per-port mass- or number-CDF Weibull/normal/lognormal/ChiSquared size distributions, deviatoric Smagorinsky, coupled restart | PeleMP latent/boiling-property fidelity, physical mixture boiling, general dispersion/collision, broader characteristic-boundary and near-wall turbulence qualification, calibrated properties and experimental validation; lateral film transport is a future extension |

MPI AMR particle/gas feedback, film ownership and changed-rank application
restart are implemented ([film integration](validation/spray-amr-film-application-20260918.md));
MPI EB spray also connects transactional regridding and continuation
([EB integration](validation/spray-mpi-eb-regrid-application-20260918.md)). These
are historical focused results, not newly rerun qualification. A pinned-source
[property audit](validation/spray-property-scope-20260925.md) identifies the
remaining PeleMP latent/boiling distinction; GCM is restricted upstream to
PeleLMeX. The current constant-latent surface cap remains an approximation.

A library API now coarsens 3D EB geometry while preserving physical moments and
face closure; it rejects disconnected or unrepresentable coarse cells
([evidence](validation/eb-geometry-coarsening-20260917.md)). Reactive volume-weighted restriction and PCM prolongation now preserve conserved
components and recover EOS temperatures transactionally
([evidence](validation/eb-state-transfer-20260917.md)). Open-area flux restriction
and interior coarse/fine reflux now have conservative transactional library APIs
([evidence](validation/eb-reflux-20260917.md)). EB Riemann fluxes now accept face
halos, and plain Euler returns accepted flux integrals
([evidence](validation/eb-face-halo-20260917.md)). Geometry extraction and PCM
coarse-donor/outflow face halos now feed fine-patch Euler updates
([evidence](validation/eb-patch-halo-20260917.md)). A serial one-patch, two-level
Euler hydro API now connects both updates, reflux, averaging and separate composite
boundary/wall inventories transactionally
([evidence](validation/eb-two-level-hydro-20260917.md)). Masked explicit/implicit
reaction halves now wrap that update with composite chemical inventories
([evidence](validation/eb-two-level-chemistry-20260917.md)). Molecular-transport
fluxes now accept selectable face donors, including the species limiter
([evidence](validation/eb-transport-halo-20260917.md)). A synchronized two-level
transport step now connects limited fluxes, reflux and averaging
([evidence](validation/eb-two-level-transport-20260917.md)). A combined
reaction/transport/hydro step now accumulates all inventories transactionally
([evidence](validation/eb-two-level-full-20260917.md)). A fixed-mechanism serial
[application case](../cases/reactive_eb_two_level_3d/README.md) now runs this path
and writes composite state/inventory CSVs
([evidence](validation/eb-two-level-application-20260917.md)). Existing validated
EB checkpoint records can now share a caller-owned formatted stream
([evidence](validation/eb-embedded-checkpoint-records-20260917.md)). The two-level
application now saves both patches and cumulative budgets together, with staged
restore and exact focused continuation
([evidence](validation/eb-two-level-restart-20260917.md)). Rejected application
steps now retry with bounded halving and transactional state/budget preservation
([evidence](validation/eb-two-level-retry-20260917.md)). A plain Euler hydro
CFL library estimate now includes cut-cell fluid volume, open faces, embedded
wall area and optional remote signal bounds
([evidence](validation/eb-euler-cfl-20260917.md)). A plain Euler diffusion
interval now scales open-face conductance by fluid volume and accepts a remote
diffusivity bound ([evidence](validation/eb-patch-diffusion-cfl-20260917.md)).
The application now optionally selects a common interval from both grids,
with donor bounds and mode-aware restart
([evidence](validation/eb-two-level-auto-dt-20260917.md)). A library hydro path
now advances ratio-based fine substeps with coarse-time face donors, accumulated
fluxes and final reflux ([evidence](validation/eb-subcycle-hydro-20260917.md)).
The transport library now uses ratio-squared fine steps and accumulates
accepted limited fluxes ([evidence](validation/eb-subcycle-transport-20260917.md)).
Reaction halves now wrap both subcycled operators, and the application connects
scaled interval selection and mode-aware restart
([evidence](validation/eb-subcycle-full-20260917.md)). Same-level peer face
halos now have a transactional overlay API with shared-metric checks
([evidence](validation/eb-peer-halo-20260917.md)). An owner/receiver peer-flux
API now corrects receiver inventories and reconciles shared physical transfer
transactionally ([evidence](validation/eb-peer-flux-20260917.md)). A common-step
hydro library now advances multiple fine patches with deterministic face
ownership, peer synchronization and coarse reflux
([evidence](validation/eb-multipatch-hydro-20260917.md)). Molecular transport
now follows the same multiple-patch transaction with limited-flux reconciliation
([evidence](validation/eb-multipatch-transport-20260917.md)). Composite-masked
reaction halves now wrap the multiple-patch hydro/transport operators with
chemical inventories and a single state/diagnostic transaction
([evidence](validation/eb-multipatch-full-20260917.md)). The serial fixed-mechanism
application now connects static multiple patches, common-interval selection,
composite CSV output and transactional checkpoint continuation
([evidence](validation/eb-multipatch-application-20260917.md)). Multiple-patch
subcycling now synchronizes peer halos/fluxes on every fine step and accumulates
fluxes for final coarse reflux, with mode-bound restart
([evidence](validation/eb-multipatch-subcycle-20260917.md)); MPI EB remains open. A selected-mechanism executable now shares
the application path and binds bundle/composition identity on restart
([evidence](validation/eb-selected-two-level-20260917.md)). Generated bundles
now expose source element composition and a generic EB inventory API consumes
it ([evidence](validation/eb-element-metadata-20260917.md)). The two-level
selected frontend now uses generic elemental checkpoint records and chemical
conservation checks, including a 21-species methanol fixture
([evidence](validation/eb-carbon-checkpoint-20260917.md)). The single-level
selected frontend now shares source-derived element inventories and uses
schema-6 restart (fixed applications use schema 4). Redistributed hydro and
coupled APIs expose actual boundary, wall and chemical exchange vectors;
single-level applications accumulate and persist them, check all component
budgets and external-exchange-corrected elemental inventories. A density-contrast
methanol fixture passes exact state/budget restart
([application evidence](validation/eb-application-exchange-20260917.md),
[API evidence](validation/eb-redistributed-exchange-20260917.md),
[checkpoint evidence](validation/eb-exchange-checkpoint-20260917.md)).
A serial fixed-geometry two-level EB remapper now conservatively moves, adds
or removes fine patches, preserves same-resolution overlap, and supports
ratio changes with coarse-donor initialization and transactional failure.
Fixed/selected serial applications now connect prescribed patch schedules,
event-aligned stepping and schedule-bound checkpoint continuation; distributed
EB evolution remains open
([remap evidence](validation/eb-remap-20260917.md),
[scheduled application evidence](validation/eb-scheduled-application-20260917.md)).
An EB temperature-tagging API now selects deterministic coarse-index tiles
from exposed coarse/fine fluid temperatures, with buffered tagging and
transactional capacity rejection ([tagging evidence](validation/eb-temperature-tagging-20260917.md)).
Multipatch update APIs now accept zero fine patches and retain exact
single-grid split behavior and exchange budgets, enabling the numerical
endpoint for full derefinement ([coarse-only evidence](validation/eb-coarse-only-20260917.md)).
Fixed/selected applications now support coarse-only initialization, scheduled
2-to-0-to-2 patch transitions and exact checkpoint continuation through the
empty topology ([application evidence](validation/eb-zero-patch-application-20260917.md)).
Temperature-policy inputs now connect exposed-cell tagging, conservative
remapping and policy-bound adaptive restart at requested evaluation times.
High-temperature and no-hot-cell cases exercise both refinement and complete
derefinement ([adaptive evidence](validation/eb-adaptive-application-20260917.md)).
A core multilevel EB tree now reuses regular AMR topology for finest-lattice
geometry extraction, recursive PCM initialization, bottom-up restriction,
exposed-fluid inventories and transactional exposed-cell chemistry. Three-level
tests cover nested conservation and rollback; multilevel hydro/transport and
application integration remain open
([tree evidence](validation/eb-patch-tree-20260918.md)).
Spatial EB tree face halos now use a private averaged snapshot with same-level
peer precedence, ancestor fallback and physical outflow, including neighbors
under different parents. Multilevel flux synchronization and stepping remain
open ([halo evidence](validation/eb-tree-halo-20260918.md)).
A serial tree flux synchronization API now reconciles peer faces, restricts
fine transfers into both incident parent-level patches, corrects exposed
states and averages down in one transaction
([flux evidence](validation/eb-tree-flux-20260918.md)). Manufactured-transfer
conservation tests cover nested and one-sided refinement. A serial common-step
Euler hydro API now connects snapshot halos, all-node updates, synchronization
and separate exposed-wall/physical-boundary budgets transactionally
([hydro evidence](validation/eb-tree-hydro-20260918.md)). A matching molecular
transport API now uses limited node fluxes and globally identified closed
physical faces with transactional synchronization
([transport evidence](validation/eb-tree-transport-20260918.md)). A common-step
R/2 T/2 H T/2 R/2 API now combines these operators with separate boundary,
wall and chemical ledgers and complete rollback
([coupled evidence](validation/eb-tree-full-20260918.md)). A tree checkpoint
stream envelope now validates topology, all embedded node records and composite
budgets before transactional restoration; fixed-mechanism continuation is exact
in the focused fixture ([checkpoint evidence](validation/eb-tree-checkpoint-20260918.md)).
A serial fixed-mechanism
[tree application](../cases/reactive_eb_tree_3d/README.md) now connects static
input, coupled steps, exposed-cell/budget CSVs and staged checkpoint saving
with exact focused continuation ([application evidence](validation/eb-tree-application-20260918.md)).
Optional automatic common intervals now combine all-node cut-cell hydro and
diffusion bounds with global donor maxima and restart-bound policy
([interval evidence](validation/eb-tree-auto-interval-20260918.md)). MPI remains open for the EB tree application.
A conservative tree remap library now supports patch movement, coarsening and
refinement on an unchanged refinement lattice, with exact same-level overlap
and transactional geometry/EOS rejection
([remap evidence](validation/eb-tree-remap-20260918.md)). The application now
connects prescribed topology changes, exact event landing and schedule-bound
restart across coarsening/refinement ([schedule evidence](validation/eb-tree-scheduled-20260918.md));
a temperature-tagging library now reuses deterministic AMR tiling while
preserving exposed fine peaks and pruning solid-only nodes
([tagging evidence](validation/eb-tree-tagging-20260918.md)). The application now evaluates this policy at scheduled times, conservatively
remaps, and restores completed generated layouts under a policy-bound V4
checkpoint ([adaptive evidence](validation/eb-tree-adaptive-20260918.md)).
A selected-mechanism executable now shares the tree application, validates
composition and initial temperature range, checks chemical element conservation,
and binds bundle/composition identity on restart
([selected evidence](validation/eb-tree-selected-20260918.md)).
Tree face halos now optionally interpolate provisional ancestor endpoints
with independent per-level fractions while preserving current peer donors
([time-halo evidence](validation/eb-tree-time-halo-20260918.md)). A serial
recursive hydro API now advances ratio-based fine substeps, reconciles peers,
accumulates effective transfers and refluxes each parent interval
([subcycle hydro evidence](validation/eb-tree-subcycle-hydro-20260918.md));
a matching transport API now uses ratio-squared substeps and limited-flux
synchronization ([subcycle transport evidence](validation/eb-tree-subcycle-transport-20260918.md)).
The coupled API and application now select recursive hydro/transport inside
reaction splitting, with mode-bound scheduled restart
([subcycled application evidence](validation/eb-tree-subcycle-full-20260918.md)).
Automatic root intervals now scale patch hydro bounds by cumulative refinement
ratios and diffusion bounds by their squares
([interval evidence](validation/eb-tree-subcycle-timestep-20260918.md)).
An MPI library now restricts EB state/temperature on child owners and transfers
coarse payloads to parent owners, with replicated geometry validation and
exposed-cell composite reductions ([MPI transfer evidence](validation/mpi-eb-tree-transfer-20260918.md)).
MPI spatial face halos now fetch only requested owner-local U/T cells, using
same-level peers, ancestor PCM donors and physical outflow, with exact focused
serial parity ([halo evidence](validation/mpi-eb-tree-halo-20260918.md)).
MPI EB flux synchronization now reconciles peer faces, restricts and mirrors fine
transfers into parent peers, and refluxes exposed owner-local cells. A common
Euler hydro API connects these operations with face exchange and separate
boundary/wall ledgers, matching serial states exactly in focused cases
([MPI hydro evidence](validation/mpi-eb-tree-hydro-20260918.md)).
MPI common-interval molecular transport and exposed-cell reaction now wrap
hydro in transactional R/2 T/2 H T/2 R/2 composition, with focused serial state
parity for explicit/implicit chemistry ([MPI coupled evidence](validation/mpi-eb-tree-full-20260918.md)).
The fixed-mechanism `pelef_mpi_reactive_eb_tree_3d` application now connects
static or scheduled common-interval stepping, CSV/checkpoints and
changed-rank/serial restart ([MPI application evidence](validation/mpi-eb-tree-application-20260918.md)).
MPI automatic intervals now combine collective remote-donor maxima and owner-local
patch bounds, with exact focused serial timestep/state parity and changed-rank
restart ([automatic MPI evidence](validation/mpi-eb-tree-auto-20260918.md)).
Sparse MPI EB remapping now supports conservative patch movement, coarsening,
refinement and changed owners on the shared immutable lattice. Prescribed
application regrids and restart across their events are connected
([scheduled MPI evidence](validation/mpi-eb-tree-remap-20260918.md)).
MPI temperature adaptation now gathers exposed fluid-cell votes, preserves fine
hot spots and removes solid-only tiles. Adaptive application events and changed-rank
restart match the focused serial case
([adaptive MPI evidence](validation/mpi-eb-tree-adaptive-20260918.md)).
MPI ancestor halo donors now support independent time fractions at each depth,
with exact focused serial parity and collective rollback
([temporal MPI evidence](validation/mpi-eb-tree-time-halo-20260918.md)).
MPI recursive hydro now uses ratio-based fine intervals, provisional peer-corrected
endpoints and accumulated effective fluxes, matching focused serial U/T exactly
([MPI subcycle hydro evidence](validation/mpi-eb-tree-subcycle-hydro-20260918.md)).
MPI molecular transport now uses ratio-squared substeps and limited-flux
accumulation; the coupled API selects recursive H/T inside reaction splitting
([MPI recursive coupling evidence](validation/mpi-eb-tree-subcycle-full-20260918.md)).
The MPI application now dispatches recursive steps and mode-bound restart.
Automatic root intervals scale hydro bounds by cumulative refinement and diffusion
bounds by its square, with focused scheduled and automatic restart parity
([MPI subcycled application evidence](validation/mpi-eb-tree-subcycle-application-20260918.md)).
The selected-mechanism MPI EB executable now shares composition/element checks,
stepping controls and identity-bound checkpoints with serial. Its dedicated
runtime excludes fixed mechanism libraries, with focused H2/O2 implicit
changed-rank/serial restart parity
([selected MPI evidence](validation/mpi-eb-tree-selected-20260918.md)).
EB tree state CSV now uses owner-local MPI-IO chunks in serial row order,
with fixed/selected byte parity and preflight output protection
([distributed CSV evidence](validation/mpi-eb-tree-csv-20260918.md)).
Checkpoint saving now stages serial patch records on their owners and writes
them directly with MPI-IO; fixed/selected patch bytes and changed-rank/serial
restart match focused references
([distributed save evidence](validation/mpi-eb-tree-save-20260918.md)).
Checkpoint loading now scans record offsets on rank zero and restores each patch
only on its owner, with collective state/ledger publication and focused
changed-rank/serial parity
([distributed load evidence](validation/mpi-eb-tree-load-20260918.md)).
Initial tree fields now allocate only on patch owners, with sparse parent-section
fetches and the serial PCM kernel, removing full-tree state replication and scatter
([owner-local initialization evidence](validation/mpi-eb-tree-initialize-20260918.md)).
Root seed U/T is now constructed only on the root patch owner and released on
return from initialization, including restart topology setup
([evidence](validation/mpi-eb-tree-seed-20260918.md)). MPI checkpoints now include
a V2 patch byte-length index, allowing direct owner reads without a full root
scan. Serial readers also accept V2; legacy V1 saves retain the bounded-memory
scan fallback ([evidence](validation/mpi-eb-tree-index-20260918.md)). Geometry
remains replicated and scaling qualification remains open.
The MPI EB application can now rebalance whole patches using sampled owner-local
reaction/hydro/transport kernel time, with a configurable interval and minimum
estimated load improvement. Accepted common or recursive steps provide samples;
failed retry samples are discarded. State and temperature migrate together
without interpolation or a root gather
([evidence](validation/mpi-eb-tree-rebalance-20260918.md)). Timing-dependent owner
maps are performance state and do not constrain checkpoint restart identity.

A new regular 3D patch-tree library represents multiple nonoverlapping children
and arbitrary depth, with parent-before-child indexing and a common refinement
ratio at each depth. It reuses existing PCM/restriction kernels and provides
transactional initialization, bottom-up averaging and composite inventories
([design](design_decisions/amr-patch-tree-3d.md),
[evidence](validation/amr-patch-tree-3d-20260917.md)). This is a geometry/state
foundation. Snapshot ghost filling now exchanges peer data across parent boxes,
uses PCM ancestor fallback and follows regular-grid periodic/physical boundary
conventions ([evidence](validation/amr-tree-ghost-3d-20260917.md)). Per-level coarse
endpoint interpolation is also available for ghost donors while retaining exact
current peers ([evidence](validation/amr-tree-temporal-ghost-20260917.md)). Limited linear coarse-donor interpolation is now available for hydro
([evidence](validation/amr-tree-plm-ghost-20260917.md)) and transport
([evidence](validation/amr-tree-transport-plm-20260917.md)); coupled/application subcycle dispatch is opt-in. A single-depth RK2 hydro API now accepts coarse
time endpoints and returns integrated face fluxes without changing inactive levels
([evidence](validation/amr-tree-hydro-level-20260917.md)). Reflux now supports
completed intervals restricted to a minimum depth and returns effective fluxes
for parent accumulation ([evidence](validation/amr-tree-interval-reflux-20260917.md)). A recursive
hydro API now connects ratio-based substeps, coarse-time donors, interval reflux,
temperature recovery and a scaled coarse CFL bound
([evidence](validation/amr-tree-subcycle-hydro-20260917.md)). Recursive transport now
uses ratio-squared substeps, coarse-time donors, interval reflux and scaled
diffusive/isothermal-wall bounds
([evidence](validation/amr-tree-subcycle-transport-20260917.md)). Fixed/selected applications
now select the recursive reaction/transport/hydro composition, scaled timestep
bounds and mode-aware checkpoint continuation
([evidence](validation/amr-tree-subcycle-application-20260917.md)). MPI tree transfer APIs now
store payloads only on assigned patch owners, with collective metadata checks and
root scatter/gather ([evidence](validation/amr-tree-mpi-distribution-20260917.md)).
Same-geometry ownership changes now move conserved state and temperature directly
between owners with exact values and collective preflight rollback
([evidence](validation/amr-tree-mpi-redistribution-20260917.md)).
A cost-driven rebalance API combines owner-reported patch costs, a minimum
predicted gain and transactional migration. Fixed/selected MPI applications now
optionally measure owner-local kernels and rebalance at a configured step cadence
([evidence](validation/amr-tree-mpi-runtime-balance-20260917.md)).
Sparse cell queries now fetch only requested owner payloads in caller order, including
empty ranks and collective rollback ([evidence](validation/amr-tree-mpi-cells-20260917.md)).
Synchronous MPI halos now share serial peer/ancestor lookup, PCM/limited PLM and
physical boundary transforms ([evidence](validation/amr-tree-mpi-ghost-20260917.md)).
Final MPI tree CSV output now writes owner-local chunks directly with MPI-IO and
uses serial-order small-vector final inventory reduction
([evidence](validation/amr-tree-mpi-csv-20260917.md)). Initial inventory and the
restart initial absolute-inventory scale now use the same distributed traversal
([evidence](validation/amr-tree-mpi-initial-inventory-20260917.md)). Checkpoint
saves now validate EOS/parent averages collectively and write owner-local payloads
with serial-byte compatibility
([evidence](validation/amr-tree-mpi-checkpoint-write-20260917.md)). Restart now
reads indexed node ranges on their assigned owners and validates candidates
without constructing or scattering full root fields
([evidence](validation/amr-tree-mpi-checkpoint-read-20260917.md)). Optional V2
checkpoints now store validated relative offsets, eliminating that reader's
root payload scan while retaining V1 input support
([evidence](validation/amr-tree-indexed-checkpoint-20260917.md)).
MPI applications now optionally assign whole patches by estimated cell-step work
at initialization, regridding and restart
([evidence](validation/amr-tree-mpi-partition-20260917.md)).
MPI bottom-up restriction now sends restricted child boxes directly to parent
owners, with optional depth limits ([evidence](validation/amr-tree-mpi-restriction-20260917.md)).
MPI shared-face flux synchronization now exchanges overlapping face rectangles
and prioritizes fine-replaced values ([evidence](validation/amr-tree-mpi-flux-peers-20260917.md)).
Distributed reflux now connects parent face averaging, uncovered-cell correction,
physical boundary accounting and bottom-up state restriction, including depth-limited
intervals ([evidence](validation/amr-tree-mpi-reflux-20260917.md)).
Owner-local temperature recovery and global minimum CFL selection now reuse the
serial thermodynamic kernels ([evidence](validation/amr-tree-mpi-thermo-20260917.md)).
Synchronous MPI SSPRK2 hydro now connects these operations for owned patches,
with PCM/PLM and physical inflow/walls ([evidence](validation/amr-tree-mpi-hydro-20260917.md)).
Owner-local exposed-cell chemistry and reaction/hydro Strang splitting now support
explicit/implicit integration ([evidence](validation/amr-tree-mpi-chemistry-20260917.md)).
Synchronous MPI molecular-transport RK2 now connects viscosity, conduction, species
and pressure diffusion with physical/moving isothermal walls
([evidence](validation/amr-tree-mpi-transport-20260917.md)).
MPI transport/wall timestep selection and the full common-step reaction/transport/hydro
composition are now connected ([evidence](validation/amr-tree-mpi-split-20260917.md)).
The fixed-mechanism MPI tree executable now connects static topology, owner-local
initialization/stepping, timestep retries and serial-format final CSV/budget output
([evidence](validation/amr-tree-mpi-application-20260917.md)). MPI recursive hydro subcycling and its scaled coarse CFL bound are now connected
([evidence](validation/amr-tree-mpi-subcycle-hydro-20260917.md)). MPI recursive transport now connects ratio-squared steps, interval reflux and scaled
diffusion/wall bounds ([evidence](validation/amr-tree-mpi-subcycle-transport-20260917.md)).
MPI reaction/transport/hydro composition now accepts recursive stepping
([evidence](validation/amr-tree-mpi-subcycle-split-20260917.md)); fixed/selected MPI application mode dispatch, scaled timestep selection and mode-aware
checkpoint continuation are now connected
([evidence](validation/amr-tree-mpi-subcycle-application-20260917.md)).
Focused hydro-only scheduled/adaptive recursive cases now cover fixed/selected
serial parity and changed-rank restart
([evidence](validation/amr-tree-mpi-subcycle-dynamic-20260917.md)). Scheduled recursive full-split
execution now also combines chemistry, molecular transport, moving walls, PLM,
measured-work rebalancing and V2 changed-rank continuation with exact serial bytes
([evidence](validation/amr-tree-mpi-dynamic-split-20260917.md)). Temperature-triggered
full-split refinement now has the same combined V2/rebalance/changed-rank check,
including no repeated tag on restart
([evidence](validation/amr-tree-mpi-adaptive-split-20260917.md)). MPI transport now advances one depth with
coarse-time donors, moving/isothermal walls and integrated face fluxes
([evidence](validation/amr-tree-mpi-transport-level-20260917.md)). MPI hydro
now advances one depth with coarse-time donors and returns integrated face fluxes
([evidence](validation/amr-tree-mpi-hydro-level-20260917.md)). MPI ghost filling now
accepts coarse time endpoints and per-depth interpolation fractions, with optional
target-depth selection ([evidence](validation/amr-tree-mpi-temporal-ghost-20260917.md)).
Fixed/selected MPI applications now connect temperature-adaptation controls,
post-step tagging/remapping and adaptation-aware changed-rank checkpoint continuation
([evidence](validation/amr-tree-mpi-adaptation-20260917.md)).
A distributed temperature-tagging API now shares only hot exposed-cell boxes and
uses the serial tile builder to produce identical geometry across ranks
([evidence](validation/amr-tree-mpi-tagging-20260917.md)).
Fixed/selected MPI tree applications now connect scheduled sparse remapping, event-time
clipping, conservation checks and schedule-aware changed-rank restart
([evidence](validation/amr-tree-mpi-schedule-20260917.md)).
A sparse MPI conservative-remap API now transfers intersecting old exposed cells
to new owners, preserves serial accumulation order, and recovers temperatures
transactionally ([evidence](validation/amr-tree-mpi-remap-20260917.md)).
Fixed/selected MPI tree checkpoint/restart now reuses the serial format and validates
context/geometry/EOS/ledger before owner redistribution
([evidence](validation/amr-tree-mpi-restart-20260917.md)). The selected-mechanism
MPI tree frontend now uses an isolated selected runtime and forwards the bundle
integrator ([evidence](validation/amr-tree-selected-mpi-20260917.md)). A transactional tree
flux-correction API now reconciles shared and periodic faces, propagates nested
fine fluxes, corrects uncovered cells and returns the physical boundary-ledger
replacement ([evidence](validation/amr-tree-reflux-3d-20260917.md)). Synchronous
SSPRK2 hydro now connects PCM/PLM face fluxes, temperature recovery and the
composite boundary ledger, with a CFL minimum over all nodes
([evidence](validation/amr-tree-hydro-20260917.md)); PLM hydro now enables limited coarse donors.
Cell chemistry now runs on exposed tree cells, with bottom-up restriction and
transactional reaction/hydro Strang splitting
([evidence](validation/amr-tree-chemistry-20260917.md)). Synchronous molecular
transport now includes viscous/thermal/species fluxes, nested reflux and physical
wall budgets ([evidence](validation/amr-tree-transport-20260917.md)). A combined
R/2-T/2-H-T/2-R/2 update now aggregates those boundary budgets transactionally
([evidence](validation/amr-tree-split-20260917.md)). A conservative overlap-remap
API now supports changed patch layouts and refinement ratios, with transactional
EOS recovery ([evidence](validation/amr-tree-remap-20260917.md)). Automatic solution-based tagging and MPI ownership are now connected as described above. The experimental
[serial tree application](../cases/amr_tree_reactive_3d/README.md) now reads explicit
patch topology, initializes analytic fields, evolves the split system and writes
exposed-cell CSV plus composite budgets
([evidence](validation/amr-tree-application-20260917.md)). Prescribed-time topology
changes now connect conservative remapping, exact event-time stepping and budget
continuity ([evidence](validation/amr-tree-regrid-20260917.md)). Configure-time
selected mechanisms now use the same tree driver through an isolated runtime
([evidence](validation/amr-tree-selected-20260917.md)). Versioned tree snapshot I/O
and physics-context comparison APIs are now available
([evidence](validation/amr-tree-checkpoint-api-20260917.md)). Fixed/selected application
stop/restart now validates schedule position, geometry and continuation context
([evidence](validation/amr-tree-restart-20260917.md)). A temperature-tagged topology
builder now creates multiple tiled children over requested depths, including
periodic buffers and fine-only hot-cell preservation
([evidence](validation/amr-tree-tagging-20260917.md)). Successful-step adaptation
cadence and fixed/selected automatic-geometry restart are now connected
([evidence](validation/amr-tree-adaptation-20260917.md)).

Sparse MPI AMR has physical ghost layers and two-level PCM/PLM hydro supporting
one coarse column per rank, with serial equality and collective rollback checks
([ghosts](validation/amr-mpi-boundary-halo-20260917.md),
[coarse hydro](validation/amr-mpi-boundary-hydro-20260917.md),
[hierarchy](validation/amr-mpi-physical-hierarchy-20260917.md)). The
[reaction/hydro split](validation/amr-mpi-physical-reaction-20260917.md) also forwards
physical boundary data and budgets. Physical molecular transport and the full
reaction/transport/hydro split are also connected
([evidence](validation/amr-mpi-physical-transport-20260917.md)). MPI application
inputs, wall/inlet timestep limits and rank-neutral boundary-ledger restart are
also connected ([evidence](validation/amr-mpi-physical-app-20260917.md)). Fixed elementary
and full-H2/O2 reacting AMR checkpoint/restart now retains canonical reaction,
composition and integrator context for serial/MPI applications
([evidence](validation/amr-fixed-chemistry-restart-20260917.md)).
The selected AMR runtimes remain independent of fixed-mechanism databases and now
include distributed checkpoint I/O in their build
([build isolation](validation/selected-amr-checkpoint-isolation-20260917.md)).

Boundary-touching AMR now has a physical flux-correction API that replaces
covered domain-surface fluxes in the boundary ledger and corrects only uncovered
interior coarse cells ([evidence](validation/amr-boundary-reflux-20260917.md)).
Serial PCM/PLM hydro now connects fine physical ghosts, subcycling, boundary
ledger replacement and average-down for such patches
([evidence](validation/amr-boundary-touching-hydro-20260917.md)). Cell-local reaction,
reaction/hydro splitting and the state-based CFL helper now accept these patches
([evidence](validation/amr-boundary-touching-reaction-20260917.md)). Molecular
transport and the full reaction/transport/hydro split are connected, including
fine wall heat/work fluxes ([evidence](validation/amr-boundary-touching-transport-20260917.md)).
Serial static application inputs, fine wall timestep bounds and checkpoint/restart
are connected ([evidence](validation/amr-boundary-touching-app-20260917.md)).
Serial temperature-driven regridding now supports boundary-touching patches with
explicit physical tagging and geometry-aware restart
([evidence](validation/amr-boundary-regrid-20260917.md)). Sparse MPI distribution,
scatter/gather and cell-local reaction accept boundary-touching geometry
([evidence](validation/amr-mpi-boundary-distribution-20260917.md)). MPI PCM/PLM
hydro and reaction/hydro splitting now connect fine physical ghosts and covered
surface flux replacement ([evidence](validation/amr-mpi-boundary-touching-hydro-20260917.md)).
MPI molecular transport and the full reaction/transport/hydro split now also
connect fine wall heat/work and physical flux replacement
([evidence](validation/amr-mpi-boundary-touching-transport-20260917.md)).
MPI static application inputs now include fine wall timestep limits and
rank-neutral restart for boundary-touching geometry
([evidence](validation/amr-mpi-boundary-touching-app-20260917.md)).
MPI physical-boundary temperature tagging, conservative remap and changed-rank
restart are also connected
([evidence](validation/amr-mpi-boundary-regrid-20260917.md)). Temperature tagging
now scans local slabs and reduces only tag bounds; unchanged patches skip full
hierarchy gathering ([evidence](validation/amr-mpi-distributed-tagging-20260917.md)).
Conservative remapping now runs within fixed coarse slabs using neighboring
coarse halos and local fine overlap reuse
([evidence](validation/amr-mpi-distributed-remap-20260917.md)). It no longer
gathers full fields during regridding. An explicit repartition API now transfers
coarse/fine states and temperatures directly between ranks, and remapping
preserves the chosen coarse partition
([evidence](validation/amr-mpi-repartition-20260917.md)). Optional cell-count
balancing is connected after initialization/restart and geometry changes
([evidence](validation/amr-mpi-load-balance-20260917.md)). General patches/levels,
measured performance/scaling and broader physical qualification remain open. MPI AMR
CSV output now uses direct collective file writes, and final inventory, synchronization
and composition diagnostics use distributed reductions
([evidence](validation/amr-mpi-distributed-csv-20260917.md)). Analytic initial
conditions are also generated on local coarse/fine slabs, with local average-down
and distributed initial inventories
([evidence](validation/amr-mpi-local-initialization-20260917.md)). Checkpoint writes
now serialize the existing format directly from local slabs after collective
metadata, temperature and synchronization checks
([evidence](validation/amr-mpi-checkpoint-output-20260917.md)). Restart now reads
canonical checkpoint cell rows directly into local candidate slabs, then commits
geometry, fields and diagnostics only after collective checks
([evidence](validation/amr-mpi-checkpoint-input-20260917.md)). Metadata remains
replicated; reformatted cell records use a bounded-memory streaming fallback.
Production I/O scaling is not yet measured.

Ahamed wall-impact classification, splash velocities and azimuthal mass
fractions, four-child construction and deposited cylindrical film patches are
connected to the regular coupled-step API with wall budgets, stationary film
heating/evaporation with per-material constant or cubic thermal conductivity,
CLI inputs/output and coupled checkpoint persistence
[evidence](validation/wall-film-conductivity-20260916.md). Lateral film transport
and physical impact qualification remain open. The gas-side film skin now uses
the pinned half-temperature/half-composition rule
([correction evidence](validation/wall-film-skin-20260916.md)).

KH shed-mass accumulation and four-child release now connect to regular CFD
and restart ([evidence](validation/khrt-children-20260916.md)); computational
splitting now partitions the KH inventories
([evidence](validation/khrt-splitting-20260916.md)). Stationary wall impact
now consumes crossed parcels before KH release, with coupled film/lineage restart
([evidence](validation/khrt-wall-20260916.md)); wet-impact correlation limits remain.

The mixed-liquid API now covers mixture properties, component transport,
regular-grid feedback, component-specific boundary budgets, representative
breakup/splitting and composition checkpoint persistence
([evidence](validation/spray-mixture-grid-20260916.md)). The combined
hydro/transport/LES/chemistry driver now retains composition and rolls it back
with all other state ([evidence](validation/spray-mixture-reactive-20260916.md)).
Mixture CLI inputs, shared/independent scheduled injection, component CSVs and
restart compatibility are connected
([evidence](validation/spray-mixture-cli-20260916.md)). Trace-species depletion and
minimum-diameter removal now return their complete inventories to gas
([evidence](validation/spray-mixture-terminal-20260916.md)). Mixed KH children now inherit composition with conservative transverse velocities
and retain their state through splitting and restart
([evidence](validation/spray-mixture-kh-20260916.md)). The mixture impact event API now partitions component masses into splash and film
with wall budgets ([evidence](validation/spray-mixture-wall-event-20260916.md));
its CLI connection is now available. Planar crossing now uses mixture properties
and the reference mass-weighted impact boiling temperature
([evidence](validation/spray-mixture-wall-crossing-20260916.md)). A supplied film-exchange API now
updates composition, temperature and height with gas/wall energy closure
([evidence](validation/spray-mixture-film-exchange-20260916.md)). Mixed film
heat/evaporation rates and a frozen-rate update are now available
([evidence](validation/spray-mixture-film-rates-20260916.md)). The checkpoint API
now retains film component masses with transactional validation
([evidence](validation/spray-mixture-film-checkpoint-20260916.md)); application wiring is connected. A transactional mixed-film grid feedback API now
conserves species and gas/film/wall energy and removes terminal patches
([evidence](validation/spray-mixture-film-grid-20260916.md)). The mixed particle
coupler now connects wall crossing, splash/deposition and remaining-interval film
evolution in one transaction ([evidence](validation/spray-mixture-wall-coupled-20260916.md)).
The reactive driver now retains mixed film state across both spray half-steps and
rolls it back after late failure ([evidence](validation/spray-mixture-wall-reactive-20260916.md)).
CLI wall inputs, film component output, conservation totals and restart are connected
([evidence](validation/spray-mixture-wall-cli-20260916.md)). Active trace-film
components now transfer their full residual inventories to gas
([evidence](validation/spray-mixture-film-trace-20260916.md)). Per-port compositions now drive initial and
scheduled injection with component budgets and restart identity
([evidence](validation/spray-port-mixture-20260916.md)). Component density/heat-capacity fits, Antoine
pressure and Clausius surface caps are available in the CLI and restart context
([evidence](validation/spray-mixture-properties-20260916.md)).

Conservation, rollback, byte-parity, and short-time H2/O2 tests do not establish
physical ignition/spray validation or complete PeleC equivalence. Choosing a
license remains an owner decision.

AMR particle ownership now locates each parcel in the finest covering regular
3D tree cell and constructs stable owner-rank bins while preserving all parcel
fields and optional component inventories. Periodic lookup, half-open patch
faces, regrid reassignment and transactional validation are supported
([evidence](validation/spray-amr-partition-20260918.md)). MPI migration now uses
this ownership map and a shared Cartesian/AMR particle exchange, preserving
component and breakup payloads with collective validation/rollback through
particle crossings, owner changes and regridding
([evidence](validation/spray-amr-migration-20260918.md)). Serial composite-field
sampling/deposition now uses a shared positive cloud-overlap stencil over exposed
cells, with volume-scaled deposition and covered-parent averaging
([evidence](validation/spray-amr-coupling-20260918.md)). A serial gas-feedback
wrapper now recovers EOS temperatures after deposition and parent averaging,
committing gas U/T only when every cell accepts the update
([evidence](validation/spray-amr-gas-20260918.md)). MPI particle-position sampling
now requests only the composite stencil cells from their owners and matches
serial accumulation order ([evidence](validation/spray-amr-sampling-20260918.md)).
MPI field deposition now routes volume-scaled particle contributions to cell
owners and refreshes covered parents without a full-field gather
([evidence](validation/spray-amr-deposition-20260918.md)). MPI gas feedback now
recovers owner-local temperatures after parent averaging and commits U/T
collectively, with initial-state validation and failure rollback
([evidence](validation/spray-amr-mpi-gas-20260918.md)). A periodic MPI AMR
mixed-particle interval now connects drag/heating/evaporation, representative
breakup and KH children, common substeps, transactional gas feedback and final
finest-owner migration ([evidence](validation/spray-amr-step-20260918.md)).
An MPI reacting-spray step now composes particle half-intervals around the
existing reaction/transport/hydro split, with common or recursive gas steps and
whole-step rollback ([evidence](validation/spray-amr-reactive-20260918.md)).
A coupled remap API now commits tree/owners, conservative gas remapping and
finest-owner particle migration together, preserving all original state on
failure ([evidence](validation/spray-amr-remap-20260918.md)).
Fixed and selected [MPI AMR spray applications](../cases/mpi_amr_spray_3d/README.md)
now connect owner-local gas setup, mixed-cone initial injection, coupled stepping,
one scheduled regrid and gas/particle CSV output
([evidence](validation/spray-amr-application-20260918.md)). Explicit single-cone
injection schedules now land on event times and commit particle payloads, IDs and
cumulative injection inventories after successful migration
([evidence](validation/spray-amr-events-20260918.md)). AMR spray checkpoints now
bind gas/particle state, topology, IDs and event position, with immutable
publication and changed-rank continuation in fixed/selected executables
([evidence](validation/spray-amr-checkpoint-20260918.md)). Up to 16 AMR spray
ports now support per-port geometry, mixture fractions, nozzle profile/swirl and
equal-time injection events with preserved port IDs across restart
([evidence](validation/spray-amr-ports-20260918.md)). Temperature-driven AMR
adaptation now updates gas/particle topology on completed-step intervals and
continues after changed-rank restart without repeating the saved step's check
([evidence](validation/spray-amr-adaptation-20260918.md)). Gas-work sampling now
drives owner reassignment with transactional gas/particle migration and
rank-neutral continuation ([evidence](validation/spray-amr-balance-20260918.md)).
The AMR particle-interval API now handles outflow and specular reflection with
a transactional global mass/momentum/energy/species boundary ledger
([evidence](validation/spray-amr-boundaries-20260918.md)). The reacting wrapper
now combines gas and particle physical-boundary inventories transactionally
for common or recursive gas steps
([evidence](validation/spray-amr-reactive-boundaries-20260918.md)). Fixed/selected
applications now connect physical gas/particle boundaries, reservoir timestep
bounds, boundary-aware migration/regridding and cumulative boundary inventory
through changed-rank checkpoint continuation
([evidence](validation/spray-amr-open-application-20260918.md)). A stationary
wall-film migration API now preserves mixed-film payloads across finest-cell
ownership changes and regridding, including exact upper-wall positions
([evidence](validation/spray-amr-film-migration-20260918.md)). A planar mixed-film
interval API now couples heating/evaporation to owner-local AMR gas with a
transactional wall-energy ledger
([evidence](validation/spray-amr-film-step-20260918.md)). The reacting wrapper
now surrounds particle/gas stepping with film half-intervals and publishes
film state and wall heat only with whole-step success
([evidence](validation/spray-amr-film-reactive-20260918.md)). Coupled topology
remapping now carries stationary films transactionally with gas and particles
([evidence](validation/spray-amr-film-remap-20260918.md)). Gas-work-based
rebalancing also carries film state atomically with gas and particles
([evidence](validation/spray-amr-film-balance-20260918.md)). AMR particle intervals
now connect the existing planar impact model to film creation, wall impulse/energy
and optional splash children, and the reacting wrapper carries those films
([evidence](validation/spray-amr-impact-20260918.md)). Film-aware checkpoint
records now restore gas/particles/films and cumulative budgets together across
rank-count changes ([evidence](validation/spray-amr-film-checkpoint-20260918.md)).
Fixed/selected applications now expose opt-in wall impact/film settings,
film CSVs, coupled adaptation/balancing and changed-rank continuation
([evidence](validation/spray-amr-film-application-20260918.md)). A serial particle
segment query now finds first contact with the original nodal level set using
the gas EB builder's tetrahedral subdivision
([evidence](validation/spray-eb-intersection-20260918.md)). A force-free serial
population drift now reflects specular EB contacts, follows the remaining
flight time and publishes particle/wall inventories transactionally
([evidence](validation/spray-eb-drift-20260918.md)). Serial single-grid EB field
sampling and integrated gas deposition now use connected fluid-volume-weighted
stencils and transactional EOS recovery
([evidence](validation/spray-eb-field-20260918.md)). A serial EB mixture interval
now connects drag/heating/evaporation, representative-drop breakup, specular
drift and conservative gas feedback with whole-interval rollback
([evidence](validation/spray-eb-step-20260918.md)). Particle half-intervals now
wrap the existing serial reactive EB gas update, preserving separate exterior,
wall and chemistry ledgers and atomic whole-step publication
([evidence](validation/spray-eb-reactive-20260918.md)). Fixed/selected
[serial EB spray applications](../cases/eb_spray_3d/README.md) now connect initial
injection, shared nodal geometry, coupled timestep control/retries, gas/particle
CSVs and cumulative inventories
([evidence](validation/spray-eb-application-20260918.md)). Ordered timed cone
events now support tied times, zero dose, per-event mass scaling and an injected
inventory ledger ([evidence](validation/spray-eb-schedule-20260918.md)). Immutable
serial EB spray snapshots now resume gas, particles, scheduled injection and
cumulative inventories with numerical/geometry context checks
([evidence](validation/spray-eb-restart-20260918.md)). A bounded MPI ownership
primitive now migrates AMR particles with collective nodal EB fluid checks
([evidence](validation/spray-mpi-eb-tree-20260918.md)). Composite EB particle
stencils now combine exposed coarse/fine support, physical fluid volumes and
centroid visibility checks
([evidence](validation/spray-eb-tree-stencil-20260918.md)). MPI sampling now
uses these stencils to fetch only requested owner-local field values, with
collective geometry checks and transactional sample publication
([evidence](validation/spray-mpi-eb-sample-20260918.md)). MPI EB deposition now
routes integrated particle feedback by owner, restricts hidden parents with
fluid-volume weights and recovers active-cell EOS temperatures transactionally
([evidence](validation/spray-mpi-eb-deposit-20260918.md)). A distributed EB
particle interval now connects mixed-liquid source evolution, specular drift,
conservative gas feedback and final owner migration with whole-interval rollback
([evidence](validation/spray-mpi-eb-step-20260918.md)). A distributed EB gas/spray
split now combines particle half-steps with MPI reaction/transport/hydro and
commits state and global ledgers atomically
([evidence](validation/spray-mpi-eb-reactive-20260918.md)). Fixed/selected MPI EB
spray executables now connect initial mixed-liquid injection, coupled stepping,
combined budgets and per-rank particle output
([evidence](validation/spray-mpi-eb-application-20260918.md)). Scheduled injection
now supports tied, zero-scale and final-time events with an explicit injection
budget ([evidence](validation/spray-mpi-eb-schedule-20260918.md)).
A particle checkpoint payload now preserves complete parcel/component and
injection metadata across MPI rank-count changes, with rollback on rejection
([evidence](validation/spray-mpi-eb-particle-checkpoint-20260918.md)).
Combined gas/particle application checkpoints now restore on a changed MPI
rank count with pending injection events and cumulative ledgers preserved
([evidence](validation/spray-mpi-eb-restart-20260918.md)).
The application also redistributes whole patches using measured gas work and
commits gas/particle ownership together
([evidence](validation/spray-mpi-eb-balance-20260918.md)).
A coupled remap library now moves refined patches, coarsens and refines with
transactional gas/geometry/particle publication and retained original nodal
geometry ([evidence](validation/spray-mpi-eb-remap-20260918.md)).
Prescribed EB spray regrids now connect to the application and changed-rank
checkpoint continuation across stored hierarchy stages
([evidence](validation/spray-mpi-eb-regrid-application-20260918.md)).
MPI EB injection now supports independent port geometry, count, mass, diameter,
temperature and composition selected per event, with restart-bound inputs
([evidence](validation/spray-mpi-eb-ports-20260918.md)). Per-port tabulated rates
now integrate explicit packet intervals with restart-bound profile inputs
([evidence](validation/spray-mpi-eb-profile-20260918.md)). External curve-file
loading and physical qualification remain open. Liquid wall-impact support is described in the subsequent evidence. Lateral film transport remains a future extension beyond the pinned stationary-film path ([source audit](validation/film-reference-scope-20260918.md)).
Coarse/fine interpolation currently uses piecewise-constant
cell values; higher-order interface treatment/convergence remain unqualified.

Particle MPI groundwork now includes Cartesian cell ownership and stable send-buffer
partitioning with component/KH state preservation
([evidence](validation/spray-partition-20260916.md)). MPI migration and receive-side ID ordering now pass 2/4-rank checks
([evidence](validation/spray-mpi-migration-20260916.md)). Distributed gas/film
feedback and application integration remain open. Collective 64-bit ID reservation
now rejects overflow and inconsistent rank histories
([evidence](validation/spray-mpi-ids-20260916.md)). Collective computational splitting
now uses these IDs while preserving component and KH inventories
([evidence](validation/spray-mpi-split-20260916.md)). Collective cone injection now
reserves IDs and retains mixed component budgets with rollback
([evidence](validation/spray-mpi-injection-20260916.md)). Distributed field sampling
now requests only interpolation cells from owner ranks and matches the serial
stencil ([evidence](validation/spray-mpi-sampling-20260916.md)). Owner-routed
feedback deposition now matches serial stencil accumulation
([evidence](validation/spray-mpi-deposition-20260916.md)). Gas feedback now recovers
owner-block temperatures and commits conserved state and temperature collectively
([evidence](validation/spray-mpi-gas-20260916.md)). Periodic mixed-particle drag/heating/evaporation now uses owner sampling, a
common substep, collective feedback retry and end-of-interval migration with
whole-interval rollback ([evidence](validation/spray-mpi-coupling-20260916.md)).
Representative ETAB and KH-RT breakup now run inside the MPI mixed-particle
candidate update with component-inventory preservation
([evidence](validation/spray-mpi-breakup-20260917.md)). Breakup CLI settings and
deformation/RT checkpoint state now support changed-rank continuation
([evidence](validation/spray-mpi-breakup-app-20260917.md)). A collective KH shedding
API now constructs children with unique IDs, inherited composition and a gas-energy
ledger ([evidence](validation/spray-mpi-kh-release-20260917.md)). MPI coupled
feedback, input and shedding/lineage restart are now connected, including ID
rollback on rejected updates
([evidence](validation/spray-mpi-kh-coupled-20260917.md)). Periodic x-slab SGS now
shares the serial Smagorinsky face kernel, with MPI SSPRK2, coupled half-steps,
CLI coefficients and restart identity
([evidence](validation/spray-mpi-les-20260917.md)). Stationary mixed wall-film
patches can now migrate to Cartesian owners while retaining geometry and component
inventories ([evidence](validation/spray-mpi-film-migration-20260917.md)). An MPI
film-coupling API now connects heating/evaporation to owner-routed gas feedback
with a global wall-energy ledger and whole-interval rollback
([evidence](validation/spray-mpi-film-coupled-20260917.md)). MPI checkpoints now
collect and redistribute film/component state and the wall budget across rank
changes ([evidence](validation/spray-mpi-film-restart-20260917.md)). A collective
mixed-wall impact event API now reserves splash IDs and commits parent, child,
film and wall-budget updates transactionally
([evidence](validation/spray-mpi-impact-20260917.md)). Distributed film-footprint
queries now sum overlapping film heights across owners
([evidence](validation/spray-mpi-film-query-20260917.md)). An MPI planar-crossing
wrapper now computes mixture impact conditions, queries film thickness and resolves
wall events ([evidence](validation/spray-mpi-wall-crossing-20260917.md)); integration
with nonperiodic gas/particle stepping and the application remains open.
A two-layer MPI ghost API now combines physical wall/inlet/outlet conditions with
neighbor exchange, including transverse corners
([evidence](validation/spray-mpi-boundary-ghost-20260917.md)). MPI hydro now uses
these boundaries in both SSPRK2 stages and returns a global physical-face flux
ledger ([evidence](validation/spray-mpi-boundary-hydro-20260917.md)).
MPI molecular transport now shares the serial physical-boundary flux kernel,
including isothermal conduction and tangential wall work, with a global boundary
ledger ([evidence](validation/spray-mpi-boundary-transport-20260917.md)).
MPI SGS now uses the same physical-boundary ghost and wall-heat rules and
returns a global surface ledger ([evidence](validation/spray-mpi-boundary-les-20260917.md)).
The MPI particle-coupling API now supports open-face removal, stationary wall
reflection and nonperiodic owner migration with a replicated boundary ledger
([evidence](validation/spray-mpi-particle-boundary-20260917.md)). Mixed wall-impact/film
creation now participates in the coupled driver transaction
([evidence](validation/spray-mpi-impact-coupled-20260917.md)). Application impact/film
inputs, per-rank film outputs, total inventories and splash-lineage restart now
connect these models ([evidence](validation/spray-mpi-film-app-20260917.md));
physical qualification remains open.
Physical boundary data now pass through the complete MPI particle/chemistry/
transport/SGS/hydro split, with a summed boundary ledger and whole-step rollback
([evidence](validation/spray-mpi-boundary-split-20260917.md)). Application inputs now
expose physical faces and wall data; history and rank-neutral restart retain the
cumulative boundary budget ([evidence](validation/spray-mpi-boundary-app-20260917.md)).
Existing distributed films can now heat/evaporate in symmetric half-steps around
the complete coupled update, with film/gas rollback and wall-budget accumulation
([evidence](validation/spray-mpi-film-split-20260917.md)).
Periodic particle/reaction splitting now connects cell-local chemistry between
particle half-steps with collective whole-interval rollback
([evidence](validation/spray-mpi-reaction-20260916.md)). Distributed wall/breakup
events remain open. Periodic
x-slab gas advection now exchanges two halo layers and uses characteristic PLM /
SSPRK2, including optional integration between reaction half-steps
([evidence](validation/spray-mpi-hydro-20260917.md)). Molecular transport now
exchanges gas and species-limiter halos and connects as half-steps around hydro
([evidence](validation/spray-mpi-transport-20260917.md)). The fixed-mechanism
`pelef_mpi_spray_3d` executable now connects initial two-component cone injection,
global time-step selection/retry, coupled stepping, budget history and rank-local
CSV output ([evidence](validation/spray-mpi-application-20260917.md)). Scheduled
independent multiport injection now merges up to 16 port schedules with deterministic
simultaneous-event order and changed-rank restart
([evidence](validation/spray-mpi-multiport-20260917.md)). Independent per-port rate
tables now override the shared profile with complete restart identity
([evidence](validation/spray-mpi-port-profiles-20260917.md)). Finite-radius nozzle,
solid/hollow/swirl and size-distribution controls now reach the MPI injector
([evidence](validation/spray-mpi-nozzle-20260917.md)). Shared piecewise-linear injection
profiles now provide per-event integrated mass, including zero-dose event handling
and restart identity ([evidence](validation/spray-mpi-profile-20260917.md)). A shared MPI application
now also serves `pelef_mpi_spray_3d_selected`, with 21-species FFCM1_Red execution
and restart evidence ([validation](validation/spray-selected-mpi-20260917.md)).
The selected case also exposed and corrected roundoff-negative final-species mass
closure ([correction](validation/reactive-zero-species-closure-20260917.md)). Repeated timed cone events now clip time steps, retain injected budgets and
resume from the saved event count ([evidence](validation/spray-mpi-schedule-20260917.md)).
MPI spray CLI save/restart
now preserves model identity and budgets, with focused same-rank and 2↔4-rank
continuation evidence ([validation](validation/spray-mpi-restart-application-20260917.md)).
The checkpoint API now collects mixed particles and owner gas blocks into a root
file and redistributes them with a different communicator size
([evidence](validation/spray-mpi-checkpoint-20260917.md)).
The checkpoint API now preserves a 64-bit ID high-water mark with strict validation,
ready for MPI reservation state ([evidence](validation/spray-id-checkpoint-20260917.md)).

The spray CVODE branch now builds after removal of misplaced wall-film setup;
focused native/CVODE chemistry tests pass
([repair evidence](validation/spray-cvode-build-repair-20260916.md)).

## Repository consolidation implemented

The README is an entry point rather than a chronological implementation log.
Contributor instructions, AGENTS.md, the documentation index and PR template
separate current capability, historical evidence and future work.

The 2613-line root CMake and 9138-line test registration were split into five
root target modules and nineteen test modules, preserving include order and
directory scope. Generated target graphs, install rules and CTest commands
matched in seven clean configurations; see [the layout record](validation/build-layout-20260909.md)
and [CMake navigation](../cmake/README.md).

Hosted workflows cover source provenance/regeneration, complete serial and MPI
suites, selected 0D CVODE with pinned SUNDIALS 7.2.0 and Cantera 3.2.0, and a
separate tests-disabled Release installation. Source identity, dependency
versions, JUnit and logs are retained. An implemented workflow is not itself a
passed qualification; completed results must be associated with its tested SHA.

Historical draft #11 was reviewed and closed as superseded: its 17 numerical,
mechanism, application/case and verification paths were byte-identical at the
later integrated full-chemistry commit. Its unique decision note is preserved
[in the historical record](history/pr11-implicit-chemistry-decision.md).
No historical branch was deleted. The audited ancestor-only bulk cleanup was
blocked by the connection's safety restriction, not performed silently.

## Evidence ledger

The saved checkpoint is `509d07929d52d9e570a6af1f6ec7e5836bd4fe5d`.
Initial PR CI failed because its vendored YAML description lacked two spaces.
Restoring the exact Cantera package bytes repaired source provenance without
changing the mechanism bundle or generated Fortran; see [the integrity repair](validation/checkpoint-integrity-20260908.md).

A later targeted numerical review found a false-success return in NASA7 molar
conversion overflow guards. Commit `bb2980e52c2e6a4bfed353dd280fe1f0e21f4893`
corrects the failure flag without changing physical formulas or tolerances.
New cp/enthalpy/entropy rejection cases fail on the parent implementation and
pass after the fix; seven focused GNU14 Debug tests passed. See [red/green evidence](validation/nasa7-conversion-status-20260909.md).
This correction is intentionally separate from the behavior-preserving layout.

**Reported historical local qualification:** the [2026-09-08 handoff](handoff-pelef-0.245.0-20260908.md)
reports 4570/4570 across Debug 524, Release 524, Debug+Cantera 529, MPI Debug 694,
MPI Release 694, CVODE Debug 533, CVODE Debug+Cantera 539, and CVODE Release 533.
It also reports tests-disabled 42/42 ELF, separate tests-enabled 23/23 installed
restart checks, and installed application/checker/byte comparisons of 32/32,
8/8 and 12/12. These belong to that frozen worktree and its local evidence;
the raw server directories were not retrieved here. Interrupted 227-test and
excluded literal-prefix exploratory results are not successful qualification.

**Independent source-integrity qualification:** commit
`566010ffc77a6de25d452d084ff61fa2bb533d95` passed
[run 34208154353](https://github.com/GTT10/Fortranslate/actions/runs/34208154353),
including 14 checker regressions, exact hashing and full bundle regeneration.
Final serial/MPI/CVODE/clean-install results belong to PR #102 and their own
source identities. They must not be replaced by older source successes.

The earlier focused 10/10 serial and 23/23 MPI results remain in
[validation/0.245.0.md](validation/0.245.0.md). Historical evidence is preserved,
not overwritten or added together as though all configurations were identical.

See [regular 3D boundary inputs](../cases/reactive_boundaries_3d/README.md) for the new
hydro/chemistry boundary scope and flux accounting.

## Remaining work

[The roadmap](completion_roadmap.md) distinguishes further solver development
from repository consolidation. General 3D AMR/EB, production detailed-fuel stiff CFD,
full LES/spray qualification and external physical validation remain genuine engineering work,
not features implicitly completed by housekeeping or large test counts.
