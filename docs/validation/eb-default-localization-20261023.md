# Frozen default-policy reference discrepancy localization, 2026-10-23

The saved default TGReact EB-AMR comparison first violates field limits at
200 ns, after passing the 100-ns comparisons. This audit reprocesses saved
native and PeleF tables; it is not a rerun of the current source. Reference,
input and limits remain those of `cases/tgreact_eb_amr_reference/acceptance.json`.
The full 300-ns frozen trajectory remains a failed fixed-reference comparison.

Run:

```
python3 build/eb-default-localization-20261023/analyze.py
```

The script checks equality of saved/current acceptance contracts, finite
values, times, unique 1024-cell hierarchy, planar volume fractions and matching
cell keys. Native CGS values are converted to SI. It recomputes all 68 maxima
(17 fields at four times) exactly as recorded in the historical checker.
`build/eb-default-localization-20261023/localization.json` stores hashes of
the eight input tables and checks.json, every recomputed maximum and location,
and volume-weighted squared-error fractions. The audit passes; this does not
change the failed trajectory verdict.

At the first failing output, 200 ns, the maxima are:

| Field | Maximum difference | Fine cell (i,j,k) | Fraction of weighted squared error in fine k=8 layer |
| --- | ---: | --- | ---: |
| Density | 3.893796188259024e-9 kg/m3 | (8,5,8) | 45.02% |
| x momentum | 3.332689058943572e-6 kg/m2/s | (5,8,8) | 37.71% |
| y momentum | 3.2811818406752735e-6 kg/m2/s | (8,4,8) | 37.49% |
| z momentum | 1.1985395956997172e-6 kg/m2/s | (3,7,8) | 61.97% |
| Total energy density | 0.013450485916109756 J/m3 | (1,5,8) | 44.33% |

Held-energy maximum also occurs at fine (1,5,8). All these maxima have
volume fraction one. The fine upper z boundary is an internal coarse/fine
interface, not the physical upper wall. At 300 ns the x/y momentum maxima
move into fine k=6, so boundary-localized maxima at onset do not imply that
all subsequent error stays on the interface. Temperature and species limits
are not violated in the saved full trajectory.

The next discriminator is the default-policy second-step coarse/fine
boundary state/flux history, especially the upper fine z interface. Check
source/binary identity before treating a new trace as the same trajectory.
This localization alone does not distinguish ancestor time interpolation,
spatial interpolation, EOS differences or synchronization. The continuous
quartic shear candidate is not adopted: its separate reacting trial already
exceeds field limits at 100 ns and timed out before the endpoint.

## Reproduced default trajectory and reacting-stage audit

The isolated frozen diagnostic executable reproduced all eight historical
ordinary CSVs byte for byte in 49.16772105358541 s. Receipt:
`build/eb-default-identity-20261023/record.json`. This confirms the executable
used for the new trace preserves the saved default trajectory; it is not
an audit of every current production module or installed executable.

A probe-only extension then records stage and synchronization observers at
the first-failure locations. The trace run completed in 50.07086436264217 s,
exit zero, again with all eight ordinary CSVs byte-identical. Build commands,
frozen probe source and source/binary hashes are retained under
`build/eb-default-stage-20261023`; binary SHA256 is
`9375e450008e788a4e8d543ad73c399b81408ad04c3089913c3759837448ce3e`.

`python3 build/eb-default-stage-20261023/audit.py` passes 216 six-face spatial
balance groups (1296 rows) and 135 pre/post synchronization pairs (1620 rows).
It checks finite data, directional momentum divergence plus wall source,
raw register mismatch, and applied state increments with the existing
1e-12 scaled tolerance. It rechecks the eight ordinary output hashes.
The stage sequence is 0/1/2, not the hydro-only two-stage shear sequence;
reaction/transport prevent importing the previous full-trajectory formula.

At step 2, the sum of measured direct closure increments at probe 13
(fine 5,8,8; x-momentum maximum) is +7.435210084283383e-7 in x momentum.
At probe 14 (fine 8,4,8; y-momentum maximum) it is -7.342264856102787e-7
in y momentum. These are nonzero, unlike the isolated shear candidate's
maximum probes. Their other component increments can be much larger:
probe 13 y increment is -7.687031474756106e-5. This is physical/numerical
state evolution within PeleF, not a measured difference from native PeleC.
Neither subtracting it from the endpoint reference error nor calling it
an erroneous correction is justified without the matched native increments.

The next evidence needed is a matching native stage/closure trace, or an
independent replay of its coarse/fine interpolation at these clocks. Preserve
the original native capture and verify any diagnostic instrumentation leaves
its ordinary outputs unchanged. Local balance success alone does not resolve
C2's reference discrepancy. No numerical policy or tolerance was changed.

## Matched native root-reflux increments

A serial read-only trace was inserted immediately before/after the pinned
PeleC `FluxRegister::Reflux` call in an isolated copy of `Source/PeleC.cpp`.
The original source hash matches the pinned manifest; only this object was
rebuilt, with every reused link object/archive hashed. The original executable
hash matches the captured reference. Original sources and capture are unchanged.
Artifacts: `build/pelec-default-sync-20261023/{build-record,record}.json`.
The diagnostic native run completed in 5.2808031011372805 seconds, exit zero.
Exported 1024-cell tables at all four times match freshly exported captured
reference tables exactly, including file hashes. The 30 trace rows cover
five full-volume fine cells, two phases, three root steps.

`analyze.py` compares 36 component increments at four matched probe locations
and three clocks, converting native CGS momentum to SI (factor 10). It pairs
native root reflux with PeleF depth-zero closure, not intermediate fine-level
closures. At the first failing output (200 ns):

| Location/component | Native increment | PeleF increment | Increment difference | Difference before closure | Difference after closure |
| --- | ---: | ---: | ---: | ---: | ---: |
| fine (5,8,8), x momentum | 7.435304252290109e-7 | 7.435210084283383e-7 | -9.416800672568115e-12 | 3.332698475744245e-6 | 3.332689058943572e-6 |
| fine (8,4,8), y momentum | -7.342303068869072e-7 | -7.342264856102787e-7 | 3.821276628457326e-12 | -3.281185661951902e-6 | -3.281181840675274e-6 |

The endpoint discrepancy is already present immediately before root closure;
the direct increment difference is roughly five orders smaller. This rules
out the final direct reflux increment as the dominant explanation at these
two points. It does not exclude previous synchronization influencing donor
states. The next diagnostic should compare pre-closure second-step evolution,
including fine-stage halo states and ancestor time interpolation, rather than
modify reflux from PeleF-only increments. `sync-comparison.json` retains all
36 comparisons. No numerical correction or qualification claim follows yet.

## Native pre-spatial stage states isolate the second-step predictor

An isolated `Advance.cpp` trace records fine-level Sborder immediately before
each of the three MOL spatial evaluations (old, predictor, post-reaction
corrector). It records the four target valid cells plus two upper-z ghosts.
The original Advance.cpp hash is checked against the pinned manifest; build
commands and hashes of all reused link inputs are retained in
`build/pelec-default-stage-20261023/build-record.json`.

The native run completed in 5.331149647012353 seconds. All four exported
ordinary tables remain byte-identical to captured reference exports. The
native reflux trace hash also matches the preceding diagnostic. There are
216 stage records: six fine intervals, three stages, four probes and three
z positions. `analyze.py` matches 72 valid-cell stage records with the saved
PeleF trace by start time, interval, stage and probe, converting momentum
from CGS to SI. The native ghost records are retained for the next comparison;
the corresponding PeleF ghost states are not yet captured.

At fine (5,8,8), x-momentum difference evolves as follows:

| Fine interval start | Stage | PeleF minus native (kg/m2/s) |
| --- | --- | ---: |
| 100 ns | old state | 2.361044992849060e-8 |
| 100 ns | first predictor | 8.232241091254977e-6 |
| 100 ns | post-reaction corrector | 3.615892539832188e-6 |
| 150 ns | old state | 3.893500770013603e-6 |
| 150 ns | first predictor | 3.303689617983885e-6 |
| 150 ns | post-reaction corrector | 3.331783934634735e-6 |

The first predictor of the second root step is the main newly localized
amplification at this point. It precedes root reflux and cannot be inferred
from endpoint state alone. The next focused comparison is its old-state
spatial evaluation and coarse/fine ghost data at start time 100 ns. Momentum
has no direct chemical production, but the spatial operator includes molecular
transport; compare complete contributions before attributing the difference
to hydro or a particular interpolation policy. `stage-comparison.json`
retains all matched component differences. C2 remains open.

## Matched valid and upper-z ghost states

The first PeleF halo trace attempt preserved all eight ordinary tables but
recorded no halos: the default caller omits optional `gas_old_halos`.
This is a diagnostic implementation failure, not evidence of absent halos.
The corrected isolated observer requests a local output and forwards it if
the caller also requests one. No production source is edited.

`build/eb-default-halo-capture-20261023/run.py` completed in
50.25350533053279 seconds, exit zero, with all eight ordinary tables
byte-identical and exactly 216 halo records. Its automatic `analyze.py`
comparison succeeded for four probes, three z locations, three stages and
six fine intervals. It matches native stage/clock/interval and converts
native conserved state and held energy to SI. Source/build/run receipts and
`halo-comparison.json` retain provenance and the matched vectors.

At start time 100 ns, stage zero, probe 13 (fine 5,8,8), differences
(PeleF minus native) are:

| Quantity | Valid cell k=8 | Upper ghost k=9 | Upper ghost k=10 |
| --- | ---: | ---: | ---: |
| Density (kg/m3) | 4.848246804023404e-11 | 9.433717695905841e-12 | same |
| x momentum (kg/m2/s) | 2.361044992849060e-8 | 8.419002051063273e-9 | same |
| y momentum (kg/m2/s) | 2.390512752678542e-8 | -3.234679571417765e-8 | same |
| z momentum (kg/m2/s) | 2.164357205638361e-9 | 2.492581030538396e-8 | same |
| Held energy density (J/m3) | 1.490121649112552e-4 | -5.424971459433436e-4 | same |

These small initial differences coexist with the much larger subsequent
predictor x-momentum discrepancy. They do not establish its source or rule
out nonlinear reconstruction amplification. This point also lies on the
upper-y coarse/fine interface: the z-only diagnostic does not capture its
complete spatial stencil. The next comparison needs the primitive/flux
reconstruction stencil in all relevant directions, including upper y,
before modifying interpolation or EOS policy. Full C2 agreement remains open.


## Full surrounding stencil comparison

The isolated diagnostic now records a 5 by 5 by 5 cube around each of the
four probes at all three stages of six fine intervals. The native run in
`build/pelec-default-stencil-20261023` preserves all four captured ordinary
output tables. The PeleF run in `build/eb-default-stencil-20261023` completed
in 50.57195531018078 seconds, exit zero, preserving all eight ordinary field
and budget tables byte-for-byte. Its binary SHA256 is
`f3da0b6c76d51bfe1ff6e3be80bf8234dba4487aac1aad7454da2ebca5e8c69f`.
Automatic analysis matches 9000 records by interval, stage, probe and cell,
with native conserved quantities converted from CGS to SI. The comparison
includes ghost and corner locations; these are stencil data, not 9000
independent physical validation samples.

At start time 100 ns, stage zero, the largest absolute component differences
in the cube around probe 13 are:

| Component | Fine index | Signed PeleF minus native, SI |
| --- | --- | ---: |
| Density | (7,9,6) | -4.516360480044668e-10 |
| x momentum | (3,8,8) | -3.179166385702814e-8 |
| y momentum | (7,9,6) | -7.549974639087242e-7 |
| z momentum | (3,6,8) | -1.494404141295851e-7 |
| Total energy density | (7,9,6) | -1.653873769100756e-3 |

These extrema alone do not identify the source of the subsequent predictor
error. The next diagnostic separates native transport flux from the added
hyperbolic flux at the six faces of each probe, before comparing direction
by direction with the saved PeleF stage trace. No numerical policy has been
changed and C2 remains open.

## Native transport/hydro face separation

`build/pelec-default-flux-20261023` compiles an isolated `Diffusion.cpp`
observer, retaining the prior stage and root-reflux observers. It records
face flux immediately before hydro addition and after hydro addition, before
EB divergence/centroid operations. The pinned original source SHA is checked
against the manifest; build receipts hash reused objects and archives.
The serial run completed in 5.380891811102629 seconds, exit zero. All four
ordinary exported tables remain byte-identical to captured reference, and
stage/reflux trace hashes match the preceding diagnostic. No pinned source
or production numerical policy was modified. Binary SHA256:
`b0c84d7c35ebc9a8a4a811da3c1fbb96035ecaad6c9b8fb4c0c70047ee4118c7`.

`analyze.py` verifies 864 unique face records (18 events, four probes,
six faces, before/after hydro), matches the PeleF stage clocks, divides native
area-integrated momentum flux by its CGS face area and converts to Pa.
Subtracting pre-hydro flux isolates the native hyperbolic contribution.
At probe 13, fine (5,8,8), old-state evaluation at 100 ns, all six apertures
are one. Its x-momentum predictor difference increment is
8.208630641326486e-6 kg/m2/s. Applying minus divergence times 50 ns to the
matched hyperbolic flux differences gives:

| Direction | Predicted x-momentum difference increment (kg/m2/s) |
| --- | ---: |
| x | 8.184388182125986e-6 |
| y | 1.156903896497141e-8 |
| z | 1.283831438058726e-8 |
| Sum | 8.208795535471545e-6 |

The remaining increment is -1.648941450589858e-10 kg/m2/s. It has not been
independently decomposed into transport and other contributions. The x-face
term accounts for approximately 99.7% of the observed increment; its low/high
normal momentum flux differences are +0.001108198644942604 and
-0.1011966536316322 Pa. This localizes the main amplification to hyperbolic
normal momentum flux, but does not yet distinguish reconstructed pressure,
velocity, sound speed or Riemann-solver effects. Next compare the actual
left/right face reconstructions at these faces. The result is a local
causal diagnostic, not C2 trajectory qualification or physical validation.

## Face reconstruction and energy basis

The isolated native `MOL.cpp` observer in
`build/pelec-default-riemann-20261023` records the actual left/right states,
cell pressures, acoustic slopes, sound speeds and flux immediately after
Riemann evaluation. It completed in 5.431222690269351 seconds, exit zero;
all four ordinary output tables remain byte-identical, and preceding stage,
reflux and separated-flux trace hashes are unchanged. Analysis matches 432
face records. Binary SHA256:
`79733b7859f6143f588f50098a3edd28ae9e24d72fb78e695575ad4cd463d8b7`.
At probe 13, 100 ns old-state evaluation, the upper-x face left/right
pressure differences are -0.09421877874410711 and -0.05952265253290534 Pa.
The cell-center pressure difference is approximately -0.047036675 Pa.

Source inspection identifies a concrete difference in the primitive recovery:
pinned `Source/PeleCUtilities.H::pc_ctoprim` uses UEDEN minus kinetic energy.
PeleF `compute_mol_eb_hydro_3d` passes its separately held internal energy
to `reactive_conserved_to_primitive`, which overrides total-minus-kinetic
recovery. This is not an inference from pressure magnitude alone.

`build/eb-default-energy-replay-20261023` recompiles a small EOS replay against
the existing selected runtime and LiDryer thermo, with the same pinned EOS
policy and unconverged inverse setting as the default probe. It reads the
saved (5,8,8) state from both solvers and evaluates both energy bases:

| Saved state | Total minus kinetic pressure (Pa) | Held-energy pressure (Pa) |
| --- | ---: | ---: |
| PeleF | 101298.78022729201 | 101298.73312515585 |
| Native | 101298.78016182895 | 101298.73306347668 |

The native observed pressure is 101298.78016183084 Pa; replay agrees within
1.9e-9 Pa. The PeleF held-energy replay agrees with its observed pressure
within 6e-11 Pa. For the PeleF state, held minus total-minus-kinetic energy
is -0.15441959074814804 J/m3. Switching only the replay energy basis changes
pressure by +0.04710213616 Pa, leaving a +6.5461e-5 Pa difference from native.
Thus energy-basis selection explains this cell-pressure discrepancy. This
replay does not prove the full reconstructed-flux or trajectory effect.

The next bounded implementation comparison should change only the primitive
recovery energy basis in an isolated MOL diagnostic, while retaining held
energy transport and transactional checks. Establish an unchanged-module
control first, then compare the candidate against the original fixed field
and conservation gates. Do not remove held-energy transport or adopt a global
EOS change from this single-cell replay. Production numerics remain unchanged;
C2 and physical qualification remain open.

## Primitive total-energy intervention

`build/eb-default-total-energy-20261023` first compiles the unchanged current
MOL hydro module against the frozen diagnostic executable. The control
completed in 50.92350962385535 seconds and all eight ordinary field/budget
tables remain byte-identical to the historical default. The candidate only
omits held energy from primitive recovery; it retains held-energy transport
and the original wall-pressure path. The original comparison checker completed
in 51.82464405708015 seconds (return code one for failed acceptance).
No original tolerances or reference data were changed.

| Time | x momentum max error | y momentum max error | z momentum max error |
| --- | ---: | ---: | ---: |
| 100 ns | 4.756288879415038e-8 | 8.108242019466694e-8 | 3.419276425911047e-7 |
| 200 ns | 1.288247959418953e-7 | 1.831170937105497e-7 | 2.815275139467425e-6 |
| 300 ns | 1.897653612026318e-7 | 2.269973269353187e-7 | 5.025297557363786e-6 |

Units are kg/m2/s; the unchanged momentum limit is 1e-6. Four violations
remain: z momentum at 200/300 ns and total/held energy at 300 ns (relative
1.4624434529211134e-8 / 1.4622681775748149e-8, versus 1e-8). All other
completed checker gates pass, including species/element inventory and the
correction-accounted energy budget (maximum residual magnitude 2.72e-17).
This is an isolated candidate, not production or full C2 qualification.

The late z-momentum maximum is fine (7,6,6), and the late energy maximum is
fine (2,6,6), both cut cells with volume fraction 0.72. Code inspection shows
native EB wall pressure uses the same total-energy-derived primitive as its
face flux. The candidate still passes held energy into both PeleF wall
contribution calls, so it mixes energy bases. The next controlled candidate
omits that argument in both wall calls while retaining held-energy transport.

Integration must preserve `test_eb_mol_hydro_3d::check_held_input`, which
explicitly tests held-pressure recovery under a large rotating potential and
late invalid-held-state rollback. An unconditional global replacement is not
justified by this diagnostic. The API/policy boundary and its lifecycle require
separate implementation and tests after the numerical comparison.

## Consistent primitive and wall energy basis

The isolated `build/eb-default-total-wall-20261023` candidate additionally
uses total-minus-kinetic energy for both EB wall contribution calls, retaining
held-energy transport. The original checker completed in 50.37314318306744
seconds, return code zero, with no violations through 300 ns. Binary SHA256:
`291b54b32f59fa4d9f46c41fe75bd155f409b80ec3e15ea3ea46df741c16cacc`.
At 300 ns, maximum x/y/z momentum errors are respectively
2.041672231278824e-8, 2.150483202623121e-8 and 8.403514433147699e-9 kg/m2/s;
maximum temperature error is 7.495609906982281e-5 K. Maximum magnitude of
correction-accounted energy-budget residual is 9.79e-17. Original density,
energy, species, inventory and H2-consumption gates also pass. `comparison.json`
retains every field error and the original failed baseline for comparison.

This establishes bounded fixed-reference agreement for the isolated frozen
serial probe. It does not establish a current-source application, MPI/restart,
spatial convergence, physical validation or complete C2 qualification.

The current MOL kernel now has an optional `primitive_from_total_energy`
argument (default false) selecting the same energy basis for face primitives
and EB wall pressure. Held-energy transport is unchanged. The focused unit
adds explicit total-energy comparisons with a deliberately different held
scalar across the existing cut geometry, rotations, orders and solver modes.
Existing held recovery and late-failure rollback checks remain intact. The
application/tree caller is not yet connected; focused build/test results
must be recorded separately when completed.

The focused current-source build completed with GNU Fortran 13.3.0 in
`build/eb-eos-application-20261021/isolated-openmpi/registration`.
`ctest --test-dir <that-directory> -R '^unit_eb_mol_hydro_3d' --output-on-failure`
passed 1/1 tests in 3.39 seconds, including the new total-energy checks and
existing held/rotating rollback checks. Build warnings remain in the log;
this was not a clean full configuration matrix.

The ordinary nonrotating `compute_mol_eb_tree_hydro_3d` caller now explicitly
selects total-energy pressure in all eight kernel call branches, including
spatial halos. Rotating paths and explicitly selected experimental
`separate_momentum` halo-policy paths retain held-energy pressure; those are
separate numerical policies. No new application option is introduced.
Current-source tree and selected-reference targets are building. Their
trajectory, MPI and restart results are not yet established by the isolated
comparison above.

The rebuilt current-source `tests/test_eb_tree_reacting_reference` also passed
the original fixed-reference checker through 300 ns with zero violations.
Its report is under `build/eb-current-total-energy-20261023/comparison/`.
This is stronger than the frozen-object intervention but remains a serial
reference probe, not full application/MPI qualification.

The separately rebuilt `unit_eb_tree_hydro_3d` failed its pre-existing
"cold Riemann star rejects below mechanism range" assertion at line 109,
in a CTU `advance_reactive_eb_tree_hydro_3d` call. Attribution is unresolved;
an unchanged-module control is required. No rejection condition or tolerance
has been relaxed. This failure prevents claiming the broader tree checks pass.
MPI direct kernel calls now select the same energy basis as serial, with
frame presence included in collective control agreement. MPI/application
targets are building; their tests have not yet passed.

## CTU rejection control

`build/eb-ctu-rejection-control-20261023` freezes the current fixed core archive
and test object, recompiles the pre-change MOL kernel plus a copy of the tree
module with only this change's energy-selection routing removed, and relinks
without modifying the worktree. Build commands and every input hash are saved.
The resulting control binary SHA256 is
`cfd426f26a8b79258efd0d835aebf4c1346c059fcdc1dcc0d40cd1f5a6a85ba6`.
It fails the identical cold CTU rejection assertion at line 109, return code
one, after 2.272438731044531 seconds. Thus this assertion failure predates the
present MOL energy-basis intervention within the current dirty source state;
it is not evidence that the present change caused it. This is not a clean
historical-commit baseline or a passing tree test.

Inspection finds current `reconstructed_pelec_state` explicitly permits
positive-temperature polynomial evaluation for transient face/star states
while retaining strict temperature checks for accepted cells. This gives a
specific candidate explanation for the old cold-star rejection expectation,
but the failing fixture's actual transient/accepted temperatures have not
been traced. Do not change that assertion merely to make the test pass.

The isolated temperature observer in `build/eb-ctu-temperature-trace-20261023`
retains the failing assertion and all numerical operations. It prints the
first 16 out-of-fit transient evaluations and accepted fluid-cell temperature
ranges immediately before the assertion. It reproduces the assertion failure
(return code one, 2.2218191102147102 seconds). The transient values are about
284.4554 K, below the common fitted minimum 300 K. The accepted fluid cells
are 365.14748189174264 to 1095.4545446548357 K, within the 300--3500 K bounds.
Binary SHA256:
`81b67d7a8ebadb4c5917f05810162d5a45ea88c449990e01d2433e55a34dba3d`.
Thus the old assertion specifically conflicts with existing transient-state
extrapolation; this observation does not show acceptance of an out-of-range
cell. A justified test correction must verify accepted-cell range and
conservation, and separately retain rejection/rollback for an actual invalid
cell, rather than removing temperature admission coverage.

## Current serial/MPI lifecycle checks

The rebuilt fixed MPI unit and selected serial/MPI application targets passed
all seven focused tests in 284.80 seconds: two- and three-rank MOL spatial
and MOL reflux checks (four tests), selected reacting restart with common
and recursive clocks (two tests), and selected held application lifecycle
(one test). Exact CTest output is
`build/eb-current-total-energy-20261023/lifecycle-tests.log`.
These exercise serial/MPI agreement, reaction, regrid and changed-rank restart
under their bounded fixtures; they are not an MPI fixed-PeleC trajectory
comparison or a full configuration/install audit.

The separate CTU test correction now checks cold-step conservation and strict
recovery of accepted fluid cells. Its rejection fixture constructs a real
cell at half the fitted lower temperature using transient polynomials, first
confirms strict EOS refusal, then verifies tree/temperature and boundary/wall
rollback. This preserves an actual cell-admission failure test while aligning
the transient-state expectation with the documented current implementation.
The corrected test's execution result remains pending.

The corrected `unit_eb_tree_hydro_3d` passed 1/1 tests in 46.42 seconds.
It now checks strict accepted-cell temperature range, cold-step conservation,
and explicit invalid-cell rejection/rollback without changing solver behavior
or numerical tolerances. The complete output is
`build/eb-ctu-temperature-trace-20261023/test-update.log`.

The checker now has a separate `--baseline-mesh-study` route for mesh/time
studies using the original mechanism and unchanged pinned native binary.
It binds the probe/exporter to a passed original comparison's recorded binary
hashes, validates the native reference manifest and every captured artifact,
and rejects experimental halo policies or extra mechanism/EOS/profile controls.
The original comparison contract remains unchanged; study reports remain
explicitly non-qualification evidence. Two negative admission checks reject
an extra mechanism control and a changed binary hash before creating a solver
work directory (`build/eb-pinned-mesh-20261023/admission-checks.json`).
An 8-cell one-step smoke comparison is running to verify that the new route
reproduces the existing first-step tables before any larger pilot.

The original-mechanism 8-cell smoke completed and passed the unchanged gates.
All six exported native/PeleF field and PeleF budget tables at 0 and 100 ns
are byte-identical to the passed ordinary comparison. The identity receipt is
`build/eb-pinned-mesh-20261023/smoke-identity.json`. A 16-cell, one-root-step
100 ns pilot is running with a 300-second per-command limit, retaining the
original frozen 33-point profile and binaries. It is a cost/error pilot,
not a completed convergence study. The mesh comparator now names this native
solver `native_pinned` and rejects mixing original and diagnostic mechanism
bases; existing diagnostic `native_held` output naming is preserved.

## Original-mechanism mesh pilot

The 16-cell one-step 100 ns pilot completed with zero original field/budget
violations. Native execution took 12.753728490322828 seconds; PeleF took
82.66807485558093 seconds. Maximum x/y/z momentum differences are
1.724881503600839e-8, 6.800146779006866e-9 and
2.734364539999645e-8 kg/m2/s; maximum temperature difference is
2.336418742743263e-5 K. The energy-budget residual is
4.435636690971158e-17. `build/eb-pinned-mesh-20261023/spatial8-16.json`
compares both solvers using conservative fluid-volume restriction at 100 ns,
separating initial sampling, endpoint and evolution increments. The y-momentum
evolution-increment RMS mesh difference is 0.0033165336613122085 for PeleF
and 0.0033165336184536763 for pinned native (kg/m2/s). Two meshes do not
establish an observed convergence order.

A 32-cell one-step spatial pilot now has a 1200-second per-command bound;
its anticipated PeleF cost is approximately 8 times the measured 16-cell cost,
subject to actual scaling. An independent 8-cell half-time-step comparison
uses two 50 ns steps to the identical 100 ns endpoint, with 180-second command
bounds. They retain fixed binaries, original mechanism and frozen profile.
The mesh comparator now reports both RMS and maximum-norm observed orders
when three grids/timesteps are available; no order is assigned for two grids.

The fixed-original 8-cell half-time-step run (two 50 ns steps to 100 ns)
passed every reference gate. Native execution took 2.3212796933948994 seconds;
PeleF took 31.36000239290297 seconds. The full-versus-half-step y-momentum
RMS difference at the common endpoint is 1.0161037910291383e-4 kg/m2/s for
PeleF and 1.0161033737871534e-4 for native; maxima are respectively
6.292002606469727e-4 and 6.292001160599914e-4. These are time-refinement
self-differences, not PeleF-minus-native errors. The two-step comparison
report is `build/eb-pinned-mesh-20261023/temporal-full-half.json`.
A four-step 25 ns run to the same endpoint is now pending for an observed
order estimate. The 32-cell spatial pilot remains separate and unfinished.

## Three time increments at 100 ns

The 8-cell quarter-time-step run (four 25 ns steps) passed every original
reference gate. Native and PeleF execution took 2.371803779155016 and
59.246976589784026 seconds. The 100/50/25 ns comparison to the same 100 ns
endpoint is `build/eb-pinned-mesh-20261023/temporal-three.json`.
All sixteen fields have positive observed RMS and maximum-norm orders.
For PeleF, momentum RMS orders are x=1.25697365, y=1.25581096,
z=1.20308521; maximum-norm orders are x=1.09149922, y=1.09344341,
z=1.15125381. Pinned native yields the same momentum orders within 2e-6.
Density RMS order is 1.12211033 and temperature 1.15680114. Species orders
span approximately 0.93--1.70 in RMS. This is decreasing temporal error
consistent with the fixed reference, not observed second-order or asymptotic
qualification. It is still only an 8-cell, 100 ns endpoint result.

A six-step 50 ns run to the original 300 ns comparison endpoint is now pending
to check whether the behavior persists over the full reference interval. The
original 100 ns-step, 300 ns baseline already exists and should be reused with
explicit provenance rather than recomputed just to supply study metadata.

The six-step 50 ns run to 300 ns passed every original reference gate
(native 5.680205721408129 s, PeleF 93.90389463864267 s). The comparator can
now reuse a passed original capture when a successful original-mechanism
study binds its exact checks hash, mechanism/source contracts, native
manifest and probe identity. It derives only comparison metadata, retains
original status and source files, and identifies the reused receipt in output.
Four focused admission checks cover valid reuse, report-hash mismatch, probe
mismatch and unequal endpoint rejection. The original 300 ns calculation was
not rerun. The two-time-step result is
`build/eb-pinned-mesh-20261023/temporal300-full-half.json`.
A twelve-step 25 ns run is pending for the full-interval order estimate.

## Time refinement over the full 300 ns interval

The twelve-step 25 ns run passed all original reference gates (native
6.7863089218735695 s, PeleF 178.06377811916173 s). The 100/50/25 ns-step
comparison at 300 ns is
`build/eb-pinned-mesh-20261023/temporal300-three.json`; it reuses the original
100 ns-step baseline with explicit provenance. All sixteen fields have
positive observed RMS and maximum-norm orders. PeleF momentum RMS orders
are 1.18885879, 1.18888698, 1.21602846 (x/y/z), and maximum-norm orders are
1.03970115, 1.03575353, 1.16496674. Native momentum orders differ by less
than 7e-7. Density and temperature RMS orders are 1.12535056 and 1.20697763.
Species RMS orders span 0.957--1.288. Thus the full reference interval shows
monotone temporal self-difference reduction closely matching the original
native solver, but not second-order or asymptotic accuracy. This remains one
8-cell geometry and does not close spatial or general C2 qualification.
The 32-cell spatial pilot is still executing; no spatial order is yet claimed.


## Completed 32-cell spatial pilot: reference failure

The original-mechanism 32-cell, one-step 100 ns pilot completed both solvers
but failed six unchanged reference gates. The authoritative capture is
`build/eb-pinned-mesh-20261023/pilot32/matched_eb_amr_niwcw8rh/checks.json`;
localization and input hashes are in
`build/eb-pinned-mesh-20261023/localization32.json`.
Momentum maximum differences (x/y/z) are 3.9535115663e-6,
1.3283232260e-5, 3.7295295489e-6 kg/m2/s against a 1e-6 limit.
Scaled density, total-energy and held-energy differences are
2.5676561413e-8, 1.9097850095e-7 and 1.9088387161e-7 against 1e-8.
The budget residual remains 4.9419963358e-17.

The y-momentum and energy maximum is at fine cell (25,32,26), x-momentum
at (25,32,25), density at (26,32,30), and z-momentum at (15,23,32).
All have volume fraction one and lie against an upper coarse/fine interface.
The boundary-layer share of volume-weighted squared y-momentum error is
84.75%, versus approximately 0.000142% at cut cells. This localizes the
failure; it does not yet distinguish interpolation from flux synchronization.
Native and PeleF executions took 83.6692492 and 536.8650333 seconds;
whole-checker elapsed time was 10:40.04 and peak RSS 1991280 KiB with no swap.
The completed 8/16-cell comparisons do not establish 32-cell reference
agreement. No accepted three-grid spatial order or C2 completion is claimed.

The probe now accepts `--default-reacting-trace` to attach the existing MOL
stage and synchronization observers to the ordinary reacting trajectory.
It preserves the existing shear probes and adds IDs 10--13 for these four
32-cell extrema. Outputs include reconstructed left/right pressure,
directional raw/centroid/tree momentum flux, local rates and synchronization
increments. This is observational instrumentation, not a numerical fix.

Instrumentation verification used GNU Fortran 13.3.0 Debug and the existing
selected LiDryer runtime archives, linked to a separate executable without
overwriting the captured reference probe. Exact compile/link commands and
source/binary SHA256 are in
`build/eb-default-trace-current-20261023/build-record.json`.
The 8-cell, one-step 100 ns traced run completed with exit zero; both state
and both budget CSVs are byte-identical to the saved untraced smoke capture.
It emitted 216 stage rows and 324 synchronization rows, recorded in
`build/eb-default-trace-current-20261023/smoke-checks.json`.
Mechanism provenance, project-contract and `git diff --check` passed.
This checks observational neutrality at 8 cells; 32-cell hotspot traces and
their native comparison are not yet collected. No new physical validation
or full configuration matrix was run for this instrumentation change.


The bounded 32-cell traced trajectory now runs separately in
`build/eb-default-trace32-20261023`, using the exact saved initial data and
verified instrumented binary. Native diagnostics were rebuilt from pinned
Advance/Diffusion/MOL sources with the four matching zero-based coordinates;
compile/link identities are under `build/pelec-default-{stencil,flux,riemann}32-20261023`.
The first native diagnostic launch omitted the original mesh CLI overrides
and failed initialization; its `record.json` and `run/run.log` are retained.
The corrected `run-matched.py` copies the complete argument tail from the
original pilot's checks record, and writes separate `record-matched.json`
and `run-matched/` artifacts. Successful completion must also establish
identity of both exported native states with the original capture.
At this checkpoint both traced calculations are pending; no new stage
comparison or neutrality result is claimed. Resume instructions and tool
handles are in `build/eb-default-trace32-20261023/handoff.json`.


The corrected native 32-cell diagnostic completed in 83.6224701 s, exit zero.
Both initial and 100 ns exports (65536 rows each) exactly match their saved
original-reference SHA256 values; see `record-matched.json` under
`build/pelec-default-riemann32-20261023`. The traces contain 3000 stencil,
288 before/after-hydro face-flux and 252 Riemann rows. MOL's local event
counter advances per tile (eight calls per spatial evaluation at this grid),
whereas the independent Diffusion counter advances per spatial evaluation.
After normalization, all 144 expected event/probe/axis/side keys are present;
repeated tile records are identical, and momentum flux reconstructed from
the independent before/after records agrees within 2.78e-17 Pa.
The normalization receipt is `riemann-normalization.json`.
PeleF's traced calculation is still pending at this checkpoint. Analysis now
requires a completed, byte-identical four-table PeleF neutrality receipt;
partial stage output alone cannot establish neutrality.


## 32-cell stage comparison completed

The traced PeleF run completed in 533.7645789 s, exit zero. Its initial/final
state and budget tables all exactly match the original pilot hashes
(`build/eb-default-trace32-20261023/neutrality.json`). Both instrumented
trajectories are therefore neutral for the saved outputs. No jobs remain
from this diagnostic pair.

The six fine spatial evaluations are now matched by interval, stage, probe,
axis and side in `build/pelec-default-riemann32-20261023/flux-comparison.json`
and `riemann-comparison.json`. First old-state hydro-flux differences are at
most 7.96e-9 Pa. At the first predictor evaluation the maximum becomes
0.0168627057 Pa at probe 10, the upper y face of cell (25,32,26).
At this face the reconstructed pressure differences are only
5.3693e-5 Pa (left) and 1.9121e-4 Pa (right). The valid-cell y-momentum at
this predictor still differs by only 3.67905e-11 kg/m2/s, from identical
initial momentum. Thus the discrepancy is already present in the predictor
spatial operator, before final coarse/fine synchronization; the final reflux
operation alone cannot explain its onset. Actual prepared halo states and
reconstruction slopes remain to be compared before assigning the cause to
interpolation or a limiter.

`momentum-balance.json` compares saved stage-state differences with the
predictor/trapezoidal hydro contributions from all six faces. At probe 10,
the first corrected-state y-momentum difference is 6.1957138e-6 kg/m2/s;
the hydro contribution is 3.5947955e-6, almost entirely from y faces.
The 2.6009183e-6 residual is unresolved and includes transport/other update
terms; this is not a closed total-source balance. Pressure then diverges
inside the cell as well (maximum across four probes 0.010606 Pa at the
first corrected evaluation). Later-stage differences must not be mistaken
for the initial cause.

The comparison/normalization/balance scripts pass Python compilation;
`git diff --check` passes. No numerical correction, new reference pass,
spatial convergence verdict or physical qualification follows from this
localization alone. The next diagnostic is the actual predictor halo and
face-reconstruction stencil, with saved native stencils available for reuse.


A separate predictor-halo diagnostic is prepared under
`build/eb-predictor-halo32-20261023`. It builds a copy of the current tree
module, captures the actual gas halo output even when the caller omits that
optional output, and forwards it unchanged when requested. The observer
writes conserved and held-energy stencils around all four extrema at the
first fine old/predictor evaluations, then deliberately stops with code 77
and `DIAGNOSTIC_PREFIX_CAPTURE_COMPLETE`. This is a bounded prefix capture,
not a full trajectory pass. Its analyzer requires exact initial-table hashes
and stage-log prefix agreement with the completed neutral run before using
the saved native stencils. Production numerical source is unchanged; build
and execution status are retained in the diagnostic directory. At this
checkpoint compilation is pending; no halo comparison is yet claimed.


## Predictor halo and sibling reconstruction

The bounded diagnostic completed its intended prefix in 305.8020874 s,
returning the documented stop code 77 and completion marker. Initial state
and budget hashes and the stage-log prefix exactly match the neutral full
trajectory. This is prefix neutrality, not another full-trajectory pass.
`build/eb-predictor-halo32-20261023/halo-comparison.json` compares 1000 rows
with the saved native stencil. Initial valid cells agree exactly; initial
halo momentum differences are below 2.8e-16 kg/m2/s. At the predictor,
maximum y-momentum difference is 4.07695e-5 in the halo versus 1.04952e-8
inside the sampled fine region. The difference thus precedes face
reconstruction and is already present in prepared coarse/fine halos.

`sibling-audit.json` reconstructs interpolation centers and directional
slopes from complete sets of eight fine ghost siblings, without another
solver execution. For coarse donor (13,17,13), adjacent to the primary
maximum, the center y-momentum difference is only 3.34476e-11 kg/m2/s.
The x-direction y-momentum slopes are 0.0096649320 (PeleF) and 0.0098216362
(native), about 1.6% apart. Density and x-momentum slopes show similar
relative scaling. Both y-direction slopes vanish. This points to spatial
slope construction/limiting rather than a large mismatch of the donor
center. It does not identify which neighbor or component causes the scale.
The default interpolator shares one limiter scale across packed components;
center/six-neighbor inputs and per-component limiting ratios are the next
required evidence. Sensitivity to small species values is a hypothesis,
not a demonstrated cause or justification for changing tolerances.
The analysis scripts pass Python compilation and `git diff --check` passes.


The next isolated diagnostic records the actual center/six-neighbor inputs
for donor (13,17,13), before spatial interpolation. Native observation is
inside the pinned `MFCellConsLinInterp::interp`, not a separately guessed
coarse snapshot; `Source/Setup.cpp` selects the EB MF interpolator, whose
regular-cell path uses this routine. Both pinned MF and current default
PeleF implementations take the minimum component limiting ratio per axis.
The new directories are `build/eb-predictor-donor32-20261023` and
`build/pelec-predictor-donor32-20261023`. Native was built and its neutral
trajectory check is running; PeleF compilation is pending. The component
ratio analyzer is prepared and passes Python compilation, but no limiting
component has yet been identified. No production numerical edit was made.


The native donor diagnostic completed in 83.2877844 s; both 65536-row
exports retain exact reference hashes, as do the existing stencil/flux
traces. `build/eb-predictor-donor32-20261023/native-limiter.json` matches
native interpolation event 2 to the saved predictor donor center within
6.94e-18 scaled difference. Recomputed component slopes agree with the
sibling-derived slopes within 7.64e-11 absolute (the largest residual is
in energy). The unique active limiter for x and z is **z momentum**, with
scales 0.9882504302745775 and 0.05756721298658442. The y scale is zero.
This identifies the native controlling component; the earlier small-species
hypothesis is not supported for this donor. PeleF donor capture is pending,
so its active limiter and the cause of any coarse z-momentum difference
remain unproven. No numerical policy change follows yet.


## Shared limiter input attribution completed

The PeleF donor prefix completed in 300.5955548 s with the expected stop77.
Its halo SHA256 is identical to the prior prefix and its initial tables and
stage logs again match the completed neutral trajectory. `matched-limiter.json`
under `build/eb-predictor-donor32-20261023` matches both solvers' centers and
reproduces their sibling slopes (largest absolute residual 7.64e-11 in energy).
Repeated candidates carry identical relevant inputs; they do not change the
identified controller. Both solvers uniquely select z momentum for x/z
limiting. Their x scales are 0.9724828980184832 (PeleF) and
0.9882504302745775 (native); z scales are 0.05822831040327126 and
0.05756721298658442.

`input-ablation.json` performs offline input substitution, with no solver
rerun or policy modification. Replacing only center/six-neighbor z momentum
with native values restores both native limiter scales exactly and reduces
the target ghost y-momentum difference from 3.75826305e-5 to
3.24430205e-11 kg/m2/s. Replacing every other component retains the PeleF
scales and a 3.75825980e-5 y-momentum difference. Thus this input component
explains the shared-limiter amplification at this donor. The largest z
momentum input differences are at the center (-3.59041651e-9) and upper-z
neighbor (-5.66426358e-10 kg/m2/s); other neighbors differ by at most
9.07e-12. The origin of those coarse-state differences remains unresolved.

The next required evidence is coarse old/predictor/corrected z momentum and
separate hydro/molecular-transport contributions at (13,17,13) and
(13,17,14), before the fine advance. Altering the shared limiter or enabling
separate momentum would mask this amplification while changing the pinned
numerical policy; no such numerical change was made. Python compilation
and `git diff --check` passed for this analysis increment. C2 remains open.


### Coarse donor source capture: native control completed

An isolated diagnostic build records coarse cells `(13,17,13)` and
`(13,17,14)` at old, predictor and corrected stages, with molecular and
hydrodynamic face fluxes recorded separately. These cells and their six
immediate neighbors have volume fraction one in the saved initial table.
The native run completed in 83.6812 s; both 65,536-row state exports retain
the original initial/final SHA256 identities. Reconstructing the predictor
and trapezoidal corrected momentum increments from these face fluxes leaves
an absolute residual at most 3.77e-15 kg/m2/s. This checks stage/units and
local flux closure, not PeleF agreement or C2 completion.

Evidence: `build/pelec-coarse-flux32-20261023/record-matched.json` and
`native-momentum-closure.json` in that directory. The companion PeleF prefix
capture is pending; resume instructions and its live handle are in
`build/eb-coarse-source32-20261023/handoff.json`. No production numerical
algorithm was changed by these isolated diagnostic builds.


### Coarse source attribution completed

The PeleF diagnostic prefix completed in 304.005 s with the intentional
capture exit 77. Its initial exports and recorded stage prefix match the
completed neutral trace; the 1,000-row halo capture retains SHA256
`ae1f8f34cead37f4c033a50227391b1f110f7e897ba78dc796e1fb4b5711e927`.
`build/eb-coarse-source32-20261023/compare-coarse.py` compares the six
coarse cell/stage samples against the native capture. Predictor and first
corrected momentum difference balances close to 3.82e-15 kg/m2/s.

At cells `(13,17,13)` and `(13,17,14)`, predictor z-momentum differences
are 2.86e-15 and -3.97e-15. At the corrected state they are 3.40e-11 and
3.60e-11, while the subsequent hydro RHS differences are -0.1436214 and
-0.0226463 kg/m2/s2. The corresponding molecular RHS differences are
4.71e-6 and -1.07e-5. Integrating the old/corrected RHS and applying the
half-time donor interpolation reconstructs the previously measured donor
z-momentum differences -3.5904165e-9 and -5.6642636e-10; see
`donor-balance.json` in that directory. Thus the coarse final hydro evaluation
is the immediate source of the donor discrepancy amplified by the shared
limiter. This does not yet identify which pressure/face reconstruction term
causes that hydro difference. Capture those coarse stage-2 faces and adjacent
states next. C2 and physical qualification remain open.


### Coarse directional faces and pressure

The isolated coarse-face prefix completed in 174.687 s (intentional exit77).
Initial exports and all sampled donor stage states match the previous capture;
all 36 sampled faces have aperture one. The native diagnostic retained both
original export hashes. Its 102 tile records reduce to 36 unique face/stage
records, with independent hydro-flux agreement within 1.11e-17 Pa.

At stage2, z faces contribute -0.1436208375 and -0.0226456206 kg/m2/s2
to the z-momentum RHS differences at the two donors; x/y contributions are
at most 5.07e-7 individually. The maximum face momentum-flux difference is
5.3124e-5 Pa; the maximum sampled cell pressure difference is already
5.1415e-5 Pa before face reconstruction. Therefore audit state-to-pressure
recovery and the chemistry-updated state next, including the temperature
initial guess. Both pinned inverse-EOS implementations stop before the final
update when abs(dT)<1e-6 K; changing this tolerance would change the reference
policy and is not justified by this capture alone.

Evidence: `build/eb-coarse-faces32-20261023/{neutrality,flux-comparison,
directional-attribution}.json` and
`build/pelec-coarse-riemann32-20261023/{record-matched,riemann-normalization,
riemann-comparison}.json`. No production numerical changes or new C2 pass
are claimed.


### Same-native-EOS state/temperature-guess ablation

A standalone executable includes the unmodified pinned LiDryer mechanism
header and calls CKUBMS, CKCVBS, GET_T_GIVEN_EY and CKPY on saved coarse
states. The 24 combinations cover both cells, three stages, and independent
PeleF/native conserved-state and temperature-guess choices. At stage2 the
native EOS reproduces PeleF pressure on PeleF state/guess within 1.46e-11 Pa.
Thus a different pressure formula is not needed to explain this center gap.

For `(13,17,14)`, the PeleF state/guess gives a first Newton correction
-9.74654e-7 K, below the pinned 1e-6 K stop-before-update threshold. The
retained guess gives pressure difference +5.14150e-5 Pa. Replacing only the
guess by the native temperature gives -2.62592e-6 Pa. Conversely the native
state with the PeleF guess gives a correction -1.04618e-6 K, crosses the
threshold and updates to native pressure within 5.83e-10 Pa. This identifies
a threshold-sensitive temperature/state interaction, not justification to
tighten the reference tolerance. Trace temperature publication after the
first chemistry correction next.

Artifacts: `build/eb-coarse-eos32-20261023/` contains `eos.cpp`, `build.json`,
`audit.py`, exact inputs/output, `audit.json` and the next-step handoff.
No solver trajectory was rerun and no production numerical change was made.


### Reaction temperature lifecycle source audit

Both tree reaction drivers start from old-state species and temperature.
However, pinned `ReactorUtils.H::fKernelSpec` integrates a temperature
coordinate alongside every species and uses that coordinate as the EOS guess;
unflatten also recovers temperature from the integrated coordinate. The current
PeleF `sundials_forced_reactor_mod` integrates Nspecies-1 independent densities,
reconstructs a dependent species, and caches the recovered temperature across
RHS and Jacobian evaluations. Its endpoint EOS uses that cache. Thus its
temperature guess can depend on callback history instead of the integrated
native temperature coordinate. Under the pinned finite EOS stopping criterion,
these are distinct numerical paths.

This source difference is a candidate cause consistent with the saved-state
ablation, not yet a demonstrated trajectory fix. Capture the donor reaction
inputs and replay a temperature-coordinate prototype before changing the
production integrator. Preserve generic reduced integration, conservation and
rollback contracts. Source hashes and next steps are recorded in
`build/eb-coarse-eos32-20261023/temperature-lifecycle.json`.


### Exact single-cell replay of coarse reaction inputs

The isolated input capture completed in171.573s (intentional exit77) and
saved four cell calls: two zero-forcing initialization calls and two first
MOL chemistry-corrector calls. The saved old state, old held energy, old
 temperature, complete nonreacting rates, dt and tolerances feed a standalone
`react_sdc_state_3d` + existing reduced CVODE replay. For both corrector calls,
all14 gas components and temperature exactly equal the previously captured
coarse stage2 state. Thus the temperature discrepancy is reproducible without
another full-grid trajectory. No tolerance or production algorithm changed.

Artifacts: `build/eb-coarse-reactor32-20261023/{record,build-replay,
replay-check}.json`, `run/reactor_inputs.csv`, `replay.F90`, `replay.csv` and
`check-replay.py`. Next compare an isolated temperature-coordinate integrator
on these exact inputs; exact replay alone is not a C2 parity pass.


### Temperature-coordinate prototype and Jacobian control

An isolated all-species-plus-temperature CVODE prototype with finite-difference
Jacobian replays the exact saved inputs. Upper-cell pressure difference falls
from5.1415e-5 to1.07383e-6Pa; gas-state changes are below3.82e-15. However,
a reduced-coordinate finite-difference-Jacobian control also shifts the upper
cell temperature by-9.74607e-7K, essentially the same change. Therefore this
improvement does not isolate temperature integration as its cause. The cached
EOS guess is sensitive to solver callback history; isolate the analytic
Jacobian callback's update of that cache next before adopting a fix.

Artifacts: `build/eb-temperature-coordinate32-20261023/{build,comparison,
fd-control,handoff}.json` and `build/eb-reduced-fd32-20261023/`. Both are
standalone diagnostic copies; production integrator and tolerances unchanged.
The first prototype compile failed due to a missing SUNDIALS_INT32_T define;
rebuilding with the installed Fortran interface's integer width succeeded.


### Callback chronology explains the temperature cache effect

An analytic-Jacobian control disabling only the Jacobian callback's temperature
cache write exactly reproduces both baseline cell outputs. That write alone
is not the cause. Neutral callback traces reproduce each baseline/FD replay
byte-for-byte. In the upper cell the analytic path evaluates at99.98588245ns,
then100ns; the final evaluation retains the previous1961.5287759100438K under
the pinned EOS stopping threshold. The FD path instead evaluates at98.92668629ns
then100ns and updates to1961.5287749354363K. There are61 versus71 logged
callbacks for this correction. Adaptive evaluation chronology therefore
changes the cached EOS guess even when gas-state differences are at roundoff.

Artifacts: `build/eb-jacobian-cache32-20261023/{comparison,
callback-comparison,handoff}.json` and the neutral trace builds
`build/eb-cache-trace-{analytic,fd}32-20261023/`. The appropriate next test is
a native-compatible temperature-coordinate path on the same saved inputs,
with matched solver settings and conservation/rollback evidence; neither
loosening EOS tolerance nor choosing a Jacobian because this cell improves
establishes fixed-reference fidelity. Production code remains unchanged.


### Pinned ReactorCvode same-input reference

A standalone driver now calls the pinned ReactorCvode box API using the
saved PeleF old states and forcing, converted to CGS. It reuses the actual
compiled PelePhysics objects, native denseAJ_direct, rtol1e-10 and CGS
atol1e-18; typical-value scaling remains disabled. Four calls completed.
Against this same-input native endpoint, baseline upper-cell temperature
difference is1.00860e-6K and the temperature-coordinate prototype differs by
3.40333e-8K. Species density differences remain approximately2.74e-12kg/m3.
The V2 PeleF tolerance is density-scaled and not exactly native CGS1e-18;
matching this and the Jacobian configurations is still required. This is
not a full-trajectory parity result or production adoption.

Artifacts: `build/native-reactor-replay32-20261023/` contains the driver,
build command record, exact input conversion, native outputs and comparison.
The first driver run failed from omitted AMReX SUNDIALS context initialization;
adding the same Initialize/Finalize lifecycle as native main fixed it. The
successful run is separate from that retained failure backtrace.


### Solver-setting controls do not explain remaining species gap

Sixteen native API replays cover four inputs, analytical/finite-difference
Jacobian, and original/matched PeleF-density absolute tolerance. For the two
corrector calls these controls change temperature by at most5.48e-11K and
species density by5.32e-15kg/m3. They do not explain the approximately
2.74e-12kg/m3 species discrepancy. A separate all-species-plus-temperature
Fortran prototype removes its max-step and explicit stop-time settings to
follow native CPU CV_NORMAL behavior; the discrepancy remains, with temperature
difference about3.40e-8K against native FD/matched-density-atol results.
Temperature absolute tolerance and further solver details are not asserted
identical. Compare same-state reaction RHS and thermodynamics next.

Evidence: `build/native-reactor-replay32-20261023/controls.py`, `controls.json`,
individual control directories, and
`build/eb-temperature-native-clock32-20261023/comparison.json`.
These are diagnostic controls, not adopted numerical fixes or C2 completion.

### Same-state kinetics and activation-basis control

Sixteen identical SI states (four initial states and native, baseline PeleF,
and temperature-coordinate endpoints) were evaluated with pinned native
CKWYR and PeleF forced_reactor_rhs, with zero forcing. Recovered temperatures
agree exactly in all sixteen states; cv and internal energy differ only by
approximately 2.3e-12 J/kg/K and 2.8e-9 J/kg respectively. Canonical chemistry
mass-rate differences reach 4.57537e-5 kg/m3/s (initial states 5.45010e-6).
Thus a species discrepancy exists before time integration.

A separate executable links the previously generated, provenance-recorded
native-activation diagnostic mechanism, without changing canonical source or
generated rates. The maximum endpoint rate difference falls to 6.61803e-9
kg/m3/s, and the initial difference to 3.60155e-11. Temperature still agrees
exactly. The native activation conversion therefore explains most of this
same-state rate gap; the residual has not been attributed. This control does
not establish endpoint, full-trajectory, convergence, or physical validation.
Production numerical sources remain unchanged.

Artifacts: `build/eb-same-state-rhs32-20261023/` contains both drivers, exact
states, build commands and binary hashes, `comparison.json`, and the separate
`activation_control.py` / `activation-control.json`. The activation source is
`build/eb-native-activation-20261020/verified-regeneration/lidryer_native_activation.F90`;
its hash and link commands are recorded in the control. Next compare saved-input
reaction endpoints with this explicit activation control before considering
any production change; keep temperature-cache chronology as a separate cause.

### Activation control at reaction endpoints and native equilibrium basis

Twelve saved-input reaction replays (four calls each for baseline reduced,
augmented temperature, and augmented native-clock controls) were linked with
the explicit native-activation variant. Baseline species-density difference
against the actual pinned ReactorCvode endpoint falls from about 2.74e-12 to
1.706e-13 kg/m3. Its upper-cell temperature difference falls from 1.009e-6 to
-1.643e-9 K. Augmented corrector endpoints differ by about 1.667e-13 kg/m3 and
-1.56e-9 K. Thus callback/temperature-cache chronology is sensitive to the
kinetics inputs; these data do not justify adopting augmented coordinates alone.
Commands, input SHA, executable SHA and results are in
`build/eb-activation-replay32-20261023/comparison.json`.

The remaining same-state rate difference was isolated separately: pinned
LiDryer mechanism.H uses refC=101325/8.31446/T (mol/m3), whereas the generic
Fortran equilibrium calculation uses the universal gas constant. An isolated
copy changes only that reference-concentration denominator to 8314.46 in
kmol units, alongside the generated activation variant. Across the same 16
states the largest mass-rate difference is then 4.761e-13 kg/m3/s and all
recovered temperatures still match exactly. No generated rates or canonical
mechanism were hand-edited. This diagnostic constant is specific to the
pinned header and must not become a generic chemistry default.

Evidence: `build/eb-same-state-rhs32-20261023/equilibrium_control.py`,
`equilibrium-control.json`, `equilibrium.csv`. Twelve additional reaction
replays with both controls leave endpoint differences near 1.706e-13 kg/m3
and 1.66e-9 K (`build/eb-native-kinetics-replay32-20261023/comparison.json`).
Matching same-state RHS therefore does not eliminate integration differences.
Next test this source-derived compatibility control in the bounded root32
trajectory with the original gates before deciding production policy.
Python syntax and project-contract checks pass; these are isolated diagnostic
builds, not a production implementation, C2 completion or physical validation.

### Bounded root32 native-kinetics trajectory control (pending)

A separate executable combines the existing ordinary probe and runtime with
only the generated native-activation mechanism and the isolated native-refC
kinetics object. It uses the exact original root32 import, one 100 ns step,
unchanged transport and solver settings. The native reference is reused;
no native solver rerun is necessary. All link-input hashes, binary/input
hashes, command, working directory and child PID are retained in
`build/eb-native-kinetics-trajectory32-20261023/run.json`. The driver enforces
a 1200-second bound; the earlier original run took 536.865 seconds.

The diagnostic comparator directly executes the authoritative checker's
loader, inventory and acceptance blocks without changing their thresholds.
It first checks that the saved original capture reproduces exactly the
original violation records, then can evaluate the new outputs. This does
not change the authoritative baseline receipt or qualify a changed mechanism.
Python syntax, canonical mechanism provenance and project-contract checks
pass. Trajectory and comparator replay were still running when this entry
was written; no new field-gate pass is claimed.

The comparator baseline replay subsequently completed with exit zero and
exactly reproduced all six original violation records. The trajectory is
still running (tool session 37698); compare its final CSVs after exit zero.

### Root32 native-kinetics control completed: trajectory error increases

The bounded run completed normally in 499.795292 seconds (binary SHA256
2d240133f55cf4d74ad4a04046220f803a25dd34dcda28cd2d7be8575886ff4b).
The original acceptance blocks still fail the same six fields. Maximum
momentum differences x/y/z are 5.84659e-6, 1.97807e-5 and 6.31907e-6
kg/m2/s, respectively 1.479, 1.489 and 1.694 times the canonical baseline.
Scaled density/total-energy/held-energy errors are 5.87696e-8,
2.76646e-7 and 2.76736e-7. The energy budget residual is 5.99204e-17;
H2-consumption relative difference is 6.68757e-7. Initial identity and
remaining original gates pass, but this is a failed trajectory comparison.

The y-momentum maximum moves from fine (25,32,26) to (25,32,30), and the
z-momentum maximum to (26,32,31). These still border the coarse/fine interface
and are above the embedded plane. Improving same-state chemistry fidelity
therefore does not improve this trajectory: chemistry changes also perturb
temperature-cache chronology and the shared interpolation limiter. The result
does not isolate either as the sole cause and does not support adoption of
the chemistry-only compatibility path.

All field errors and worst cells are in
`build/eb-native-kinetics-trajectory32-20261023/diagnostic-comparison.json`;
`summary.json` retains baseline/candidate ratios. Source and input identities
remain in `run.json`. Next use the new maxima and their coarse donors for
stage-pressure/limiter comparison, retaining original gates. Production code
is unchanged, and C2 convergence/trajectory qualification remains open.

### New-maximum donor capture launched

The failed native-kinetics control moves the y-momentum maximum to fine
(25,32,30). A bounded isolated prefix now records the coarse donor
(13,17,15), its six neighbors and the ghost (25,33,29), plus coarse-stage
source terms at z=15,16 and fine predictor halos around the new maxima.
It deliberately stops after the first fine predictor capture (exit77), with
a 600-second bound. A matching serial pinned-native donor trace records the
same donor, with a 180-second solve bound and exact exported-field neutrality
checks against the saved native run. Production sources are unchanged.

Artifacts: `build/eb-native-kinetics-donor32-20261023/` and
`build/pelec-new-donor32-20261023/`. Build/run records retain source and binary
hashes; handoff.json records sessions87514/87684 and the analysis command.
Both captures are pending. The reused limiter audit is prepared for the new
donor and all five Python scripts pass syntax checks; no capture result is
claimed yet.

The new-donor native run completed in 83.149351 seconds. Initial and final
65536-row exports have exactly the original hashes; diagnostic insertion is
neutral for this capture. Native donor (13,17,15) events2--6 have an x-axis
shared slope scale of 0.9019448035085728, limited by mz; z slopes are limited
to zero by held energy. The recorded native events still need matching to
PeleF interpolation time before a difference is attributed. Native-only
analysis is `build/eb-native-kinetics-donor32-20261023/native-limiter.json`.
PeleF prefix session87514 remains running; the two-solver comparison is pending.

### New-maximum capture completed: shared-limiter amplification persists

PeleF stopped at the intended prefix marker (exit77) in 284.147125 seconds.
At interpolation fraction0.5 its duplicate donor records agree exactly.
Native events2/3/4 are identical and match the midpoint much more closely
than the initial/final events (maximum scaled center difference1.826e-9,
versus6.937e-5). The audited interpolation with offset(-.25,-.25,-.25)
reconstructs the captured PeleF ghost gas at (25,33,29) exactly.

At donor(13,17,15), the x shared scale is0.882320524528229 for PeleF and
0.9019448035085728 for native; mz limits both. The donor-center mz difference
is-1.8257944774e-9 kg/m2/s. Offline substitution of only mz at the center and
six neighbors changes the target ghost my difference from6.0360128455e-5
to2.2400914457e-11 kg/m2/s. Substituting every other component leaves
6.0360105996e-5. This isolates the same cross-component amplification at the
new maximum, despite the improved same-state chemistry evaluator. It is not
a production modification or a proof that altering the limiter is appropriate.

Evidence: `matched-limiter.json`, `input-ablation.json`, and recorded hashes
in `build/eb-native-kinetics-donor32-20261023/handoff.json`. Both jobs have
finished. Next compare coarse-stage pressure and directional flux at
(13,17,15) and (13,17,16), upstream of this momentum discrepancy. Python
syntax, diff whitespace and project-contract checks pass. C2 remains open.

### New-donor coarse pressure/flux capture launched

An isolated coarse-stage probe now records pressure, reconstructed face
states, directional hydro fluxes and surrounding 5x5x5 states at coarse
(13,17,15) and (13,17,16), using the same native-kinetics diagnostic objects.
It stops after coarse stage2 (expected exit77;400-second bound). The matching
pinned-native instrumentation records transport and total face fluxes plus
Riemann inputs at the same cells; normal exported trajectories must remain
identical. Native building precedes its bounded run. No production changes.

Artifacts: `build/eb-native-kinetics-faces32-20261023/` (session7247) and
`build/pelec-new-coarse-riemann32-20261023/` (session40221; intermediate flux
build in `build/pelec-new-coarse-flux32-20261023/`). Prepared analyses verify
PeleF initial/center-state neutrality against the completed donor capture,
then compare pressure and flux and independently normalize native tile events.
Python syntax passes; captures and numerical comparisons remain pending.

The native new-coarse-flux run completed in84.175670 seconds with initial and
final65536-row exports exactly identical to the saved reference. The144
Riemann records normalize to36 unique stage/probe/face rows; all duplicate
values agree. Their momentum flux agrees with the independently recorded
transport-subtracted total face flux to4.500e-18 Pa. Native normalization is
complete; PeleF session7247 remains running, so cross-solver pressure and
flux differences are not yet available.

### New-donor pressure/flux attribution completed

PeleF completed the intended coarse-stage capture (exit77) in164.580189 s.
Initial exports and12 repeated coarse-center states exactly match the prior
donor capture; all36 selected faces have unit aperture and full-volume cells.
Maximum hydro-flux differences over stages0/1/2 are1.528e-9,2.731e-6 and
5.131e-5 Pa. Corresponding maximum cell-pressure differences are2.911e-11,
2.372e-6 and4.596e-5 Pa, so the discrepancy is already present before face
reconstruction at the final coarse evaluation.

Stage2 mz RHS differences at coarse z=15/16 are-0.0729891/+0.0973426
kg/m2/s2. The z-face contributions alone are-0.0729887/+0.0973423;
x/y contributions are at most3.587e-7. Transport differences are
-4.272e-5/+6.827e-5 kg/m2/s2. Combining old and final flux RHS with the
actual100 ns trapezoidal update and midpoint interpolation reproduces the
observed donor mz differences (-1.82579e-9/+2.43527e-9 kg/m2/s), with
residuals1.07605e-15/-1.94750e-16. Independently reconstructed PeleF flux
RHS matches the recorded source to2.843e-14 kg/m2/s2.

Evidence: `build/eb-native-kinetics-faces32-20261023/flux-comparison.json`,
`flux-summary.json`, `directional-attribution.json`, and the native
`riemann-comparison.json` in `build/pelec-new-coarse-riemann32-20261023/`.
This closes the momentum-source balance for these two diagnostic cells, not
C2 qualification. Next replay temperature recovery on the saved final-stage
states and guesses; the limiter remains unchanged. Python syntax, whitespace
and project-contract checks pass. No jobs from this capture remain running.

### New-donor saved-state EOS replay and bounded temperature-coordinate test

Thirty native EOS evaluations cover the two cells, three coarse stages and
five temperature guesses per state. For stage2, replay with the saved guess
reproduces PeleF pressure within2.911e-11 Pa. First Newton corrections at
z=15/16 are-9.01690e-7/-8.26340e-7 K, smaller than the pinned1e-6 K
stop-before-update threshold. The pressure differences against native are
4.59598e-5/4.23949e-5 Pa. Offsetting the same state's initial guess by
+/-0.001 K reduces those differences to approximately-6.151e-7/-2.877e-7 Pa;
+/-1 K gives-6.177e-7/-2.904e-7 Pa. Native states/guesses at these new cells
were not captured, so this establishes saved-guess sensitivity, not an exact
native/F-state swap or full attribution of the residual.

Artifacts: `build/eb-new-donor-eos32-20261023/audit.py`, exact inputs, outputs,
and `audit.json` with binary/input hashes. The canonical EOS is unchanged.

A bounded follow-up links the previously tested all-species-plus-temperature
CVODE prototype under the normal reactor module name into the same isolated
coarse-face probe, retaining the explicit native-kinetics controls. It uses
native CV_NORMAL clock behavior and finite-difference Jacobian, as in the
single-cell prototype; these are coupled changes and are not claimed to
isolate temperature coordinates alone. Only the coarse prefix is run
(expected exit77;400-second limit), to test pressure improvement before any
full trajectory or production decision. Session42715 and records are in
`build/eb-temperature-kinetics-faces32-20261023/`. No production source edits.

### Integrated-temperature coarse-prefix result: pressure improves

The candidate reached the intended coarse stop (exit77) in173.098968 seconds.
Stage2 maximum cell-pressure difference falls from4.59597e-5 to6.12345e-7 Pa;
maximum face hydro-flux difference falls from5.13051e-5 to7.31699e-7 Pa.
Stage1 differences remain approximately2.37413e-6 Pa in cell pressure and
2.73307e-6 Pa in hydro flux. This is improvement at the diagnosed final
coarse stage, not complete stage agreement or whole-trajectory qualification.

The first initial-CSV byte-identity assertion failed because this candidate
changes the independently computed initial reaction history. Inspection of
all65536 rows showed exact gas, geometry and temperature identity; budget CSV
is also exact. Only rho_omega fields and rhoe_dot differ (maximum species
source difference2.42508e-6 kg/m3/s and energy-source difference5.82077e-4
J/m3/s). The diagnostic comparator now checks gas/geometry/temperature/budget
identity explicitly and records every history difference in initial-identity.json,
consistent with the original reference checker's separate history reporting.
No original field acceptance threshold has changed.

Evidence: `build/eb-temperature-kinetics-faces32-20261023/flux-summary.json`,
`flux-comparison.json`, `initial-identity.json` and build/run records.
The candidate combines integrated temperature, all species, finite-difference
Jacobian and native-clock handling with native-kinetics controls; their
individual effects are not separated by this result.

A full root32 one-step100 ns candidate run was launched (session77873),
using the same import and saved native reference, with a1200-second limit.
`build/eb-temperature-kinetics-trajectory32-20261023/` records every link input,
binary/input SHA and PID; compare.py applies unchanged original gates after
normal completion. Syntax, whitespace and project-contract checks pass.
Production sources remain unchanged; C2 is still open.

### Integrated-temperature full trajectory: rejected MOL step

The full candidate solver exited with code 1 after 613.542532 seconds:
`ERROR STOP Reference-state MOL step rejected` at probe line 329. The Python
driver returned zero, so its process exit is not evidence of solver success;
`run.json` records the authoritative solver failure. No `fortran_001.csv`
was produced. The attempted comparator stopped with FileNotFoundError,
therefore no final-field, inventory, or budget acceptance result is available.
Binary SHA256: `0f69d0aad124bc188e784232ff09cb42ab796a6a120d9cd16cea8e5fa2a965f1`.
The coarse-prefix pressure improvement does not establish a usable full step.

An isolated rejection diagnostic in
`build/eb-temperature-rejection32-20261023/` preserves this candidate's
numerical operations and adds logging to rejection branches of
`advance_forced_reactor_cvode`. It records source-line identifiers, CVODE
status, trial coordinates, and original replay inputs. Session 26120 ran the
same root32 input with a 1200-second bound; no production source or tolerance
was changed. If no reactor rejection is logged, subsequent MOL admission
remains the next localization target. Compiler warnings concern existing real
equality and unused arguments; compilation and Python syntax checks passed.
The fixed H2/O2 mechanism provenance preflight passed. These checks do not
qualify the diagnostic numerical candidate or close C2.


The reactor diagnostic completed with solver code 1 and zero instrumented
reactor rejection events. The subsequent MOL return diagnostic
`build/eb-mol-rejection32-20261023/` completed in 605.990426 seconds, also
with solver code 1. Its 67,688 return events include normal early returns;
only the failure propagation is used here. Fourteen attempts returned from
`recover` at original tree line 11382 after `reset_held_eb_tree_family_3d`
reported failure, then from the predictor at line 8035. Final propagation
was through lines 6561, 6745 and 4761. The evidence localizes admission
failure to held-energy reset/recovery; it does not yet identify the failing
cell or establish that the reset policy itself is defective.
`return-analysis.json` verifies the original tree source SHA before mapping
logged line numbers. Diagnostic binary SHA256:
`c6738f363ec4775a0eed21d4e477a67ad69a70b67d8aed019201a59f5521423b`.

The next isolated diagnostic, `build/eb-held-rejection32-20261023/`, saves
state, held energy and temperature before the cell reset and stops at its
first failed admission with diagnostic exit 77. This is intended to enable
single-cell replay without repeating the full retry sequence. Production
sources and acceptance thresholds remain unchanged.


The held diagnostic terminated at its intended exit 77. The first rejected
cell is node 2, local (8,16,22), selected depth 1, at original tree line
2988 (`reset_nasa7_dual_energy_cell`). All captured values are finite and
all species densities positive (minimum 2.0042658941e-9). The relative
species-sum minus density residual is -1.2925524828069153e-11, exceeding
the unchanged 5e-12 reset admission bound. The standalone reset replay
reproduced rejection and exact preservation of state, held energy and
energy-change output. Two separate saved-state controls, setting density
to the species sum or closing only the largest species to density, both
produce zero residual and successful reset. These controls isolate mass
closure as a sufficient cause of this rejection; neither is an adopted
repair. Evidence: `build/eb-held-rejection32-20261023/rejection.json`,
`reset-controls.json`, and the associated replay inputs/logs. Log SHA256:
`e60e0ad633b24f0bf4b1d5e2bc4e7ad2ea1022db7e2168c94fe06c0712cf515f`.

Next trace the all-species reactor endpoint and reaction-source sum into
the fine predictor. `react_sdc_state_3d` sets endpoint density from the
species sum but returns species reaction sources independently; the
predictor's mass closure must be checked before choosing a correction.
The temperature-coordinate candidate still has no accepted full root32
step. No production tolerance or conservation requirement was relaxed.


The predictor mass diagnostic (`build/eb-predictor-mass32-20261024/`)
completed with intended diagnostic exit 77 in 403.762777 seconds. Its
last candidate vector exactly equals the previous rejected saved state.
Two fine predictor events use dt=50 ns. In both, the old-state and hydro
species-sum minus density residuals are zero. The combined-rate dt
contributions are -1.4089447831e-15 and -1.2518563018e-12 kg/m3; the
second candidate residual is -1.2518597270e-12 kg/m3. Reconstructing
all 14 candidate components from old+dt*(rates+hydro) gives zero error.
This isolates the incoming combined rate, not yet its reaction versus
spatial constituents. `combined_rates` adds IR to species without a mass
source; the next diagnostic must separate those constituents and capture
the preceding reactor input before assigning the defect to integration.
Evidence: `mass-analysis.json`; diagnostic binary SHA256
`a4c7551157b83966f03bf16ccfaf915e101b9940ddf763fc43e4eda3a65ffb8a`.

### Local reactor mass-invariant isolation (2026-10-24)

The constituent diagnostic completed with intended solver exit 77 in
407.943701 seconds (`build/eb-mass-reactor32-20261024/run.json`). Both
captured predictor events have zero E0/E1/D0/D1 mass residual, and neither
optional auxiliary source nor admitted mean rate is active. IR sums are
-2.8178894152e-8 and -2.5037126038e-5 kg/m3/s. Three saved reactor calls
at node 2, local (8,16,22), reproduce the last IR sum in standalone replay
(-2.5037126038e-5). This is local reproduction, not a completed root32 step.

`build/eb-reactor-invariant32-20261024/` instruments the isolated
all-species plus temperature prototype before endpoint clipping. The three
raw endpoint mass errors are -1.41553e-15, -1.95455e-13 and -1.25186e-12
kg/m3. Across 61, 112 and 61 RHS calls (40, 61 and 40 CVODE steps), the
largest absolute chemical mass-rate sum is 4.39926e-15 kg/m3/s. Thus the
large residual is already present at the integration endpoint; it is not
created by SDC source unpacking or endpoint clipping. This evidence does
not yet distinguish nonlinear convergence, finite-difference Jacobian,
or dense-output contributions within CVODE.

An isolated prototype in `build/eb-reduced-temperature32-20261024/`
combines eight independent species densities with an integrated temperature,
reconstructing the initially largest species from the analytic forced mass
balance at every RHS evaluation. It keeps the same tolerances and finite-
difference Jacobian. All three local calls finish successfully; endpoint
mass errors are 1.38778e-17 kg/m3, with IR sums between 1.91e-10 and
3.48e-10 kg/m3/s. Relative to the all-species candidate, the largest endpoint
temperature change is 2.32149e-10 K. These are same-input candidate
comparisons, not independent pinned-native agreement or element-conservation
qualification. No production source or admission threshold was changed.

Source base: main `0a0e00ebe9138145fbb0e0e745665d0d768df009` with the
existing dirty worktree. Each prototype `record.json` records exact compile
and link commands, input/source/binary/log hashes and replay exit 0; the
compiler is `/usr/bin/gfortran`, using the existing isolated OpenMPI selected
runtime and pinned SUNDIALS 7.2 static libraries. This turn independently
passed H2/O2 provenance preflight, Python syntax for the three new scripts,
and project contract 0.246.0. No MPI or full application test was rerun.
Next compare the reduced-temperature prototype with pinned-native local
reactor results and the saved coarse-stage pressure evidence before a full
trajectory or production adoption. C2 remains open.

The reduced-temperature prototype was subsequently compared with the pinned
native `ReactorCvode` replay for all three fine-cell inputs. Six fresh native
runs (analytic and finite-difference dense Jacobians for each input) exited
successfully, using the captured relative tolerance and species-density
absolute tolerance converted to CGS. Maximum temperature and species-density
differences are 2.41812e-9 K and 2.47670e-13 kg/m3. The native executable
hash was checked against its existing build receipt before execution.
`native-comparison.json` records the outputs, tolerances and input/binary
hashes in `build/eb-reduced-temperature32-20261024/`.

Four saved coarse-cell inputs also completed with the new prototype. Against
the previously recorded native controls at matched tolerance, the largest
temperature difference is 1.69030e-9 K and species-density difference is
1.74862e-13 kg/m3 (`coarse-comparison.json`). Those coarse native controls
were reused, not independently rerun this turn. These local comparisons
support proceeding to the stage pressure/flux diagnostic; they do not prove
root32 trajectory, MPI lifecycle, or physical convergence. A bounded coarse
prefix diagnostic is in progress; its result must be read from `record.json`
and checked with `compare-flux.py` before it is reported as captured/passed.

The reduced-temperature coarse prefix subsequently completed in 167.249208
seconds with intended diagnostic exit 77. The original comparison checks
initial gas/geometry/temperature identity and budget-byte identity before
comparing the pinned native face trace. Stage 0/1/2 maximum cell-pressure
differences are 2.91038e-11, 2.37356e-6 and 6.12825e-7 Pa; maximum momentum
hydrodynamic face-flux differences are 1.52795e-9, 2.73250e-6 and 6.79167e-7
Pa. The prior all-species temperature candidate's stage-2 pressure/flux
values were 6.12345e-7/7.31699e-7 Pa. Thus restoring mass closure preserves
the observed coarse pressure improvement, without establishing fine-step
or complete-trajectory agreement. Candidate binary SHA256:
`a55d22b50821eaaf00207baf1e77a74a7199542aabe58831034a6272cc2c0214`.

A bounded original root32 100 ns trajectory is now running under
`build/eb-reduced-temperature-trajectory32-20261024/` (1200-second cap).
Its comparison script requires terminal solver success and unchanged binary
and input hashes before applying the original field and budget acceptance
blocks. No full-trajectory result is claimed here. The prototype remains
isolated from production; C2 and physical qualification remain open.

While the trajectory runs, an independent inventory calculation checked H,
O and N for all seven saved local inputs (three fine, four coarse). It uses
the mechanism's atomic composition and molecular weights to sum
`atomic_count/MW * (endpoint_rhoY - (old_rhoY + dt*forcing))`, normalized by
the initial inventory of the same element. The reduced-temperature candidate's
largest absolute relative residual is 1.89659e-14 (H, first coarse input).
The corresponding maxima over saved native analytic/FD outputs are
2.27095e-14 and 2.97856e-14. These are bounded local inventory measurements,
not a newly defined acceptance gate or whole-trajectory conservation verdict.
`build/eb-reduced-temperature32-20261024/element_balance.py` and
`element-balance.json` retain the definition, mechanism hash and all 63
per-element measurements. The coarse native outputs remain historical controls.

The isolated reduced-temperature reactor also passed a forced-failure/context
reuse check while the trajectory was pending. Each of the three captured fine
inputs was limited to one CVODE step, repeated 65 times (exceeding the 64-slot
context registry capacity if slots leaked). All 195 calls rejected with zero
completed steps and exactly unchanged species densities and temperature.
An ordinary call after each sequence succeeded and reproduced the successful
SDC replay species endpoint exactly. This is a bounded transactional/reuse
check, not an MPI test or full generic-mechanism regression. Commands, binary
and source hashes are in `build/eb-reduced-temperature32-20261024/rollback.json`;
expected CVODE failure messages are retained in `rollback.stderr`.

The reduced-temperature root32 100 ns trajectory finished with solver exit 0
in 526.751158 seconds. It avoids the preceding all-species prototype's mass-
closure rejection, but **fails seven original reference field gates**:
density 3.13670e-6 (limit 1e-8), x/y/z momentum 5.46901e-4,
2.85604e-4 and 1.19453e-4 (each limit 1e-6), scaled total/held energy
4.31533e-6/4.31440e-6 (each limit 1e-8), and temperature 0.00180884 K
(limit 0.001 K). The authoritative checker's correction-accounted energy
budget relative residual is 5.73512e-17. Successful execution and bounded
budget closure therefore do not establish reference agreement.

The candidate y-momentum maximum is approximately 21.5 times the original
baseline value 1.32832e-5. Localization places it at fine cell (32,1,28),
with 98.5666% of volume-weighted squared y-momentum error in the fine
boundary layer and only 0.02174% in cut cells. The density/temperature maxima
are at fine (32,7,27); x/z momentum maxima are at (32,8,27). All seven field
maxima are in full-volume cells. The prior two-probe coarse pressure/flux
improvement does not cover this newly dominant upper-x interface region.
The candidate remains unadopted. Next inspect the coarse donors and capture
fine halo/face pressure and momentum stages at these cells before choosing
another numerical change or launching another full trajectory.

Evidence: `build/eb-reduced-temperature-trajectory32-20261024/` contains
`run.json`, `diagnostic-comparison.json`, `localization32.json`, `summary.json`
and the completed CSVs. The comparator checks unchanged binary/input hashes
and uses the original field/budget acceptance blocks; localization verifies
that its extrema exactly reproduce the comparison report. Binary SHA256:
`d708e76fb2928753007149b68b46585f1473180e927ea5051aed4de12ccc43aa`.
No MPI, broad configuration matrix, installation or physical qualification
is claimed. Production source and acceptance limits remain unchanged.

The upper-x prefix diagnostic completed for baseline and reduced-temperature
candidate in 368.440507 and 355.061753 seconds, respectively. Both intentionally
stopped with exit 77 after the first fine stage 2 (`UPPERX_PREFIX_COMPLETE`).
`build/eb-upperx-stage32-20261024/` retains build/run identities, pressure
stencils, stage-face CSVs, `comparison.json`, and `origin.json`. These are
candidate-minus-baseline measurements, not a new pinned-native comparison.

At fine probe 20, (32,1,28), stage 1 has exactly equal cell momentum in both
runs, but the candidate's y-momentum spatial rate differs by -4168.80615.
The x/y/z directional contributions to that rate difference are -0.885398,
-4167.92076 and -6.64409e-11. The low-y face raw y-momentum flux differs by
-0.488830023 Pa, whereas its reconstructed left pressure differs by only
-0.000277123007 Pa. Across the captured fine stage-1 faces the largest right
pressure difference is 0.000483880052 Pa. Halo-inclusive pressure stencils
have a larger maximum difference of 0.00245618305 Pa, at (34,6,26).
Stage 2 then has a maximum cell y-momentum difference of -0.000120788161
and a maximum sampled cell pressure difference of -0.276288719 Pa.
This establishes that a large normal momentum-flux difference precedes the
observed cell-momentum difference; pressure differences alone do not identify
its cause. Density, normal velocity, sound speed and reconstruction inputs
are not captured by this pressure-only halo diagnostic.

For all 18 captured cell/stage samples per variant, reconstruction of all
three momentum rates from aperture-weighted directional centroid fluxes
passes the roundoff-scaled check (54 component checks per variant, 108 total).
At these sampled full-volume faces raw and centroid y-momentum fluxes are
identical; interval-normalized tree flux differs by at most 1.45519e-11 Pa.
Thus the stage-face representation does not introduce the observed difference.
This is not a test of later reflux/synchronization: capture ends before that
fine-step completion. The next discriminating diagnostic is the full halo
primitive and reconstructed Riemann states at the low-y face of (32,1,28)
and upper-x faces of (32,7,27)/(32,8,27), before changing interpolation or
flux synchronization. No production numerical change is justified yet.

The full Riemann-input capture in `build/eb-upperx-riemann32-20261024/`
completed baseline/candidate prefixes in 308.072594/301.906466 seconds,
each intentionally stopping with exit 77 at fine stage 1. All 90 stage-face
rows per variant exactly match the preceding diagnostic prefix, so adding
input/reconstruction capture did not alter those recorded results. Independently
reconstructed face pressures also matched the actual solver pressure pair
exactly at runtime. `comparison.json` records input hashes and all deltas.

At fine stage 1, valid-cell stencil velocities agree exactly, while ghost
velocities already differ. At upper-x faces (lower cell (32,7,27) and
(32,8,27)), ghost velocity differences are (0.0174077349,-0.145197701,
8.14181e-7) m/s. At the low-y face (lower cell (32,0,28)), ghost differences
are (0.0110855273,-0.0110980582,4.01291e-7) m/s. This locates the velocity
difference upstream of face reconstruction, rather than creating it solely
inside the Riemann calculation. The exact halo interpolation source remains
to be isolated; these captures do not yet identify which limiter or donor
operation introduces it.

A local production-routine replay (`replay.F90`, `replay.py`, `replay.json`)
reproduces all three momentum components on all three faces for both variants
exactly with the ordinary Riemann branch. The alternate Lax-Friedrichs branch
has maximum residual 35.8106468 Pa and is not the branch of these captures.
Replacing only candidate stencil velocities with baseline velocities changes
the low-y normal momentum-flux difference from -0.488830023 to -0.000173667373
Pa (approximately a factor of 2815 smaller). Replacing only pressure leaves
-0.488656377 Pa; replacing only cell sound speeds leaves -0.488830023 Pa.
At upper-x (32,7,27), velocity replacement changes the normal flux difference
from -0.764033551 to 0.000249641555 Pa; pressure replacement leaves
-0.764283224 Pa. Thus velocity input differences dominate these local face
flux changes. These controlled mixed-input replays diagnose local sensitivity,
not a conservative or qualified trajectory correction. Next trace the ghost
velocity interpolation/limiting and its coarse donor inputs. No production
numerical change or C2 qualification is claimed.

The coarse-donor capture in `build/eb-upperx-donor32-20261024/` completed
both intentional exit-77 prefixes (303.218702/297.339396 seconds). Each variant
contains 12 interpolation events and preserves the preceding stage-flux CSV
byte for byte. `analyze.py` reproduces all captured common-limiter results
within 128 epsilon scaled by component magnitude. `match.py` matches the
three actual fine stage-1 ghost primitive inputs within 1e-12; each has two
identical repeated donor records, retained without assigning a unique call.
The captured policy is the default shared limiter, not separate momentum.

At both upper-x ghosts, the z-axis shared scale changes from 1 to 0, limited
by H2O, H, O and OH. At the low-y ghost it changes from 1 to
0.505412358950422, limited by O. The x/y scales remain zero. Offline
`ablate.py` substitutes baseline species donor inputs into the candidate
while preserving its density, momentum and energy inputs. This restores the
z scale to 1 and reduces upper-x velocity differences from
(0.017407735,-0.145197701,8.14181e-7) m/s to at most 9.56e-11 m/s;
the low-y maximum residual becomes 1.60e-10 m/s. Replacing just H2O/H/O/OH
has the same effect. Replacing O alone does not restore the upper-x limiter,
and restores the low-y scale only to 0.9210803694992324. Maximum neighboring
species-input differences are 1.218e-12 and 1.224e-12, respectively.

This identifies amplification of small species donor differences by the
shared limiter as the source of the sampled candidate-versus-baseline ghost
velocity differences. Together with the preceding local Riemann replay,
it explains the sampled normal momentum-flux change before fine-step reflux.
It does not establish which chemical integration trajectory matches pinned
PeleC, nor justify removing species from its reference interpolation policy.
The next reference comparison must inspect the same native coarse donors
and limiter factors; a policy change alone would not prove C2 parity.
No production numerical change, new MPI evidence or physical qualification
is claimed. Artifacts: `run.json`, `limiter.json`, `matched.json`,
`ablation.json`; source HEAD remains
`0a0e00ebe9138145fbb0e0e745665d0d768df009` with existing working-tree changes.

The matched native donor probe in `build/pelec-upperx-donor32-20261024/`
completed with exit 0 in 83.886919 seconds. The initial and 100 ns final
exports each contain 65,536 rows and match the saved pinned reference exactly.
The diagnostic binary SHA256 is
`2d356ece57771f7b637bda348c45f48baf20cd794e1fc43cdf22f88bef5845c0`.
Its AMReX interpolation source matches the preceding native donor capture.
The probe includes upper-x parent (17,4,14) and both periodic representations
of the low-y parent; captured low-y records use the outside-domain image.

`comparison.json` retains all native candidates ranked by scaled donor-center
difference. Events 2/3/4 tie and have identical centers and neighbors. Their
center differences from baseline are 4.385e-11 at upper x and 6.960e-11 at
low y; their z shared-limiter scale is 1. The next distinct events differ
by approximately 6.237e-5 and 3.444e-5. This strongly separates the nearby
state candidates but is not an independently timestamped stage assignment.

`gradients.json` records the trace-species differences in those nearest native
donors. At upper x, native H2O z-minus/plus differences from center are
(-3.615e-15,+1.616e-15), whereas the reduced-temperature candidate has
(-2.270e-14,-1.889e-14). H, O and OH also change to same-sign differences
in the candidate. Thus the candidate develops trace-species extrema that
zero the shared slope. This is consistent with the offline attribution,
and favors diagnosing chemical endpoint/spatial differences before altering
reference interpolation. An exact stage/time correspondence and the origin
of those small donor differences remain open. No production correction or
C2 trajectory qualification follows from these local observations.

The timestamped native interpolation probe in
`build/pelec-upperx-time32-20261024/` completed exit 0 in 83.370170 seconds.
Initial/final exports remain identical to the pinned reference, and the entire
native donor CSV is byte-identical to the preceding capture. The diagnostic
binary SHA256 is
`5bbbf242737d2402c1f7b4e6c8c9289856e46d348a30f9f57a851c2bf598fcb0`.
Instrumentation brackets the three MOL FillPatcher calls and resets context
after each call. Events 2/3/4 all request 50 ns at level 1 with a 50 ns step:
event 2 is first fine-step stage 1, event 3 is its MOL iteration 2, and event 4
is the next fine-step stage 0. Thus event 2 provides the relevant timestamped
native comparison for the first fine stage-1 ghost; the repeated identical
records reflect reuse at the same coarse interpolation time, rather than
three unresolved physical times. `time-assignment.json` preserves all events.

`endpoint_check.py` also checks whether 50 ns donors can be reconstructed as
half of the saved initial plus final coarse plotfile fields. This does not
hold at the sampled interface parents: maximum neighboring species residuals
are 3.4793e-6 and 4.6250e-6 kg/m3. Selected farther neighbors agree to
6.94e-18. The final export is after coarse/fine synchronization and is not an
authoritative replacement for the pre-synchronization coarse endpoint used
by FillPatcher. This check rules out using final saved plotfiles for the next
chemical endpoint attribution; it does not identify a chemistry defect.
Next capture the pre-synchronization coarse reaction inputs/endpoints at
these parents and z neighbors. Existing C2 gates and limiter policy remain
unchanged.

The six-cell native coarse reaction capture in
`build/pelec-upperx-react32-20261024/` completed exit 0 in 84.656118 seconds;
initial/final fields and prior donor trace retain their recorded hashes.
It yields 18 before/after pairs: six initialization records and two subsequent
coarse reaction calls for each parent and its two z neighbors. Diagnostic
binary SHA256 is
`9317a002ec20702d5c9cdd007f4f81d3687d3c950d67a05b657f7aa4824fbd67`.
`endpoint_match.py` reconstructs all nine species of the event-2 50 ns donors
from initial species and the final coarse reaction endpoint (reaction event 4)
with exactly zero residual at all six cells. This identifies the relevant
pre-synchronization chemical endpoint without using the final plotfile.

A locally relinked reduced-temperature candidate replay uses the native
physical old states, nonreacting sources and held energy, with CGS-to-SI
conversion and the native density-scaled absolute tolerance. All 12
non-initial reaction records complete; maximum species difference is
1.7160e-13 kg/m3 and maximum temperature difference 1.694e-9 K.
The comparison preserves conversion-induced energy-source residuals; it is
not a claim of bit-identical arithmetic across the two unit systems.
Despite small endpoint differences, the replayed 50 ns O gradients already
form extrema: upper-x minus/plus differences are
(-6.064e-15,-1.165e-15) kg/m3, and low-y differences are
(+6.189e-16,+7.075e-15) kg/m3. Native gradients at these points straddle zero.
Thus even replaying native physical inputs through this candidate can alter
the limiter branch; upstream flow-input discrepancies are not necessary for
this local effect. The integrator/Jacobian/roundoff route needs isolation
before any production adoption. No tolerance relaxation or limiter change
is made; C2 remains open. Scripts/artifacts: `analyze.py`, `reaction-inputs.json`,
`replay.py`, `replay-build.json`, `replay-comparison.json`, `endpoint-match.json`.

The same 12 captured native reaction inputs now have isolated native controls
(`native_controls.py`, `native-controls.json`). Native dense analytic-Jacobian
replay reproduces all recorded species and temperatures exactly. Native dense
finite-difference Jacobian replay differs by at most 3.9961e-14 kg/m3 and
6.658e-10 K, but all four sampled trace-species z gradients retain opposite
signs around both parent centers. Thus replacing an analytic Jacobian with
finite differences alone does not reproduce the candidate's local extrema.
These are native controls, not a test of the Fortran candidate's Jacobian.

`roundtrip.py` tests native analytic-Jacobian inputs after CGS/SI/CGS roundtrip:
endpoint species remain exactly equal. Using the energy forcing reconstructed
through the Fortran wrapper's SI arithmetic instead gives maximum species
difference 1.356e-17 kg/m3 and preserves the O gradient signs at both parents.
The reconstructed energy-source residual is at most 5.239e-4 W/m3, but this
control does not create the observed limiter-branch change. These bounded
controls direct the next comparison toward the candidate's reduced state
coordinates, RHS/temperature evaluation and CVODE trajectory rather than
attributing the extrema to unit conversion or Jacobian choice alone.
No tolerance, golden reference or production numerical policy is changed.

The full-species temperature-coordinate replay (`full_state.py`) completes
all 12 native-input cases with maximum species difference 1.72291e-13 kg/m3.
It retains an O extremum at the low-y parent, although the upper-x O
minus/plus differences straddle zero. Reduced coordinates therefore affect
the local gradient but are not necessary for every observed extremum.

Code inspection also found that the prior temperature-coordinate diagnostic
still called the energy-inversion RHS using integrated temperature as its
guess. `direct_temperature.py` builds isolated full/reduced variants that
evaluate production rates and thermal capacity directly at the integrated
temperature, retaining initial/final thermodynamic admission. Both complete
12/12 cases. Maximum species errors remain 1.72291e-13 and 1.71598e-13 kg/m3,
respectively; this change does not remove the low-y O extremum. These local
controls do not qualify either integrator or justify production adoption.
They narrow the remaining investigation to numerical scaling/error control,
initial-temperature handling and RHS arithmetic beyond this temperature
inversion choice. Compile/link commands and source/binary hashes are stored
in `direct-temperature.json`; original tolerances remain unchanged.

CVODE error-weight controls (`weights.py`, `weights.json`) compare the full
integrated-temperature diagnostic with native physical tolerances. The native
reactor sets vector absolute tolerances to 1e-18 in its CGS species-density
and Kelvin coordinates (without typical-value scaling in these replays).
The Fortran diagnostic previously applied 1e-15 to both SI density and Kelvin.
A vector tolerance keeps species at 1e-15 kg/m3 and sets temperature to
1e-18 K; separate variants retain the supplied initial temperature and combine
both changes. All three complete 12/12 cases, with maximum species difference
1.72291e-13 kg/m3 and the same low-y O extremum. This corrects a diagnostic
weight mismatch but does not explain the local branch change.

`cgs_coordinates.py` additionally scales internal full-species CVODE
coordinates and derivatives to CGS while retaining SI rate evaluation and
physical endpoint admission. It completes 12/12 cases with maximum species
difference 1.72468e-13 kg/m3; low-y O minus/plus gradients remain same-sign
(+1.34133e-15,+6.85880e-15) kg/m3. The sampled extremum therefore persists
under these weight, initial-temperature and coordinate-scaling controls.
The next useful comparison is the actual same-state species and thermal RHS
and CVODE configuration, not another whole-flow run. These are isolated
build-directory diagnostics; no production solver or acceptance tolerance
has changed.

Inspection identified a substantive CVODE configuration difference: pinned
PelePhysics `ReactorCvode.H` defaults `m_cvode_maxorder` to 2 and explicitly
sets it, while the candidate diagnostic omitted this call. Isolated maximum-
order-2 controls (`order2.py`, `order2.json`) complete 12/12 native-input cases
for full, reduced and CGS full coordinates. All restore opposite-sign O
z-neighbor differences at both parents. Maximum species differences become
1.65146e-15, 8.95200e-15 and 2.44249e-15 kg/m3, respectively.

A separate control changes only the maximum order in the original reduced-
temperature candidate (`original_order2.py`, `original-order2.json`), without
the direct-temperature, vector-weight or CGS alterations. It completes 12/12
and reduces maximum species error to 8.95201e-15 kg/m3; upper-x O gradients
are (-2.21817e-15,+9.92822e-16), low-y (-3.07383e-15,+3.38252e-15) kg/m3.
This supplies a concrete configuration cause for the sampled branch change.
A same-condition 32-cell/100-ns trajectory test is launched under
`build/eb-order2-trajectory32-20261024/`; its results are pending and no
trajectory gate is claimed passed here. Production source remains unchanged.

### 2026-10-24: maximum-order cost and optimized diagnostic trajectory

The original order-2 32-grid trajectory above terminated at its 1200-second
limit (1200.0705 s). Its wrapper killed and reaped its own solver; only the
initial field and budget were written, with no restart artifact. This is an
incomplete run, not a failed numerical acceptance check. See
`build/eb-order2-trajectory32-20261024/timeout-analysis.json` and `run.json`.

Under `build/pelec-upperx-react32-20261024/`, `profile-order.json`,
`profile-full-order.json` and `profile-native.json` record the same 12 captured
reaction inputs. Reduced-temperature default/order-2 internal step totals are
587/5543; full-species direct-temperature order-2 uses 5404 steps and 5608 RHS
evaluations. Pinned native analytic-Jacobian replay uses 5400/5604; native
finite-difference replay uses 5394/5599 plus 600 linear-solver RHS evaluations.
All order-2 controls have zero error-test failures. Native verbose replay
outputs are byte-identical to the corresponding earlier native controls.
These counters cover only the captured cells, not the full trajectory.

The 240-case replay (20 repetitions) costs 8.190 s for the original reduced
order-2 candidate, 5.832 s for full-species direct temperature, and 0.165 s for
native. This is not a matched release-performance comparison: native uses
O3, whereas the Fortran libraries were built for diagnostics. Optimizing only
the forced reactor wrapper leaves 5.832 s; also optimizing kinetics gives
3.873 s, and NASA7 thermo gives 2.170 s. These isolated relinks preserve all
240 numeric output rows exactly. Commands and source/binary hashes are in
`benchmark_replay/`, `optimized_full_wrapper/` and `optimized_dependencies/`.

The isolated kinetics diagnostic then removes unused concentration-derivative
work and reuses Gibbs thermo values within each production-rate invocation.
No cache survives a call. The existing public interfaces are retained using
private diagnostic helpers. Three repeated timings give medians 2.064 s
without thermo reuse and 1.276 s with it; every output remains identical.
See `no_unused_derivative/` and `cached_thermo/`. Invalid-input semantics and
the wider mechanism/backend matrix have not been qualified; these changes
have not been adopted into production source.

The combined full-species/direct-temperature/order-2/thermo-reuse diagnostic
completed the saved 16-grid, 100-ns trajectory in 80.96345 s. The authoritative
comparison blocks, with unchanged case acceptance, report zero violations.
Maximum absolute momentum errors (x,y,z) are
(5.04846e-8,2.10474e-8,4.73162e-7) kg/(m2 s); maximum temperature error is
2.33174e-5 K. Density scaled error is 6.80171e-10 and total-energy scaled
error is 1.37072e-9. This is a serial diagnostic trajectory result, not MPI,
restart, convergence or physical qualification. Artifacts and exact linked
object identities are under `build/eb-cached-order2-trajectory16-20261024/`.
The checkout remains dirty main at
`0a0e00ebe9138145fbb0e0e745665d0d768df009`; the binary includes isolated
build-directory diagnostic overrides and is not a pristine HEAD build.

The same-candidate 32-grid, 100-ns comparison completed in 550.03351 s
under `build/eb-cached-order2-trajectory32-20261024/`. The unchanged
acceptance check fails seven field gates: density scaled error 3.08015e-6,
absolute x/y/z momentum errors 5.46904e-4 / 1.07307e-4 / 1.19024e-4,
total/internal energy scaled errors 4.23280e-6 / 4.23498e-6, and temperature
error 1.77513e-3 K. The y-momentum maximum is at fine cell (32,8,28);
96.5185% of its weighted absolute error lies in the fine boundary layer.
The local reaction replay improvement therefore does not establish trajectory
agreement or resolve the shared interpolation limiter sensitivity. The next
isolated stage capture targets this new maximum and adjacent (32,7,27) and
(32,8,27) cells under `build/eb-order2-max-stage32-20261024/`.
C2 and the later completion gates remain open. The existing probe object also retains its linker warning about
an executable stack; this sequence does not constitute an installed release
audit.

The targeted stage capture in `build/eb-order2-max-stage32-20261024/`
completed both baseline and combined-candidate prefixes (366.84574 s and
394.86836 s). Exit 77 is the instrumented stop after fine stage 2, not a
complete trajectory. The 108 full-cell directional-divergence checks pass.
At fine (32,8,28), stage 0 is identical and stage 1 still has identical
cell y-momentum. Its candidate-minus-baseline y-momentum rate difference
is -52.60926527, with directional contributions
(-52.60926564, 3.72529e-7, -5.03311e-10). Cell pressure differs by only
-7.19243e-7 Pa; x-face y-momentum flux differences are -0.0013710883
and +0.0068491095 on the lower and upper faces. At stage 2, cell pressure
difference is +0.27632818 Pa and the rate difference is +1245.43564,
now dominated by the y-direction contribution +1280.19810. Raw and centroid
fluxes agree exactly at the sampled full cells; tree-normalized fluxes
agree within 1.45520e-11. This locates an early x-face contribution followed
by pressure/transverse feedback before later synchronization. It does not
establish pinned-native agreement or identify the responsible donor/limiter
branch. The next capture must inspect actual ghost donors and reconstructed
Riemann inputs at this new maximum. See `origin.json`, `comparison.json`
and each variant's `directional-balance.json` for the recorded evidence.

### Coarse donor species-source decomposition (2026-10-24)

The isolated candidate prefix in `build/eb-order2-source32-20261024/`
completed in 358.727 s with the intentional `UPPERX_PREFIX_COMPLETE` exit 77.
It records 12 combined-source events and 135 spatial-input states around
coarse donor (17,4,13:15). The six distinct iteration/cell species forcings
match the captured reactor inputs exactly; each combined-source event is
recorded twice with identical values. For H2O, H, O, OH, HO2 and H2O2,
the old spatial and external terms are zero, so the forcing is exactly half
the new-time spatial term. The final-corrector spatial input state and
temperature equal the preceding predictor reaction output exactly.
The first-iteration energy rate changes between source combination and
reactor admission; species equality must not be reported as equality of
all conserved-variable rates. The stage-flux CSV matches the preceding
reaction-input diagnostic exactly. See `run.json`, `source-parts.json`,
`attribution.json`, `analyze.py` and `attribute.py` in that directory.

This narrows the next diagnostic to the advective and diffusive parts of the
new-time spatial operator. It does not establish which part causes the
pinned-native discrepancy. The separate capture under
`build/eb-order2-spatial-split32-20261024/` records both rates and six face
fluxes immediately before addition; no completed result is claimed here.
No production numerical policy or tolerance is changed by these captures.
C2 convergence and trajectory qualification remain open.

The split capture completed in 350.9524 s with the intentional exit 77.
`analyze.py` confirms nine stage/cell records, exact species-source sums and
an unchanged stage-flux CSV. `directional.py` verifies 162 species divergence
balances from time-integrated intensive face fluxes using the input time
step and spacing; all 17 target/neighbor stencil cells have unit volume
fraction. These are local operator checks, not trajectory acceptance gates.

At final-corrector coarse (17,4,14), the H2O diffusive rate is
7.3887446431e-6 versus an advective rate of 1.9315310006e-9. However,
the half endpoint difference of the forcing across k=13 and k=15 is
4.6680125307e-9 from advection and -6.8365558588e-12 from diffusion.
Their sum is 4.6611759749e-9, compared with pinned-native 4.6131828197e-9.
The absolute source magnitude therefore cannot identify the operator
responsible for the native discrepancy. A matching native split or
same-state operator comparison is still needed. `gradient.json` records
this distinction for all six trace species. No numerical fix is adopted.

The matching native split in `build/pelec-spatial-split32-20261024/`
completed in 82.8186 s. Both exported 65536-row tables are byte-identical
to the pinned reference, and the reactor trace is byte-identical to the
previous native reaction capture. `analyze.py` matches three coarse spatial
events to candidate stages 0/1/2 and verifies the native reaction context.
For final-corrector (17,4,14), candidate-minus-native H2O spatial-rate
difference is -2.3331228526e-10: diffusion contributes -2.3394171054e-10
and advection +6.2942527620e-13. H, O and OH show the same diffusion-dominated
discrepancy at this cell. Native trace-species forcing minus the averaged
old/new face divergence is at roundoff (below 1.5e-19 at this cell), so
source recombination does not explain this discrepancy.

`faces.py` resolves the diffusion difference into x/y/z contributions and
checks their sum against the independently computed operator difference
for 12 species/cell pairs. At k=14 the H2O contributions are
(-2.3416470677e-10, 1.6500794829e-13, 5.7988282964e-14); at k=15 the dominant
x contribution changes sign to +1.8648166951e-10. This locates the difference
at x-direction diffusion, but does not distinguish an input-state or
gradient difference from a transport coefficient/operator difference.
The next same-state/stencil comparison must resolve that distinction.
These captures retain `qualification=false`; no production fix is adopted.

The native state/transport capture in
`build/pelec-diffusion-state32-20261024/` completed in 84.8783 s. Both
65536-row reference tables and the preceding native reactor/flux traces
remain byte-identical. `analyze.py` compares 135 spatial-input states.
`reconstruct.py` independently reproduces all nine species fluxes at the
six sampled native x faces, with maximum absolute residual below 2.63e-17
kg/(m2 s). It uses the pinned molecular weights, harmonic face rhoD,
mole-fraction and pressure gradients, and correction velocity.

Holding native rhoD and pressure fixed while substituting candidate
composition reproduces the eight dominant H2O/H/O/OH x-diffusion-rate
differences at k=14/15 within relative residual 0.016 percent. All 12
comparisons, including the smaller k=13 differences, are recorded in
`attribution.json`. This attributes the dominant local diffusion discrepancy
to its input composition, not a demonstrated transport-coefficient error.
The neighboring x=16/18 predictor reaction inputs must now be compared
before changing integration or interpolation policy. Local flux replay
is not trajectory or physical qualification; C2 remains open.


The neighboring candidate capture in
`build/eb-reaction-neighbors32-20261024/` completed in 353.9168 s
(intentional prefix exit 77). Analysis verifies 27 paired reaction records,
unchanged prior center-cell inputs and stage-flux tables, and exact linkage
of all nine predictor outputs to the saved stage-2 spatial inputs.
`replay_substitutions.py` runs 90 isolated predictor replays (nine cells,
ten substitution modes); the candidate mode reproduces the captured state
and temperature exactly. Initial conserved states match the native inputs.

Keeping candidate inputs and changing only the mass-fraction absolute
integration tolerance from 1e-14 to native-equivalent 1e-15/rho reduces the
maximum H2O/H/O/OH density difference from 1.71234e-14 to 3.33108e-16 kg/m3
across these nine cells (about 51 times). At x=16, k=13/14/15, the H2O
errors change from [+6.48782e-15, -1.71234e-14, +1.61668e-14] to
[+3.33108e-16, +3.27871e-16, +3.32259e-16] kg/m3. Replacing all forcing
rates while retaining the candidate tolerance leaves a maximum trace error
of 2.33685e-14 kg/m3. This implicates integration error control in the
composition pattern entering diffusion; it is not evidence that a full
trajectory or interpolation limiter has been repaired. A density-scaled
native-equivalent tolerance needs an isolated full-trajectory comparison
against unchanged reference gates before any production adoption.
All results retain qualification=false; C2 remains open.


The isolated density-space absolute-tolerance trajectory candidate in
`build/eb-native-atol-trajectory32-20261024/` failed after 606.7685 s:
solver exit 1, `Reference-state MOL step rejected`. The wrapper process
itself returned zero, so the authoritative solver return code in `run.json`
is required. The unchanged comparison script correctly refused to evaluate
an incomplete trajectory. Its nine local replays exactly reproduce the
previous native-equivalent tolerance substitutions, but this does not prove
full-step admissibility. No numerical policy is adopted. Next evidence must
identify the refusal site without relaxing endpoint or rollback checks.
