# Reacting spatial-gradient AMR reference, 2026-10-19

## 2026-10-28 current-source rerun

The current GNU Release probe was rebuilt after the shared chemistry fixes.
Both imported and independently computed initial reaction histories pass the
original recursive four-time comparison through300 ns against the SHA-verified
saved pinned reference. Native simulation was reused, not rerun. Maximum
temperature error is1.209263e-5 K; maximum total-energy scaled error is
8.627650e-9 against the unchanged1e-8 limit. Relative hydrogen-consumption
error is2.243646e-5; maximum energy-budget residual is4.591537e-17.
Both comparisons execute all four times with no violations. Source and
executable hashes were independently checked before the next checker edit.
Commands, toolchain build and result logs, source hashes and verification are
in `build/tgreact-amr-current-20261028/verified.json` and `record.json`.
This is a current serial operator comparison, not MPI/application,
transport, convergence or physical qualification.

## 2026-10-28 common-clock comparison

The checker now accepts `--subcycling 0` and drives matching common level
steps in the native input and production probe. Recursive stepping remains
the default. Captured-reference manifests record the mode; a mode mismatch
is rejected before execution. Original manifests without this field retain
their recursive interpretation. Case inputs and acceptance limits are unchanged;
the copied common-clock input differs only in `amr.subcycling_mode = None`.

An unchanged pinned native TGReact run supplies the new common-clock reference.
Both imported and computed initial histories pass all four times through300 ns.
Maximum temperature error is1.209276e-5 K and maximum total-energy scaled
error is8.639131e-9, below1e-8. Maximum energy-budget residual is8.180522e-17;
hydrogen-consumption relative error is2.246787e-5. No field/inventory/activity
gates fail. The default recursive regression also passes, and deliberately
reusing the recursive reference for common stepping is rejected as expected.
Source/executable and native artifact hashes, the actual native timestep input,
all four completed times and command exits were independently verified in
`build/tgreact-amr-common-20261028/verified.json`. Exact commands are in its
`record.json`. Python syntax, project contract and diff whitespace checks pass.
This adds one serial reacting MOL2 clock combination without transport; it
is not MPI, refinement convergence, physical validation or C2 completion.

## 2026-10-28 reacting Simple-transport comparison

The regular production-integrator probe now optionally reads V3 input with
pinned LiDryer Simple viscosity, conductivity and binary-diffusion fits. The
coefficient importer checks SHA256, Git blob and species order against the
existing pinned transport manifest; its header hash must also match the
regular TGReact source contract. V2 remains the no-transport default. The
checker `--transport` sets all four native molecular-transport flags and the
same four production-integrator flags. Transport mode and fit source are
recorded in new manifests; mixed-mode reference reuse is refused.

Two new unchanged pinned-native simulations cover common and recursive clocks.
Both imported and independently computed initial histories pass in each mode:
four comparisons, each with all four physical times through300 ns. Maximum
errors across these four cases are1.234492e-5 K for temperature and8.823670e-9
for scaled total energy (unchanged limit1e-8). Maximum corrected energy-budget
residual is5.617193e-17; hydrogen-consumption relative error is2.246801e-5.
No field, positivity, species-closure, element or activity gate fails. The
rebuilt V2 default also passes, while deliberately mixing a no-transport
reference with transport enabled fails the intended mode check.

`build/tgreact-amr-transport-20261028/verified.json` independently verifies
all five completed comparisons, source/probe/native hashes, every native
artifact hash and actual diffusion flags. `record.json` retains commands and
exit codes; `core-source.json` records the core and CMake cache. GNU Release
build, Python syntax, mechanism provenance, project contract and diff checks
pass. This is serial operator reference agreement on the declared periodic
4-cubed hierarchy; it does not establish application/MPI equivalence,
mesh/time convergence, physical-boundary coverage or physical validation.

## 2026-10-28 eight-cell hierarchy: original field gates fail

The probe/checker now supports root sizes4/8/16, preserving the physical
lower-x-half refined region, original33-point profile,100-ns root step and
300-ns duration. V4 carries mesh size and transport enablement; V2/V3 remain
compatible. New manifests record root size and reject mismatched reference
meshes before execution. The rebuilt original4-cell default passes, and the
intentional4-to8 reference mismatch is refused.

The8-cell recursive, pinned-Simple-transport comparison completes all four
times for both imported and computed initial history, but each fails13 original
field checks. At100 ns, total/internal energy exceed1e-8 (total energy
1.149105e-8); at200 ns density, x/z momentum and energy fail; at300 ns all
three momenta, density and energy fail. Maximum total-energy scaled difference
is8.443683e-8 and x-momentum absolute difference4.813204e-6 kg m^-2 s^-1.
The imported/computed-history results are essentially unchanged. All inventory,
species-closure, positivity and chemical-activity gates pass; maximum corrected
energy-budget residual is7.580749e-17. Failed gates remain failed.

At100 ns, the worst total-energy cell is coarse(5,1,4), next to the refinement
boundary. At200/300 ns it is fine(1,1,12), next to the periodic coarse/fine
connection. These locations motivate transport-off and common-clock controls;
they do not prove a particular interpolation or synchronization defect.
`build/tgreact-amr-grid8-20261028/verified.json` records completed failures,
source/executable identities and original gate violations. All native artifact
hashes and four completed times were checked. No16-cell or convergence claim
is made. Controls completed under `build/tgreact-amr-grid8-controls-20261028/`.
Recursive transport-off still fails10 field gates; common-clock transport-on
fails13. Both complete all four times with inventory gates passing. At100 ns
total-energy absolute errors remain0.00270017/0.00273621 J m^-3 respectively
versus0.00273766 for recursive transport-on. At300 ns x-momentum maxima are
2.15609e-6/7.12819e-6 kg m^-2 s^-1, versus4.81320e-6. Neither removing
transport nor removing subcycling eliminates the discrepancy. All native
artifact and executable hashes were independently checked in `verified.json`.

The regular MOL wrapper passes held internal energy to the shared EB hydro
kernel without selecting its total-energy primitive option. This differs
from the already-corrected ordinary EB path. An isolated object under
`build/regular-mol-pressure-control-20261028/` selects total-energy pressure
only for nonrotating regular MOL and replays n4/n8 captured references. This
is pending diagnostic evidence, not an adopted numerical change.

## Original 2026-10-19 evidence

The production Fortran tree integrator agrees with unchanged pinned PeleC
TGReact for a nonuniform LiDryer mixture over three recursive MOL2 steps.
All declared field, reaction and inventory criteria pass after importing the
complete initial gas and reaction history. This advances C2 for a frozen
serial operator case; application initialization and the broader gate remain
open. No production numerical kernel or acceptance tolerance was changed.

## Case and result

The periodic 1 cm cube has a 4³ root and one ratio-two fine patch covering
half the domain: 320 stored cells, 288 exposed cells. Coarse steps are 100 ns,
fine steps 50 ns. The profile has a 2000 K mean temperature, 50 K cosine
variation, H2/O2/N2 mass fractions near 0.06/0.5/0.44, 101325 Pa pressure and
Taylor–Green velocity amplitude 10 m/s. Both sides use LiDryer and CVODE;
transport is disabled. Actual initial temperature span is 92.388 K and
maximum speed is 8.001 m/s. The probe imports only time zero, including
conserved gas, actual UEINT, temperature and old reaction estimates, then
calls the production AMR integrator for every later state.

| Quantity, maximum over 0–300 ns | Observed | Predeclared limit |
| --- | --- | --- |
| Initial imported state difference | exactly zero | `1e-11` scaled; exact full-state import also asserted |
| Temperature difference | `1.210e-5 K` | `1e-3 K` |
| Species mass-fraction difference | `7.231e-10` | `1e-7` |
| Total-energy difference / reference field maximum | `8.628e-9` | `1e-8` |
| Momentum difference | `1.076e-7 kg m^-2 s^-1` | `1e-6` |
| Mass and H/O/N inventory drift, either solver | `1.333e-14` relative | `1e-10` |
| Fortran total-energy residual including recorded numerical correction | `4.791e-17` relative | `1e-8` |

Reference hydrogen consumption is `4.629433349e-5` of its initial inventory;
Fortran gives `4.629329483e-5`. Their consumed-amount relative difference is
`2.244e-5` (0.002244%), below the declared 0.1% limit. Both exceed the minimum
activity requirement `1e-6`. Positive species and density/species closure
also pass. This short trajectory is not ignition or physical validation.

## Initial reaction history correction

The first comparison initialized Fortran's old reaction estimate to zero.
It completed, but total/internal energy errors reached about `1.338e-8`,
exceeding the unchanged `1e-8` gate. The pinned `Source/PeleC.cpp` post-init
path calls `react_state` for the initial level timestep; `Source/Advance.cpp`
then uses that estimate in the first MOL predictor. Gas fields alone are not
the complete initial state.

The fixture now exports the native initial `rho_omega_*` and `rhoe_dot`,
converts them from CGS to SI, and imports them with the gas state. Native
initial history is nonzero on the 32 exposed coarse cells and all 256 fine
cells; covered coarse history is zero in this pinned initialization order.
No future reaction estimate is imported. The final check asserts exact
initial identity on all 320 cells. The native gas trajectories before and
after enabling these extra output variables are numerically identical at
all four times (`history-diagnosis.json`). Thus the corrected comparison
retains the original solver and physical inputs, rather than altering the
reference trajectory or acceptance criteria.

The ordinary Fortran tree application currently initializes its fresh-start
reaction history to zero in `src/driver/amr_patch_tree_application_3d_mod.F90`.
Connecting an independently computed initial estimate, with explicit
coarse/fine ordering and preservation of restored checkpoint history, remains
application work. This operator pass must not be presented as an application
initializer pass. The previous nonreacting MOL2 register/conservation
difference also remains open.

## Reproduction and identity

Inputs and instructions: `cases/tgreact_amr_reference/`. The checker is
`tools/check_tgreact_amr_reference.py`; the explicit SUNDIALS-enabled build
target is `test_amr_tree_reacting_reference`. It uses a separate generated,
source-checked LiDryer fixture. `references/pelec_tgreact_amr.json` binds the
pinned TGReact, solver and mechanism/reactor source bytes. TGReact files were
checked against the official pinned archive before execution.

PeleF: `main`, base `0a0e00ebe9138145fbb0e0e745665d0d768df009`, with existing
uncommitted work. Reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` and
the baseline's pinned dependencies. GNU Fortran/C++ 13.3.0, CMake 3.28.3;
Fortran Debug with SUNDIALS 7.2, native C++ Release. Runs are direct serial
executions, with no MPI or performance qualification.

Artifacts under `build/tgreact-amr-reference-20261019/`:

- `runs/matched_tgreact_di4bo23x/checks.json`: final passing comparison,
  commands, binary/checker hashes, runtimes, fields and independent inventories.
- `runs/matched_tgreact_bobq9iv2/`: completed but failed zero-history attempt.
- `history-diagnosis.json`: unchanged native trajectory and comparison maxima.
- `source-identity*.json`, `source*/`, build caches, patches and build logs:
  exact source/toolchain and the pre-correction snapshot.

Executed builds:

```sh
cmake --build build/pelec-pinned-spray --target PeleC-TGReact --parallel 2
cmake --build build/reacting-hot-diagnosis-20261019/registration --target test_amr_tree_reacting_reference --parallel 2
```

Two native simulations and two Fortran simulations completed: one failed
comparison and one passed comparison. The probe was rebuilt after adding
initial reaction history. No cases were skipped inside the final comparison;
common clocks, MPI, regridding, EB, transport, refinement studies, full CTest
and installed-binary audits were not run in this increment.

Completed checks, with logs in the artifact directory:

```sh
build/cantera-venv/bin/python tools/check_mechanism_provenance.py --bundle mechanisms/lidryer.json --source mechanisms/lidryer.yaml --regenerate
python3 -m py_compile tools/check_tgreact_amr_reference.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All four passed. Acceptance limits were fixed before the first comparison
and remain identical in the failed and passed runs; only the input-output
hash changed when reaction history was added to native plotfiles.

The first isolated pressure-control build failed Fortran line-length checking;
its log remains under `regular-mol-pressure-control-20261028`. The continuation-
line-only repair runs separately under `regular-mol-pressure-control-v2-20261028`.


The v2 pressure-only diagnostic completes all three replays. The n4 original
case still passes, but n8 fails11 fields without transport and13 with transport.
At100 ns the energy discrepancy is essentially unchanged (1.133361e-8 and
1.149104e-8 scaled respectively). At300 ns x-momentum errors increase to
2.565055e-6 and5.228117e-6 kg m^-2 s^-1. Inventory gates pass. Source, linked
archives, binary and completed comparison identities are verified in its
`verified.json`. The pressure-only change is not adopted: this does not isolate
it as the cause of the observed original mismatch.

A separate explicit native-activation diagnostic is running under
`build/regular-amr-activation-control-20261028/`. It uses the pinned existing
source generator and canonical spatial algorithm; it does not combine the
pressure intervention with the chemistry intervention or replace the canonical
mechanism/reference contract.


The explicit native-activation diagnostic completes, with all23 activation
paths verified by the pinned importer. n4 passes; n8 still fails10 gates
without transport and13 with transport. At100 ns scaled total-energy errors
are1.133353e-8/1.149096e-8, essentially unchanged. At300 ns x-momentum errors
are2.156130e-6/4.813238e-6 kg m^-2 s^-1. The original spatial algorithm and
canonical source remain unchanged. Generated source, linked archives and
binary identities are verified in its `verified.json`. This control does
not remove the observed initial interface error.

A matched reaction-off diagnostic is now available in the checker/probe.
It requires zero imported reaction history and uses hydrogen conservation
rather than minimum chemical consumption; all field and other inventory
limits remain unchanged. Both `reaction_enabled` and `mol_iterations` are
recorded in reference manifests, with mismatch refusal. The original reacting
MOL2 default remains intact. The n8 no-transport reacting/nonreacting comparison
uses MOL1 to avoid the separately established pinned nonreacting-MOL2 flux
registration conflict. Its results are pending under
`build/regular-amr-reaction-controls-20261028/`.


The reacting MOL1 control completes but still fails field gates: total-energy
scaled errors are1.024599e-8 at100 ns and4.700579e-8 at300 ns. The first
nonreacting attempt completes native evolution but stops in the reader because
PeleC omits `rho_omega_*`/`rhoe_dot` diagnostics when chemistry is disabled.
This is a checker failure, not evidence of a failed Fortran evolution.
The reader now supplies zero only for known inactive reaction-history fields
and asserts any supplied history is zero. Other missing fields still fail.
A cached-native replay and default regression are pending under
`build/regular-amr-inactive-replay-20261028/`; the failed initial log is retained.

The reader repair reaches the probe but reveals a second diagnostic wiring
error: the probe supplied reaction history and callback with chemistry disabled,
which the production API correctly refuses. The probe now omits both arguments
when inactive. Kernel admission is unchanged. Rebuild/default/inactive replay
is pending under `build/regular-amr-inactive-replay-v2-20261028/`.


The corrected inactive probe completes all four times. The original reacting
n4 regression passes. The n8 nonreacting MOL1 case passes field gates at100 ns
but fails later field gates: scaled total-energy differences are3.497208e-8
at200 ns and4.259121e-8 at300 ns; x-momentum differences are2.545084e-6 and
3.109202e-6 kg m^-2 s^-1. Source/executable identity and completed command/time
coverage are verified in `regular-amr-inactive-replay-v2-20261028/verified.json`.
This establishes that chemistry is not necessary for the later discrepancy;
it does not attribute all reacting error to hydro or prove a specific coarse/
fine operator defect. The next discriminator is same-stage hydro reconstruction
and face-flux capture at the upper-x coarse/fine interface in this inactive
control (worst fine-cell x index8), distinct from the reacting periodic-side
maximum. Inactive hydrogen inventory is exactly conserved in both solvers;
maximum corrected energy-budget residual is6.259958e-17. Original failed gates remain open.


## 2026-10-28 inactive stage/face capture

Isolated candidate/native builds capture six faces at fine cells(8,2,13),
(8,1,13),(8,2,4) in the n8 nonreacting MOL1 control. The candidate uses the
existing shared MOL diagnostic; native adds observation only to its original
Diffusion source and links all other original objects. All four candidate CSVs
and budget CSVs are byte-identical to the uninstrumented baseline; all native
exported numeric rows are identical. The original comparison still fails.
No production source is modified for this capture.

The first native capture retained an unnamed time argument, causing the new
column to print a function-pointer conversion instead of physical time. It is
retained as an invalid-time diagnostic. v2 explicitly names that argument and
records12 fine-stage evaluations at0,50,...,250 ns (two evaluations per fine
step). Candidate/native event counts and dt match; their stage correspondence
is inferred from this sequence, not an independent candidate physical-time
stamp. Native phase0/1 records bracket hydro flux addition; transport is off.

At event1 the three-cell maximum x-momentum-rate difference is4.66e-8
kg m^-2 s^-2. At event2 it is0.0436967, almost entirely the x-face divergence;
the selected cell-center pressure difference remains below8.74e-11 Pa.
At event3 y-momentum-rate difference reaches0.321667, dominated by y faces.
Thus equal center pressures alone do not explain the evolving face fluxes.
The candidate records both reconstructed face pressures, but native face
reconstruction is not yet captured. Next compare event2 reconstruction and
coarse donor/ghost inputs before choosing a numerical fix.

Fluxes are instantaneous: candidate dt-integrated intensive flux is divided
by dt; native extensive CGS flux is converted to SI and divided by face area.
The candidate divergence reproduces its recorded14-component spatial rate
to below5e-15 when normalized by the sum of absolute opposing-face terms
(to account for cancellation). This is spatial hydro-rate closure, not a
full-step or reflux attribution. Data, schemas in `analyze.py`, comparison,
source/binary hashes and output-identity verification are under
`build/regular-amr-stage-capture-v2-20261028/verified.json`.


## 2026-10-28 reconstruction input localization

Native MOL observation now records left/right pressure, density, normal
velocity and sound speed plus four primitive stencil cells. The isolated
candidate records the corresponding existing diagnostic. All four native
outputs remain numerically identical and candidate fields/budgets byte-identical
to the uninstrumented inactive baseline. There are216 matched faces across
12 events and3 probe cells; duplicate native face entries, if present, are
required to agree exactly. Native CGS quantities are converted explicitly.

At event2, probe(8,2,13), the interior stencil cells differ only at roundoff,
but the exterior coarse/fine ghost differs by-2.133417e-4 Pa and-3.197201e-6
m/s in x velocity (candidate minus native). Its density difference is
-2.174959e-10 kg m^-3. This input difference reaches the two reconstructions
of the last valid cell with opposite pressure changes of about5.33354e-5 Pa
and normal-velocity changes of7.99300e-7 m/s.

A432-record same-input replay sends each captured candidate and native stencil
through the unchanged production MOL face routine. Candidate inputs reproduce
all eight recorded face quantities exactly. Native inputs reproduce pressures
within1.455192e-11 Pa, densities within2.775558e-17 kg m^-3, normal velocities
within3.726880e-14 m/s, and sound speeds exactly. This localizes the selected
mismatch upstream of face reconstruction, in the supplied ghost state. It
supports no claim yet about temporal interpolation versus spatial interpolation
or about unobserved cells. The next audit captures coarse start/end donors and
interpolation fractions for this event. Production numerics remain unchanged.

The replay initially produced a trailing CSV delimiter; its output formatter
was corrected and the complete replay rerun. Commands, source/binary hashes,
stencils, units/conversion code and output identity are retained under
`build/regular-amr-reconstruction-capture-20261028/verified.json`. Python syntax,
project contract and diff whitespace checks pass.


A coarse-donor observer now records candidate start/end states and time fractions
for seven-point stencils around coarse(5,1,7) and(5,1,2), with periodic wrapping.
The native observer records temporally prepared coarse interpolation inputs
and the actual FillPatcher time/stage. Two preparation attempts failed because
insertion markers matched multiple routines; their directories contain only
stage-capture reruns and are not donor evidence. The successfully prepared
version is pending under `build/regular-amr-donor-capture-v3-20261028/`.
Production sources remain unchanged; output-identity checks precede analysis.


## 2026-10-28 coarse endpoint localization

The v3 donor capture completes with all four outputs unchanged. All168
seven-point donor rows across12 stages match coordinates after periodic
wrapping and explicit CGS/SI/component mapping. At event2 the candidate
fraction0.5 corresponds to native FillPatcher time50 ns; event3 uses the same
time as its old state. Candidate interpolation starts from recorded coarse
start/end states exactly once (`set_ancestors` copies the start state; the
reader applies the blend). There is no observed double time interpolation.

Differences already exist in coarse states. At event1 the donor center(5,1,7)
is equal to roundoff, but its covered neighbor(4,1,7) differs in total energy
by-0.0884135273 J m^-3 and temperature by-0.00665788757 K; held energy differs
only by-2.91e-11 J m^-3. At event2 the center's blended density and x-momentum
differences are-2.174959e-10 kg m^-3 and-3.307274e-7 kg m^-2 s^-1, consistent
with the previously observed ghost primitive discrepancy. The coarse endpoint
update therefore contributes before spatial interpolation.

Source inspection shows candidate wrapper-entry and old-stage energy resets,
whereas native MOL fills its old-state border before its first hydro evaluation.
A diagnostic object bypasses those two resets only in the matched nonrotating
MOL fixtures and compares n4 reacting/n8 inactive/n8 transport cases. It is
pending under `build/regular-amr-old-state-control-20261028/`. This broad bypass
is not a production change or a valid general admission policy; any adoption
must preserve validation, rollback and held/rotating contracts. Donor evidence
and hashes are in `regular-amr-donor-capture-v3-20261028/verified.json`.


The old-state-only control completes all three cases: n4 passes, n8 inactive
still fails9 original gates with essentially unchanged errors, and reacting
n8 transport fails11 (previously13). The latter's100-ns scaled energy error
falls to7.302260e-9 and300-ns x-momentum error to3.945380e-6 kg m^-2 s^-1;
conservation still passes. This is partial diagnostic improvement, not an
adopted reset policy. Source/archive/binary evidence is retained in its
`verified.json`.

The covered initial donor difference corresponds to the difference between
restricted total kinetic energy and kinetic energy computed from restricted
momentum; a held-based reset removes that difference. Ordinary candidate
pressure still uses held internal energy, so preserving total energy alone
need not change hydro. A combined isolated control now pairs old-state
preservation with the previously tested total-energy primitive option under
`build/regular-amr-old-total-factorial-20261028/`. It uses the same three
references and unchanged gates. Both individual controls remain separately
recorded; no production bypass is adopted.

The combined control has now completed all four output times in each of the
three cases. Independent verification confirms unchanged production split
source, diagnostic source, pressure object, linked archives and probe hashes.
All inventory gates pass, but all three field comparisons fail: n4 has two
violations at 300 ns (scaled total energy 1.652521e-8 against 1e-8); n8 inactive
has nine (300-ns x-momentum absolute error 3.588465e-6 kg m^-2 s^-1); reacting
n8 transport has six (300-ns x-momentum 1.202249e-6 against 1e-6 and scaled
total energy 2.466402e-8). The latter improves from thirteen baseline failures
and eleven old-state-only failures, but the newly failing n4 case rules out
adopting this combined bypass as a fix. The wrapper exits zero to finish the
experiment; each comparison checker exits one on its original field gates.
Evidence is in `build/regular-amr-old-total-factorial-20261028/verified.json`.

The next diagnostic should compare coarse old/end states across the first
coarse-step synchronization and the second-step entry, where errors reappear
after the much closer first output. Preserve separate total and held energies
and the restricted kinetic-energy difference. Admission, rollback, rotating
contracts and acceptance thresholds remain unchanged; C2 is still open.

## 2026-10-28 synchronized-state entry capture

An isolated observer around coarse synchronization and wrapper entry preserves
all 24 field/budget files byte-for-byte against the combined control. Source,
object, archive and executable identities pass independent verification.
The original field failures remain 2/9/6; this is observational evidence only.
At the second step's entry, repeated restriction changes covered-cell total
energy by up to 0.4836103262 J m^-3 (n4), 0.2654831543 (n8 inactive), and
0.2656706295 (n8 reacting/transport). Held energy and temperature are unchanged;
all three quantities are unchanged on coarse leaf cells. These energy changes
reverse the preceding covered-cell synchronization reset. The third step shows
the same pattern. Capture, output identity and analysis are retained under
`build/regular-amr-sync-capture-20261028/` (`analysis.json`,
`entry_comparison.json`, `verified.json`).

Pinned native `post_timestep` refluxes, averages down, then computes temperature;
`do_mol_advance` swaps time levels and fills the old border without repeating
restriction. A diagnostic omitting wrapper-entry restriction in addition to the
previous old-state preservation/total-pressure control is running under
`build/regular-amr-entry-preserve-20261028/`. It keeps the same three cases and
original gates. No production admission or rollback bypass is adopted.

The entry-preservation control passes all four times and original inventory
checks in all three cases. At 300 ns, scaled total-energy errors are
3.075171e-12 (n4), 1.067618e-10 (n8 inactive), and 2.649424e-10 (n8 reacting
with transport); x-momentum absolute errors are respectively 1.600537e-10,
3.103386e-9 and 1.843448e-8 kg m^-2 s^-1. Source/object/archive identity was
verified before editing production files; see its `verified.json`.

A production candidate now preserves synchronized old data in the serial
nonrotating stored-energy MOL path with total-energy updating enabled. It
validates both total and held thermodynamic states without modifying accepted
old data; the ordinary regular MOL face adapter uses total-energy pressure.
Rotating and independent-held update policies retain the previous entry path.
A focused regression exercises invalid density, total energy, held energy and
temperature with unchanged state and diagnostic outputs on rejection. Debug
split/rotating/Riemann tests and Release native comparisons are pending under
`build/regular-amr-admission-20261028/`. Mechanism provenance, project contract,
Python syntax and diff whitespace checks pass. MPI's corresponding entry
restriction/reset remains to be aligned and tested; this is not a completed
serial/MPI correction or C2 qualification.

The production serial candidate passes the rebuilt Debug split (including four
invalid-input rollback checks), rotating-reaction and MOL Riemann executables.
The rebuilt Release n4/n8-inactive/n8-transport comparisons pass all original
field and inventory gates at four times each. All 24 field/budget outputs are
byte-identical to the successful entry-preservation diagnostic. Current source
hashes and results are verified in
`build/regular-amr-admission-20261028/verified.json`.

The MPI counterpart now preserves the same ordinary old state and performs
owner-local thermodynamic validation followed by collective logical agreement.
A regression injects invalid density or held energy on one owner and checks
refusal/rollback on every rank, including an empty rank with three processes.
Six rebuilt 2/3-rank stored-energy/MOL/rotating tests and four application
lifecycle modes are pending under `build/regular-amr-mpi-admission-20261028/`.
No broader MPI, mesh/time or physical qualification is claimed yet.

The first MPI verification completes with four application lifecycle modes
passing (20 runs) and both rotating MPI tests passing. Four exact comparisons
fail: stored-energy 2/3-rank cases disagree at the second reacting MOL iteration,
and ordinary MOL diagnostics disagree despite exact state and inventory output.
Failure artifacts remain in `regular-amr-mpi-admission-20261028`.
Inspection identifies serial auxiliary-source reevaluation at iteration 2,
where the pinned native and MPI paths rebuild only spatial rates. Serial now
uses the same ordering. MPI reflux also recomputed non-species correction
magnitudes by subtracting large updated/input energies; its diagnostic now
retains the directly accumulated correction and accounts for repaired species
separately. These changes do not loosen comparison gates. Serial/MPI/restart
rechecks are pending in `regular-amr-mpi-parity-fix-20261028`.
A separate current Release rebuild reruns the three cached references and adds
n16 reacting with transport under `regular-amr-grid16-current-20261028`.

The current-source n4 and both n8 reference comparisons pass again. The first
n16 attempt completes the native reference capture but the candidate subprocess
hits the checker's fixed 120-second wall-time limit before field comparison;
this is neither a numerical gate failure nor a completed trajectory. The timed
out child is confirmed absent. The checker now accepts a finite positive
`--command-timeout` (default unchanged), records it per command and marks
subprocess timeouts. Zero, negative, NaN and infinite values are rejected before
filesystem work. The captured n16 reference is reused with a 900-second limit
under `build/regular-amr-grid16-resume-20261028/`; numerical gates are unchanged.

The parity recheck is now independently verified: one serial split test, six
MPI tests (2/3 ranks, retained-energy, ordinary MOL, rotating MOL), and all four
lifecycle modes pass. The latter comprise 20 application runs; saved serial,
MPI2, MPI2-to-MPI3 restart and serial-restart field/budget snapshots are exact
within each mode. Source hashes remain unchanged. Evidence and snapshot hashes
are in `build/regular-amr-mpi-parity-fix-20261028/verified.json`. The initial
failed attempt remains separate evidence, not counted among these passes.

The n16 rerun completes in 324.79 seconds including exports and comparison,
passing all four field times and original inventory gates. Maximum scaled
energy error is 5.181104e-9 (200 ns); maximum x-momentum absolute error is
3.722402e-7 kg m^-2 s^-1 (200 ns), at fine(1,13,14). The maximum energy-budget
relative residual is 6.835686e-17. Current source hashes are verified in
`build/regular-amr-grid16-resume-20261028/verified.json`. These reference errors
are larger than n8 and do not demonstrate mesh convergence.

A matched n8 temporal study now uses 100/50/25-ns coarse steps through the same
300-ns endpoint. The checker records a time-refinement factor, rejects cached
references with a different factor, and retains the original scientific gates.
The unchanged default and a deliberate mismatch refusal are checked alongside
the new native/candidate runs in `build/regular-amr-time-refinement-20261028/`.
Its prepared analysis compares both solvers at common 100/200/300-ns times,
using composite leaf volumes for L1 and maximum absolute differences for Linf;
these self-differences are distinct from fixed-reference agreement.

## 2026-10-28 bounded regular-AMR convergence evidence

The n8 transport temporal study completes: dt50 ns passes all seven native
field times; dt25 ns passes all thirteen, with original inventory gates.
The unchanged default passes and a wrong temporal-factor cache is refused.
At common 100/200/300-ns times, all 16 gas/species/temperature observables
have decreasing volume-weighted L1 and Linf self-differences in both solvers.
Observed candidate orders span 0.9268–1.4604 (L1) and 0.8719–1.5280 (Linf);
native spans 0.9268–1.4589 and 0.8719–1.5280. These are approximately first-order
bounded temporal observations, not a claim of second-order accuracy. Exact
formulas, all 96 solver/time/field rows, input hashes and gate verification are
in `build/regular-amr-time-refinement-20261028/{temporal_convergence,verified}.json`.

The current n4 transport rerun also passes all four native field times.
Conservative factor-two, same-level volume restriction compares n4/n8/n16
composite leaf states at dt100 ns. Both solvers show decreasing state L1/Linf
for every field at 100/200/300 ns. To expose initialization dependence, a
second comparison restricts each mesh's U(t)-U(0). All increment L1 values
still decrease, but x-momentum increment Linf increases at 100 and 200 ns
in both solvers (candidate orders -0.30190 and -0.10307; native -0.30190 and
-0.10307). State differences alone would hide this limitation. Temperature
is volume-averaged as a diagnostic, not recovered from restricted energy.
This fixed-dt three-grid result does not prove coupled asymptotic space/time
convergence. All 224 rows, definitions and 24 source CSV hashes are in
`build/regular-amr-spatial-convergence-20261028/{spatial_convergence,verified}.json`.
C2 spatial maximum-norm qualification and the separate EB open gates remain
unclosed; no physical validation is claimed.

The regular-AMR spatial increment audit localizes the n8/n16 x-momentum
self-difference: all top 144 points at 100/200/300 ns lie in the immediate
one-cell layer at a coarse/fine interface in both solvers. At 100 ns the
candidate Linf is 0.01636147 kg m^-2 s^-1 there, versus 0.00721463 in the
remaining 1664 comparison cells. Unweighted spatial error-pattern cosines
across momentum components, both mesh pairs and all three times exceed
0.9999999999787. This supports shared interface sensitivity but is not a
causal operator identification or an exemption from convergence requirements.
One/two/three local-cell bands and definitions are preserved in
`build/regular-amr-spatial-convergence-20261028/localization.json`.

The matched regular probe now has a shared-source MPI build target and checker
launcher options. Owner-local stepping and reaction-history initialization call
the production MPI routines; output gathering alone replicates fields. Serial,
MPI2 and MPI3 (one empty rank), each with imported/computed initial histories,
are being compared directly to the same original n8 transport reference under
`build/regular-amr-mpi-reference-20261028/`. This is pending direct reference
evidence, separate from the completed application/restart matrix.

## 2026-10-29 MPI reference build isolation

The first direct MPI attempt builds but does not enter the numerical probe:
all four MPI launches fail to load `libmpi.so.12`. Inspection finds that the
old Release CMake cache combines OpenMPI Fortran interfaces with an Intel MPI
core library and Intel hwloc. Serial results in that attempt pass, but no MPI
scientific result is claimed. The failed logs and cache are preserved.
A fresh GNU/OpenMPI Release directory is configured with a minimal environment
and explicit system MPI/hwloc libraries. A clean-environment `ldd` audit must
pass before the six reference comparisons. This rerun is pending under
`build/regular-amr-mpi-clean-reference-20261029/`; it reuses the same frozen
native reference and unchanged numerical gates.

The fresh OpenMPI build completes all six n8 recursive-transport comparisons:
serial/MPI2/MPI3, each with imported and independently computed initial reaction
history, pass all four original field times and inventory gates. For each
history mode, all eight CSV/budget files are byte-identical across ranks.
The clean-environment link audit resolves OpenMPI only. Source, cache,
executable and output hashes are independently checked in
`build/regular-amr-mpi-clean-reference-20261029/verified.json`. This closes the
bounded direct MPI comparison for this case, not the whole CPU/MPI matrix,
restart qualification, spatial convergence or physical validation.

An extension uses the same verified clean binaries for common-clock n8
serial/MPI3 comparisons with both history modes, plus MPI3 at n16 recursive
and n8/dt25 ns. Original native references and gates are retained. These six
runs are pending under `build/regular-amr-mpi-reference-extension-20261029/`;
the verifier also checks matching serial CSV/budget artifacts, including all
13 dt25-ns outputs. This extends bounded reference coverage and does not
replace the separate spatial maximum-norm or EB obligations.

The six MPI extension cases complete and pass their original field/inventory
gates: common-clock n8 serial/MPI3 with both history modes, n16 recursive
MPI3, and n8/dt25-ns MPI3 through 300 ns. Matching serial field and budget files
are byte-identical, including all thirteen dt25-ns times. Source, binary and
cache identities were verified before the subsequent boundary-probe edit;
see `build/regular-amr-mpi-reference-extension-20261029/verified.json`.

The next independent C2 boundary comparison disables x periodicity and sets
both x faces to native `SlipWall` / candidate `slip_wall`, leaving y/z periodic.
The checker/probe record the new mode and reject inconsistent cached references.
Fresh n4 transport-off/on serial and MPI3 comparisons, plus the unchanged
periodic default and intentional cache refusal, are pending under
`build/regular-amr-slip-reference-20261029/`. No existing periodic pass is
counted as physical-boundary agreement.

The n4 slip-wall matrix completes: transport off/on, serial/MPI3, all four
times pass the original native field and inventory gates. Native copied inputs
are verified to have x periodicity disabled and both x boundaries set to
SlipWall. All eight candidate field/budget files match exactly between serial
and MPI for each transport mode. The periodic default still passes and a
wrong-boundary reference cache is refused. Source and artifact verification
is in `build/regular-amr-slip-reference-20261029/verified.json`.

A follow-up records cumulative candidate boundary flux and external-source
integrals for mass, three momenta and total energy. The checker compares these
with composite inventory changes (including numerical energy correction),
normalizing signed-momentum cancellation by absolute inventory and using the
existing inventory tolerance. This adds a wall-impulse ledger check rather
than inferring impulse closure from field agreement. Periodic and n4 slip-wall
controls are rerunning, and a new n8 transport slip-wall serial/MPI3 pair is
pending under `build/regular-amr-boundary-ledger-20261029/`. No production
flux or update formula is changed by the observer/checker addition.

The boundary-ledger extension passes all seven completed comparisons: periodic
n4 default, slip-wall n4 transport off/on serial/MPI3, and slip-wall n8 transport
serial/MPI3. All original field gates and the new mass/momentum/energy ledger
checks pass at four times each. Maximum normalized gas-budget residual is
3.431661e-15; all matching serial/MPI field and budget files are byte-identical.
Sources and artifacts are verified in
`build/regular-amr-boundary-ledger-20261029/verified.json`.

The next boundary control selects stationary adiabatic x `NoSlipWall` in
native and `no_slip_wall` in the candidate, preserving periodic y/z, initial
profile and original field/ledger tolerances. The probe rejects conflicting
wall modes; reference manifests distinguish slip and no-slip. N4 transport
off/on serial/MPI3 comparisons are pending in
`build/regular-amr-noslip-reference-20261029/`. This tests tangential reflection
and molecular stresses separately from the completed slip-wall family.

The n4 no-slip controls pass all four times and original field plus gas-ledger
gates for transport off/on in serial/MPI3. Boundary input mapping and exact
serial/MPI field/budget identity are independently verified in
`build/regular-amr-noslip-reference-20261029/verified.json`. At 300 ns the
inviscid no-slip/slip fields are identical; with molecular transport, maximum
no-slip/slip differences are 2.022956e-4 kg m^-2 s^-1 in y momentum,
0.3936591 J m^-3 in total energy and 7.445917e-4 K in temperature. Thus the
viscous wall distinction is exercised, not merely a renamed periodic case.
These cross-boundary differences are diagnostics, not reference errors.
An n8 transport extension with common/recursive clocks in serial/MPI3 is now
pending under `build/regular-amr-noslip-grid8-20261029/`.

The n8 no-slip extension passes all four original field times and gas-ledger
gates for common/recursive clocks in serial/MPI3. Corresponding outputs match
byte-for-byte. Maximum normalized gas-budget residual is 6.057225e-15;
wall mass and total-energy flux integrals are zero, while the nonzero normal
wall impulse closes the momentum inventory. Verification is retained in
`build/regular-amr-noslip-grid8-20261029/verified.json`.

Because regular AMR MOL entry/source ordering is shared with experimental
spray/LES callbacks, a clean Debug/OpenMPI spray-enabled regression build is
now exercising the existing AMR particle-reactive, source/force, coupled LES
and serial/MPI common/recursive LES tests. Original test names and timeouts
are retained; the selected CTest inventory and JUnit results are recorded in
`build/regular-amr-coupling-regression-20261029/`. This is an impact check of
shared code, not a physical spray result or completion of C3.

The clean Debug/OpenMPI coupling-impact build completes all 11 selected
original tests, with zero failures, disabled tests or skips (523 s CTest time).
They cover AMR particle reaction/split/fixed/force paths, coupled LES at 2/4
ranks, and serial/MPI common/recursive LES-SDC. Independent JUnit verification
matches all planned names, verifies all recorded source hashes and records
cache/binary hashes. Clean-environment loader audits exclude missing libraries
and Intel/OpenMPI mixing. See
`build/regular-amr-coupling-regression-20261029/verified.json` and `tests.xml`.
This completes the focused shared-code regression check; it does not establish
native particle trajectories, physical LES/spray validation, or C3/C5 completion.

### 2026-10-29: first-order extrapolation boundary comparison

The regular-AMR reference checker and shared serial/MPI probe now expose
`--x-boundary outflow`, mapped to pinned native `FOExtrap` and PeleF `outflow`
on both x faces, with periodic y/z. This is zero-gradient extrapolation, not
an imposed-pressure or nonreflecting boundary. The native parser mapping was
checked in `Source/PeleC.cpp`, and its scalar/normal/tangential components all
map to `foextrap` in `Source/Setup.cpp`.

Open-domain probe budgets additionally export every species' boundary/source
integral. Candidate mass, momentum, energy and elemental inventory changes
are checked against supplied amounts, with the existing limits. Native mass
and elemental inventories are compared with candidate inventories; native
boundary flux is not exported, so no independent native conservation claim
is made. H2 inventory consumption is adjusted for candidate H2 boundary supply
for both solvers, explicitly labeled in the report. Original closed-domain
checks and all field thresholds are unchanged.

The clean GNU/OpenMPI Release build and four n4 recursive comparisons pass:
serial and MPI3, each with molecular transport disabled/enabled, at all four
output times through 300 ns. Each serial/MPI pair has byte-identical field and
budget CSVs. The maximum candidate gas-budget relative residual is
9.778e-16; the elemental residual maximum is 3.629e-16. Net mass supply at
300 ns is nonzero: -2.9270e-15 kg without transport and -3.1877e-15 kg with
transport. The default periodic comparison also passes after rebuilding;
reuse of a periodic native capture as outflow is correctly refused.

Artifacts and reproducible verification:
`build/regular-amr-outflow-foextrap-20261029/verify.py`, `verified.json`,
`record.json` and per-case `checks.json`. These record source and binary hashes,
commands and result identities. Dirty main remains HEAD
`0a0e00ebe9138145fbb0e0e745665d0d768df009`. The first attempt under
`build/regular-amr-outflow-reference-20261029/` used the internal name
`Outflow` as a native input keyword; both native runs rejected it before
calculation. Those failures remain recorded and are not counted as comparisons.
The wrapper's completed status does not override their nonzero command codes.

Mechanism provenance, Python compilation, project contract and diff whitespace
checks pass. No production numerical kernel was changed. This closes only the
n4 recursive extrapolation comparison; finer grids, common clocks, other open
boundary prescriptions, moving/isothermal walls and C2 convergence remain open.

### 2026-10-29: n8 extrapolation with both AMR clocks

The n8 molecular-transport comparison now passes common and recursive clocks,
each in serial and MPI3, through all four original field-comparison times to
300 ns. The four checker runs returned zero (46.05, 33.54, 62.71 and 49.92 s,
including native generation/export where applicable). Both serial/MPI pairs
have byte-identical field and budget CSVs: eight files per pair. All candidate
open gas/element ledgers pass; maxima of their relative residuals are
3.2474e-15 and 2.9439e-16 respectively. Net mass supplies are -3.91223e-15 kg
(common) and -3.90867e-15 kg (recursive). Native inventory comparison remains
distinct from an independent native boundary-ledger audit.

`build/regular-amr-outflow-grid8-20261029/verify.py` independently checks all
command return codes, four successful reports, exact case modes, four output
times, eight field/budget files per pair, nonzero mass transfer and source/
binary identities. `verified.json` stores its output and artifact hashes.
The source was unchanged during these runs. No production kernel, input
physics, acceptance tolerance or reference solver was changed. These are
bounded trajectory and MPI checks, not spatial/temporal convergence or C2
completion. Temporal refinement with the open boundary is the next distinct
check; finer-grid agreement alone does not establish convergence.

### 2026-10-29: extrapolation-boundary temporal refinement

For the n8 recursive molecular-transport case, dt50ns and dt25ns serial runs
pass every original field and open-inventory gate through300ns (7 and13 output
times). MPI3 dt25ns independently passes all13 times and matches all26 serial
field/budget CSVs byte for byte. Maximum gas and elemental budget residuals
are4.346e-15 and3.211e-16. Checker elapsed times were100.80,166.53 and141.86s.
Artifacts: `build/regular-amr-outflow-time-20261029/`; `verify.py` checks source
identity, commands, modes, output counts and exact MPI equality, then runs the
saved-field `analyze.py`. `verified.json` and `temporal_convergence.json` retain
results and input hashes.

At100/200/300ns, all16 observables have decreasing L1 and Linf self-differences
for both solvers over100/50/25ns steps. L1 is volume-weighted absolute difference;
Linf is the maximum absolute difference. Observed order is
log2(D100_50/D50_25), separate from fixed-reference agreement. Candidate order
ranges are0.96750..1.79149 (L1) and0.88688..1.96789 (Linf); native ranges are
0.96750..1.78980 and0.88688..1.97423. These results do not demonstrate uniform
second order, spatial convergence, other open-boundary prescriptions or C2
completion. The n4/n8/n16 spatial comparison is the next distinct check.

### 2026-10-29: n16 extrapolation and nonmonotone spatial differences

The n16 recursive transport serial and MPI3 comparisons pass all four original
field/open-inventory times through300ns. Checker times were454.75 and388.90s.
Their eight field/budget CSVs match byte for byte. Maximum gas/element relative
budget residuals are4.245e-15/1.646e-16. Artifacts are under
`build/regular-amr-outflow-grid16-20261029/`; `verify.py` checks modes, counts,
source/binary identities, command results and MPI identity before analysis.

The n4/n8/n16 same-level conservative-restriction analysis does **not** establish
spatial convergence. Both solvers show the same nondecreasing set: six state
Linf entries and23 evolution-increment Linf entries, plus z-momentum L1 at300ns
for both state and increment. Increment is U(t)-U(0) on each mesh before
restriction. Thus the initial profile alone does not explain the issue.
Candidate z-momentum L1 grows from0.0013063941 to0.0013621379 kg m^-2 s^-1
(order -0.06028), while its Linf decreases slightly. No failed norm is waived.
These fixed100ns comparisons are not a coupled space-time asymptotic study.

`localize.py` checks384 regional field/time/scope/grid decompositions against
the complete saved global norms. Regions are the one-cell physical-x layer,
the one-cell coarse/fine-x layer and the remaining cells on each coarser
comparison mesh; regional L1 contributions use the full-domain volume. This
is localization, not evidence that a particular boundary operator is wrong.
For n8/n16 increments, density and total-energy maxima occur in the physical-x
layer at100/200/300ns. The x-momentum maximum moves from that layer to the
coarse/fine layer at300ns. The z L1 growth is not confined to either boundary:
physical-layer contributions decline0.000597717 to0.000545903, coarse/fine-layer
contributions decline0.000452682 to0.000327755, and the remaining contribution
increases0.000255996 to0.000488480 kg m^-2 s^-1. These mesh-dependent layers do
not define a fixed physical subvolume, so this is not a support-growth proof.

Across the compared signed spatial error fields, the volume-weighted cosine
between candidate and native is at least0.999999999934. This confirms a shared
pattern, not physical correctness or permission to change reference gates.
`spatial_convergence.json` and `localization.json` retain definitions, all
locations, regional norms and input SHA256 identities. Python compilation,
mechanism provenance, project contract and diff whitespace checks pass.
No production numerical change is justified by localization alone. The next
cause-discriminating comparison must separate boundary/profile compatibility,
transport and time-step contributions before a correction; further grid
refinement alone is not the next action. C2 remains open.

### 2026-10-29: initial compatibility of the extrapolation fixture

Before changing a boundary operator, `audit_initial_boundary.py` in
`build/regular-amr-outflow-transport-control-20261029/` checks the actual initial
boundary-adjacent states against pinned TGReact initialization. Across80/320/1280
cells on n4/n8/n16, normal velocity matches the analytic TG expression within
4.0e-15m/s and temperature matches the piecewise-linear profile exactly. Native
`prob.H`, `prob.cpp`, `prob_parm.H`, boundary-mapping source and all consumed
CSV/profile hashes are retained in `initial-boundary.json`.

The normal TG velocity has nonzero normal derivative at x=0,L (sampled maxima
5363,6044,6223s^-1), whereas `FOExtrap` imposes equal ghost/interior values.
The fixed33-point profile also has endpoint temperature slopes of approximately
-3074.36/+3074.36K/m and H2 mass-fraction slopes of -0.122974/+0.122974m^-1.
Consequently this periodic-profile adaptation is not a smooth zero-gradient-
compatible boundary convergence fixture. It remains a valid bounded translation
comparison, and its measured nonmonotone norms remain recorded. This audit
alone does not prove that the compatibility issue causes the later error.
A matched transport-off n8/n16 control is pending to separate molecular
transport sensitivity before modifying any boundary formula or adding a
separate compatible convergence fixture. No baseline or gate is replaced.

### 2026-10-29: transport-off control retains nonmonotone spatial differences

Both n8/n16 recursive transport-off comparisons pass all original field and
open-inventory gates through300ns. Initial candidate CSVs are byte-identical to
their transport-on counterparts. The n4 transport-off member comes from the
previous completed extrapolation comparison. `verify.py` under
`build/regular-amr-outflow-transport-control-20261029/` verifies modes, initial
identity, command status and all reports, then recomputes the224 spatial rows
and records paired transport-on/off orders in `verified.json`.

Removing molecular transport does not remove the spatial nonmonotonicity.
At300ns candidate z-momentum increment L1 order is-0.08652 without transport
versus-0.06028 with it; the n8/n16 L1 magnitude ratio is0.98642. Total-energy
increment Linf order is-0.01183 versus-0.01210, with magnitude ratio1.00803.
At100/200ns, transport-off total-energy Linf is also nonmonotone. Thus molecular
transport is not necessary for the observed issue; this is not evidence that
transport has no effect or that the boundary operator is defective.

A separate `--initial-velocity zero` checker mode now changes native `prob.u0`
to zero and imports the corresponding state, requiring exactly zero initial
momentum. It records/checks the mode in cached-reference manifests. The ordinary
TG minimum-speed guard remains intact. No production kernel or comparison
threshold changed. This control removes all initial TG flow and cannot by
itself attribute differences solely to the boundary derivative. The n4/n8/n16
transport-off zero-flow campaign is pending under
`build/regular-amr-outflow-zero-flow-20261029/`. Its profile remains the original
33-point temperature/composition table, including the endpoint slopes audited
above. It does not replace the original failed convergence evidence.

### 2026-10-29: zero-flow control and an unresolved momentum ledger failure

The zero-initial-flow, transport-off campaign completed its numerical outputs.
n4 passes every original check. n8/n16 pass native field comparisons but fail
the unchanged candidate x-momentum budget limit1e-10: n8 step1 is1.154e-10;
n16 steps1..3 are4.726e-10,4.198e-10,3.052e-10. These runs are **failed**
trajectory checks, not successful comparisons. Absolute x-momentum residuals
are8.845e-22kg m/s at n8 step1 and3.907e-21..2.126e-20kg m/s at n16.
No denominator or tolerance was changed. Separate sums of opposing face fluxes
in the serial/MPI regular-tree ledger are a cancellation candidate, not yet a
proven cause; actual face values must be captured before a production change.

`build/regular-amr-outflow-zero-flow-20261029/audit_controls.py` verifies the
completed outputs, original failures, zero initial momentum and matching
density/temperature/species. The initial total-energy difference agrees with
removed kinetic energy within8ULPs (maximum1.162e-10J/m3). Covered coarse cells
require the average of child kinetic energy, not kinetic energy of averaged
momentum. The first audit incorrectly used the latter and failed; correcting
the accounting did not alter any solver output. Initial reaction histories are
exact on n4/n8; n16 has an energy-reaction-source difference4.657e-4J/m3/s,
so complete history identity is not claimed. Species histories remain exact.

The separate saved-field diagnostic explicitly permits only the recorded
x-momentum budget violations and labels its result unqualified; it does not
change `checks.json` or reclassify failed trajectories. With zero initial flow,
z momentum remains exactly zero, and the n8/n16 total-energy increment Linf
at100/200/300ns is7.192e-5/1.307e-4/1.941e-4 times the TG transport-off value.
Its spatial Linf orders become0.671/0.967/1.202. This shows dependence on the
initial TG flow, not specifically on its boundary derivative: the control
removes all TG advection and strain as well. The original spatial failures and
the new ledger failures remain open.

`control-audit.json` contains initial-state checks, all original budget failures,
paired spatial diagnostics and hashes. `analyze.py` retains the fixed100ns
space/time limitation. The checker default periodic regression passes, and
captured TG references are refused for zero-flow reuse. Python compilation,
LiDryer provenance, project contract and diff whitespace checks pass. No
production numerical kernel was modified in this diagnostic.

### 2026-10-29: captured root-face summation explains n8 ledger discrepancy

An isolated copy of the serial split module records actual opposing root
boundary faces and the separately applied reflux correction. Its n8 zero-flow
run reproduces all eight original field/budget CSVs byte for byte, including
the failed budget. Capture, compile/link commands and binary/source identities
are in `build/regular-amr-ledger-faces-20261029/record.json`.

`analyze.py` recomputes each root boundary integral with80-digit Decimal
arithmetic from the captured binary64 faces (192 pairs per step). Only root
ledger arithmetic is replaced; the solution and reflux corrections remain
frozen. The cumulative residual changes from8.845e-22/1.109e-21/1.204e-21 to
1.714e-23/6.850e-23/-1.037e-23kg m/s at steps1..3. A paired-face `math.fsum`
replay agrees with each high-precision root result to displayed binary64
precision. This directly identifies the loss from subtracting separately
accumulated opposing boundary fluxes in this captured case; it does not prove
all larger-grid or MPI effects are resolved.

A scoped serial/MPI implementation now subtracts each opposing face pair
before summing, leaving fluxes, state updates, reflux, normalization and gates
unchanged. Pre-edit copies are in `build/regular-amr-ledger-paired-20261029/`.
Current-source n8/n16 serial/MPI3 reruns and focused Debug regressions are
pending. Adoption requires original field/budget gates, exact field identity
against the failed baselines, and exact serial/MPI field and budget agreement.
The existing failure records remain authoritative until this verification.

### 2026-10-29: paired boundary-face ledger correction verified and adopted

The serial and MPI root boundary budgets now sum opposing face differences
rather than subtracting two independently accumulated face sums. The actual
face fluxes, gas state updates, reflux corrections, energy policy, acceptance
normalization and tolerances are unchanged. The task-only patch and pre-edit
copies are in `build/regular-amr-ledger-paired-20261029/`.

All four n8/n16 zero-flow serial/MPI3 comparisons pass the original native
field and candidate open-inventory gates at all four times through300ns.
Every solution CSV is byte-identical to the corresponding pre-change failed
baseline. Each serial/MPI pair also has byte-identical budget CSVs. The maximum
x-momentum relative residual is3.098e-11, below the unchanged1e-10 criterion;
the earlier maximum was4.726e-10. The failed baseline records remain preserved.
This closes the scoped root-ledger cancellation defect without changing the
computed trajectories or reclassifying their spatial-convergence evidence.

The clean GNU/OpenMPI Release build completed in35.26s. Four checker times
were42.53/87.91/329.63/677.45s for n8 serial/MPI3 then n16 serial/MPI3. The
existing clean GNU/OpenMPI Debug build completed in14.70s; all five focused
tests pass in336.45s, with zero skips: serial split and SDC, MPI split ranks2/3,
and MPI SDC ranks4. Exact commands, JUnit, source/binary identities and output
hashes are in `record.json`, `tests.xml` and independently generated
`verified.json`; reproduce with `python3 build/regular-amr-ledger-paired-20261029/verify.py`.

LiDryer provenance, project contract, Python compilation and diff whitespace
checks pass. This is a bounded numerical-ledger correction, not a full release
matrix, general spatial-convergence result or physical validation. C2 remains
open, including the original TG-flow spatial maxima and the separate EB and
nonreacting-MOL2 obligations.

###2026-10-29: common slope sign-test regression

Following the separate extreme-value limiter correction, fresh n8 zero-flow
outflow serial/MPI3 runs pass all original four-time field and inventory gates.
Their state and budget CSVs are byte-identical to the paired-ledger baseline.
See [sign-test evidence](slope-sign-underflow-20261029.md) for the captured
failure, new scale-invariance regressions and source identities. This does not
change the existing TG-flow spatial-convergence result or close C2.

###2026-10-29: explicit isothermal-wall reference input

The checker and serial/MPI probe now accept an explicit positive finite
x-wall temperature. It is recorded in reports/native reference manifests;
reusing a reference at a different temperature is rejected. Native inputs
set `do_isothermal_walls` and both x-face temperatures. The initial pilot
incorrectly left candidate wall types adiabatic while supplying temperatures:
all four serial/MPI wall cases fail18 field gates each, with zero candidate
wall-energy ledger. These failures remain in
`build/regular-amr-isothermal-20261029/` and are fixture-routing failures,
not evidence about the correctly selected isothermal operator.

The candidate now selects existing `isothermal_wall`, which includes no-slip
velocity reflection. The explicit temperature option therefore requires
`no_slip_wall` and transport; isothermal slip remains unsupported by this
fixture and an open coverage item. NaN/zero/negative temperatures and missing
wall/transport combinations reject before work-directory creation. Default
adiabatic options retain an absent optional temperature array. No production
flux or tolerance changes accompany this probe extension.

The n4, 500 K, no-slip recursive serial/MPI3 comparison passes the original
four-time field and candidate boundary-energy ledger gates through 300 ns.
All eight state/budget CSVs are byte-identical between serial and MPI3. The
cumulative outward energy transfer is 1.2580890167552766e-5 J at 300 ns
(candidate signed boundary energy: -1.2580890167552766e-5 J). Reusing the
500 K reference with a 600 K request rejects before work-directory creation.
Commands, identities and checks are retained under
`build/regular-amr-isothermal-routed-20261029/`, including `verified.json`.
This checks candidate ledger closure and native field agreement; an independent
native boundary-energy ledger has not been exported. No spatial-convergence
or physical-validation claim follows from this bounded comparison.


The n8 extension also passes all four original field/inventory times through
300 ns for common and recursive clocks, each in serial and MPI3. Each clock
has eight byte-identical serial/MPI state/budget CSVs. A periodic-default
control passes and preserves all eight CSVs from the earlier boundary-ledger
baseline. The two modified probe/checker source hashes remained unchanged
through the run. Provenance, GNU/OpenMPI link audit, Python syntax, project
contract and whitespace checks pass. The initial project-contract invocation
omitted its required `--source` argument; the corrected invocation passed.
Reproduce the output/hash verification with
`python3 build/regular-amr-isothermal-grid8-20261029/verify.py`; exact commands,
source/binary identities and timings are in the adjacent `record.json`.
This extends bounded stationary no-slip thermal-wall coverage only. Wall
mesh/time convergence, isothermal slip, moving walls and the other C2 gates
remain open.


The n8 recursive 100/50/25 ns study is complete under
`build/regular-amr-isothermal-time-20261029/`. The 50 ns serial and 25 ns
serial/MPI3 cases pass all 7/13/13 original field and inventory times through
300 ns. All 26 finest-step state/budget CSVs are byte-identical serial/MPI3.
The maximum gas-budget relative residual is 2.12e-16; closed H/O/N inventory
drift across candidate and native is at most 2.89e-15.

At 100/200/300 ns, both solvers have decreasing composite-volume-weighted L1
and Linf temporal self-differences for every one of 16 observables (96
solver/time/field rows). Observed log2 ratios range from 0.9385 to 1.7059 for
L1 and 0.8490 to 1.9424 for Linf. This is decreasing temporal error, not
uniform second-order convergence or a spatial-convergence claim. Exact
commands, source/binary hashes, norms and reproduction scripts are in
`record.json`, `temporal_convergence.json`, `verify.py` and `verified.json`.
The verifier initially expected the open-boundary element-ledger key; closed
wall cases instead use the checker's closed element inventories. The verifier
was corrected to measure those drifts; no solver or acceptance gate changed.


### 2026-10-29: thermal-wall spatial convergence remains open

The n16, 500 K recursive serial/MPI3 comparison passes all four original
reference times through 300 ns, with all eight state/budget CSVs byte-identical.
The n4/n8/n16 spatial audit nevertheless fails to show convergence for several
observables in both solvers. At 300 ns, candidate evolution-increment Linf
increases from 0.03796 to 0.05823 kg m^-2 s^-1 for x momentum, from 160.41 to
315.88 J m^-3 for total energy, and from 0.58415 to 1.73894 K for temperature.
Energy and temperature L1 also increase. Thus reference agreement is not
spatial qualification. Conservative restriction, separate initial/increment
norms, commands and hashes are retained in
`build/regular-amr-isothermal-spatial-20261029/` (`verify.py`,
`spatial_convergence.json`, `localization.json`).

For both solvers at all three evolved times, the n8/n16 energy and temperature
maxima lie in the first wall-adjacent cell layer, seven comparison-cell layers
from the coarse-fine interface. At 300 ns all top 144 energy/temperature and
x-momentum differences lie in that wall layer; transverse momentum maxima
remain at the coarse-fine interface. The n16 initial wall-adjacent gas spans
2049.04--2049.52 K against the 500 K wall. This motivates a wall-temperature-only
2050 K control at n4/n8/n16, now running under
`build/regular-amr-wall-jump-control-20261029/`. It tests sensitivity to the
initial thermal mismatch; it does not replace the 500 K baseline, prove an
asymptotic limit, or justify changing production fluxes or acceptance gates.


### 2026-10-29: wall-temperature intervention separates two error regions

The 2050 K wall-only control completes n4/n8/n16 serial trajectories and passes
all four original reference times on every mesh; source hashes remain fixed.
At 300 ns, n8/n16 evolution-increment temperature Linf decreases from 1.73894 K
(500 K wall) to 0.21495 K, and total-energy Linf from 315.884 to 118.374 J m^-3.
Their maxima move from the wall-adjacent layer to the coarse-fine interface.
The corresponding y/z momentum maxima change by only 0.0142%/0.0030%; both
remain at the interface. Thus changing wall temperature strongly affects the
thermal wall extrema but does not remove the transverse-momentum interface
error. This is a controlled sensitivity result, not a decomposition of all
coupled causes or justification for altering conservative fluxes.

The 2050 K energy/temperature L1 and Linf now decrease across n4/n8/n16, but
only weakly: at 300 ns the temperature orders are 0.0811/0.2067 and energy
orders 0.1236/0.1936. Both solvers still have nonmonotone x-momentum increment
Linf at 100/200 ns and slightly increasing z-momentum L1 at 300 ns. Eight
solver/scope/time/norm entries remain nondecreasing. Neither wall temperature
therefore supplies general spatial qualification. The original 500 K failed
spatial result remains unchanged. No production numerical fix is adopted.

`build/regular-amr-wall-jump-control-20261029/verified.json` records the three
reference checks, source identities and temperature intervention ratios;
`localization.json` separates physical wall and coarse-fine cell layers.
Reproduce with `analyze.py`, `localize.py`, then `verify.py` in that directory.
Python syntax, project contract and diff whitespace checks pass. This control
is serial; the earlier 500 K MPI equivalence does not establish MPI equivalence
for every new wall temperature. The remaining interface investigation must
compare stage pressure/flux and synchronization increments, retaining the
separate wall sensitivity and original reference gates.


An isolated current-source MOL observer is running in
`build/regular-amr-interface-stage-20261029/` for the 2050 K n8/n16 cases.
It observes the n8 transverse-momentum comparison cells (coarse 5,3,1 and
fine 8,15,12) and their eight n16 children. Per-call records contain cell
pressure, left/right reconstructed face pressures, all directional hydro
fluxes and the returned hydro rate. The observer must preserve all 16 baseline
state/budget CSVs. `analyze.py` checks directional divergence closure and
compares volume-averaged child rates. These are hydro-only instantaneous
observations; call order is not yet a physical stage/time mapping, and
transport, reaction, reflux and synchronization remain separate contributions.
Results are pending. The production source is unchanged.


The first observer completes with all 16 baseline CSVs unchanged. However,
reconstructing the divergence from dt-integrated flux divided back by dt fails
the diagnostic 1e-12 criterion: normalized residuals reach 3.69e-11 (n8) and
1.33e-9 (n16). These failed checks are retained in `analysis.json`; they are not
promoted to a successful closure result. The existing MOL diagnostics retain
centroid flux before multiplication by dt. A second isolated observer in
`build/regular-amr-interface-unscaled-20261029/` records that value directly,
with the same closure criterion and output-identity requirement. Its results
are pending. Source call order maps each nine-event coarse interval to three
root evaluations followed by three evaluations for each of two fine intervals:
old state, predictor state, then state after the first corrector. Only old and
last evaluations enter the accepted hydro trapezoid; this mapping does not
include transport or reflux contributions.


The second observer also preserves all 16 outputs, but its closure residual
is order one. Source inspection identifies an observer-only indexing error:
`centroid_flux` already occupies the physical box and `crop_valid_flux` uses
`g=0` in this call, whereas reconstructed pressure occupies the larger margin
box. Applying the pressure array's +1 offset to centroid flux sampled the
neighboring face. These failed records remain preserved; no physical
attribution uses them. A third observer in
`build/regular-amr-interface-physical-flux-20261029/` uses direct physical-box
face indices and retains the original 1e-12 diagnostic criterion. It is running.
`attribute.py` requires that closure before integrating accepted hydro stages
and comparing them with the actual conserved-momentum evolution. Its residual
is explicitly the unresolved transport/synchronization/arithmetic contribution,
not automatically a reflux error. Production source and reference gates remain
unchanged.


The physical-box observer completes with all 16 baseline CSVs byte-identical.
The original 1e-12 divergence criterion now passes: maximum scaled residuals
are 2.65e-16 (n8) and 4.41e-16 (n16), across 27/216 observed cell-events.
The preserved earlier failures distinguish dt multiply/divide cancellation
from the subsequent observer indexing mistake. Source identities and output
checks are in `build/regular-amr-interface-physical-flux-20261029/verified.json`.

Integrating only accepted old/last MOL evaluations gives, at 300 ns, total
coarse-grid-minus-restricted-fine y-momentum evolution difference -0.01647056
kg m^-2 s^-1 at coarse cell (5,3,1), with hydro contribution -0.01638242 and
unresolved remainder -8.81424e-5. The y-direction hydro contribution is
-0.01621480. At fine cell (8,15,12), z-momentum total difference is -0.00612363,
hydro -0.00589073, remainder -2.32903e-4; z-direction hydro contributes
-0.00589317. These signed components describe specific targets, not a global
causal fraction. `attribution.json` retains all three times and directions.
The numerical residual must still be split into transport and synchronization;
it is not valid to identify it as reflux alone.

A separate isolated synchronization observer is now running under
`build/regular-amr-interface-sync-20261029/`. It captures those cells before
reflux, after reflux and after thermodynamic recovery, again requiring all
16 baseline CSVs to remain unchanged. The existing pressure/flux observations
and this direct synchronization increment are complementary. No production
fix or relaxed criterion has been adopted.


The synchronization observer completes with all 16 baseline CSVs unchanged
and source hashes verified. At 300 ns, observed coarse-cell y-momentum reflux
contributes +8.06576e-5 kg m^-2 s^-1 to the n8-minus-restricted-n16 difference,
opposite to the total -0.01647056. The fine-cell z-momentum reflux contribution
is exactly zero. Thermodynamic recovery changes neither target's momentum.
After hydro and synchronization, unresolved y/z differences are respectively
-1.687999e-4 and -2.329030e-4 kg m^-2 s^-1. Thus direct reflux cannot explain
the dominant transverse difference at these targets; it remains coupled to
the earlier pressure/flux histories rather than being an independent cause.
The initial analyzer misread interval fractions 0..1 as seconds; it now checks
those fractions and the recorded phase order, using the known coarse-step
sequence. No numerical output changed.

Evidence is in `build/regular-amr-interface-sync-20261029/analysis.json` and
`verified.json`. An accepted-rate observer is running in
`build/regular-amr-interface-transport-20261029/` to capture old/new molecular
transport rates, accepted hydro rates and actual momentum updates. Its closure
checks must pass before assigning the unresolved part to transport. The
production source and original comparison criteria remain unchanged.


The accepted-rate observer preserves all 16 CSVs and closes the target-cell
momentum accounting. Per-update absolute residuals are at most 5.47e-17
kg m^-2 s^-1; rate-decomposition scaled residuals are at most 1.01e-16.
Across both targets, three times and three components, hydro + molecular
transport + reflux + recovery reproduces the observed grid-evolution
difference to 8.43e-17 kg m^-2 s^-1. At 300 ns, transport supplies -1.687999e-4
to coarse-target y momentum and -2.329030e-4 to fine-target z momentum,
matching the previously unresolved contributions. Exact commands, identities,
checks and signed components are in
`build/regular-amr-interface-transport-20261029/analysis.json` and `verified.json`.
This is arithmetic contribution accounting, not independence or causality of
operators that evolve coupled states.

A further isolated observer in `build/regular-amr-interface-pressure-20261029/`
now captures actual Riemann interface pressure/velocity to split diagonal
momentum flux into pressure and mass-flux-times-velocity. It intentionally
repurposes two diagnostic-only kinematics slots in the isolated object; the
production diagnostic interface and solver are unchanged. All 16 baseline
CSV identities and the flux decomposition are required before interpreting
these records. Results are pending. Reconstructed left/right pressures remain
in the same trace; Riemann interface pressure must not be equated with their
arithmetic average or attributed solely to interpolation.


The Riemann observer preserves all 16 baseline CSVs. On every observed face,
normal momentum flux equals mass flux times interface velocity plus interface
pressure exactly in the saved binary64 values. The integrated pressure/advection
split closes to below 1e-12 kg m^-2 s^-1. At 300 ns, coarse-target y-direction
hydro is -0.01621480, comprising pressure -0.01591030 and advection -0.00030450.
Fine-target z-direction hydro is -0.00589317, comprising pressure -0.00589090
and advection -2.27726e-6. These are signed differences, not independent causal
percentages.

An algebraic pressure split further distinguishes the targets. Coarse-target
y pressure has reconstructed left/right mean contribution +0.00513018 and
Riemann departure -0.02104048; fine-target z pressure has mean -0.00603344 and
departure +0.00014254. Units are kg m^-2 s^-1 after accepted-stage integration.
The departure includes velocity/impedance/wave selection and cannot be called
interpolation error. Records and hashes are in
`build/regular-amr-interface-pressure-20261029/{analysis,pressure_split,verified}.json`.
No production change follows from this decomposition alone.

A same-input replay capture is running under
`build/regular-amr-interface-stencil-20261029/`, saving the four primitive-cell
states and sound speeds feeding each observed face. The physical-box face
index maps to the prepared primitive box by +3, with stencil offsets -1..2.
Replay must first reproduce actual flux and left/right pressure before any
order or component intervention is interpreted. Output identity remains a
required gate, and the prior spatial-convergence failures remain open.


The stencil observer completes with all 16 baseline CSVs unchanged. Replaying
all 1,458 observed faces through the current clean-build library reproduces
both full flux vectors and reconstructed left/right pressures exactly. The
same-input first-order control is not an improvement: at 300 ns, the coarse
y-direction integrated difference changes from -0.0162148 to -0.0298386
kg m^-2 s^-1, while the fine z-direction changes from -0.00589317 to -0.00503179.
This rules out adopting first-order reconstruction based on this diagnostic;
it is not an evolved first-order convergence result.

Additional frozen-stencil component controls keep sound speeds and all other
primitives fixed. Flattening normal velocity to the central-pair mean changes
the coarse y-direction difference to +0.00182661, whereas flattening pressure
to its central-pair mean changes it to -0.02191190. For fine z the corresponding
values are -0.00605549 and -0.00511319. These controls reveal different algebraic
sensitivities, consistent with the earlier pressure decomposition. They are
not thermodynamically consistent physical interventions. In particular,
flattening pressure within each four-cell stencil retains differences between
neighboring face means, so it does not remove all pressure-gradient effects.
No interpolation, Riemann, limiter or synchronization modification is justified
by these controls alone. General C2 spatial convergence remains unproved.

Exact commands, source/binary hashes and data are under
`build/regular-amr-interface-stencil-20261029/`: `replay.py` verifies baseline
reproduction and first-order controls, `component_replay.py` performs component
controls, and their analysis scripts retain the accepted-stage sums.
`verified.json` records identities and the complete 1,458-face exact baseline.
The evidence now supports tracing the coarse-target normal-velocity stencil
history and the fine-target pressure-mean history separately, rather than
attributing both maxima to a single coarse-fine synchronization defect.


### 2026-10-29: initial-versus-evolved stencil history

A further replay restores only pressure or normal velocity from each face's
first old-stage stencil, retaining current sound speeds and other primitives.
The untouched branch again exactly reproduces all 1,458 flux/pressure records.
At 300 ns, restoring initial pressure changes the fine-target z-direction
integrated difference from -0.00589317 to +0.00049034 kg m^-2 s^-1; restoring
initial normal velocity gives -0.00605474. Coarse-target y instead retains
-0.01751747 when initial normal velocity is restored (baseline -0.01621480),
while initial pressure gives -0.02358032. The coarse normal-velocity sensitivity
therefore exists in the initial stencil too; the fine z pressure effect is
strongly tied to evolved pressure. These frozen-input replacements are not
consistent fluid states or evolved trajectories and do not establish causality.

The captured fine z-face stencil pressure is initially uniform on each mesh
(event 4, span zero). At the first predictor evaluation (event 5), the union
pressure span is 24.427 Pa for the n8 target and 14.701 Pa across its eight n16
children. This span is neither a derivative nor a pointwise intermesh error.
It localizes the first observed nonuniform pressure to the predictor update,
prior to the fine-step chemistry corrector, rather than an initial pressure
reconstruction artifact. The predictor includes old reaction/transport/hydro
sources, so those contributions still require separate examination.

`history_replay.py`, `analyze_history.py`, `history_control.json`,
`pressure_history_ranges.json` and `history_verified.json` under
`build/regular-amr-interface-stencil-20261029/` retain these controls and input
hashes. Python syntax, project contract and whitespace checks pass. The next
bounded discriminator is the first fine predictor's conserved-state source
contributions to pressure; coarse initial velocity/Riemann sensitivity remains
a separate question. No production numerical change is adopted, and C2 remains
open.


A first-fine-predictor observer is now running under
`build/regular-amr-first-predictor-20261029/`. It records all valid cells used
by the target z-face stencils (five n8 cells and twenty-four n16 cells), before
and after predictor energy recovery. Records include old/predicted conserved
states, old/new temperature and held energy, hydro, combined predictor rate,
molecular transport, external source and reaction history. The same-build
EOS replay executable is prepared with the probe's pinned EOS policy.
Before source substitutions, the replay must reproduce event-5 stencil
pressures and the observer must preserve all 16 original CSVs. Pressure is
nonlinear in conserved variables; additive source accounting does not imply
additive pressure contributions. Results are pending; no production source
has changed.


### 2026-10-29: first predictor pressure is generated primarily by hydro

The first-predictor observer preserves all 16 baseline CSVs. Both the saved
predictor update and its hydro/transport/external/reaction decomposition close
exactly in the evaluated arithmetic. Same-library EOS replay reproduces all
29 event-5 cell pressures exactly. All frozen-source substitution states are
accepted without clipping. Evidence and identities are in
`build/regular-amr-first-predictor-20261029/{analysis,verified}.json`.

For the five n8 stencil cells, hydro alone produces pressure range
101302.9461--101334.8435 Pa (span 31.89745 Pa). Reaction alone is spatially
uniform at 101324.5307 Pa; transport alone spans only 0.0002594 Pa. Combining
sources before recovery gives span 31.89720 Pa; after energy recovery it is
24.42724 Pa. Across the twenty-four n16 cells, hydro-only span is 19.19876 Pa,
reaction-only 0.005864 Pa and transport-only 0.050175 Pa; after recovery the
span is 14.70099 Pa. These union spans are not pointwise errors or additive
pressure fractions. They support hydro as the dominant local generator of
this first pressure pattern, while recovery modifies its amplitude.

A frozen-hydro EOS scaling check evaluates Uold+a*dt*H for a=1,1/2,1/4,1/100.
After dividing each pressure change by a, reducing a from1 to1/100 changes
individual n8/n16 pressure responses by at most 0.000333/0.000151 Pa, compared
with full-response magnitudes 22.054/14.309 Pa. The pressure response is thus
already present in the frozen spatial hydro rate; finite-step EOS nonlinearity
is small in this control. This is not a coupled temporal-convergence result.
`hydro_scale.py` and `hydro_scale.json` preserve the definition and all cells.

The next source-level discriminator is the pressure response to the hydro
material/energy/momentum rates and their directional fluxes, retaining the
native reference, conservative update and shared-limiter behavior. A source
replacement, tolerance relaxation or lower-order scheme is not justified by
these observations. Python syntax, project contract and whitespace checks
pass; C2 spatial convergence and the broader completion gates remain open.


### 2026-10-29: normal ghost velocity controls the initial fine pressure mismatch

Frozen-hydro EOS accounting separates material (density plus all species),
momentum and total-energy increments using all six group orderings. Material
and energy pressure contributions mostly reinforce here; momentum contributes
at most 0.00605 Pa on n8 and 0.00155 Pa on n16. The grouping is nonlinear
algebraic accounting, not causal independence (`group_pressure.json`).

At the target fine cell and its eight n16 children, six-face hydro divergence
closes to the captured rate. Splitting by direction before the same EOS gives
n8-minus-mean-n16 pressure difference -7.95662 Pa: x contributes -8.10802 Pa,
y +0.15140 Pa and z zero (`directional_pressure.json`). The later z acceleration
therefore need not originate in a z-directed initial flux: x/y updates first
generate nonuniform pressure. This is the first predictor only.

For the first fine old-stage x faces, a replay substitutes exact initial TG
normal velocity only in ghost cells outside the fine patch. Interior values
match the analytic field within 1e-12 m/s; all other primitives and sound
speeds remain captured. The resulting hydro-only predictor pressure difference
changes from -7.956622 Pa to +0.000263393 Pa. This controlled sensitivity
identifies the observed ghost normal velocities as important for this local
mismatch. It is not a proposed interpolation: analytic ghosts are unavailable
in general, and evolved/native/conservation gates are not replaced.

Inputs, exact baseline replay and controls are in
`build/regular-amr-interface-stencil-20261029/analytic_ghost.py` and
`analytic-result.json`; EOS results and hashes are in
`build/regular-amr-first-predictor-20261029/analytic_ghost_pressure.json` and
`directional_verified.json`. A donor observer is now running in
`build/regular-amr-initial-donor-20261029/` to record each component's centered
and limited slope and the shared coefficient for the relevant coarse cells.
The pinned shared-limiter policy must be distinguished from any implementation
defect before changing production interpolation. Results are pending.


The first donor observer preserves the 16 solution/budget CSVs but produces
no donor records: its filter used level0, whereas this reader uses one-based
source levels (root=1). That run is not successful donor evidence. A corrected
observer under `build/regular-amr-initial-donor-level1-20261029/` limits execution
to one coarse step by changing only the requested step count, and requires all
eight step0/step1 CSVs to match the existing trajectories.

The completed n8 first-step record identifies density as the first x-direction
shared-limiter restriction: center 0.1035354713 kg m^-3 exceeds both neighbors
(0.1034611260 and 0.1012834210), so its limited slope and the common coefficient
are zero. The x-momentum component's own slope ratio is1, but the common zero
coefficient suppresses it too. These are ordinary-sized gradients, not a
subnormal product ambiguity. The pinned `mf_cell_cons_lin_interp_llslope`
source computes the same per-direction minimum across components; this is a
source-formula comparison, not an independently rerun native donor replay.
The n16 first-step observer is still running. `analyze.py` will verify the
saved coefficient arithmetic and exact outputs before finalizing the result.


### 2026-10-29: shared donor limitation reproduced in pinned AMReX

The corrected root-level observer completes for both meshes with all eight
step0/step1 CSVs unchanged. The saved directional coefficient arithmetic is
reproduced exactly. Both n8 and n16 donors are density extrema, which sets the
shared x coefficient to zero while the x-momentum component's own ratio is1.
For n16, density center 0.1037767180 exceeds neighbors 0.1037367598 and
0.1031455340 kg m^-3. The forward/backward gradient products are -6.697e-7
(n8) and -1.009e-7 (n16): ordinary extrema, not product underflow. The n8 shared
y/z coefficients are also zero; n16 y is1 and z is zero.

`native_replay.py` compiles a small caller of the actual pinned AMReX
`mf_cell_cons_lin_interp_llslope` kernel, without copying its arithmetic. Using
the identical captured SI donor values, ratio2 and interior boundary codes,
all 96 slopes (two donors, three directions, sixteen components) match exactly.
The kernel header, AMReX archive, probe binary, commands and inputs are hashed
in `native-verified.json`; combined checks are in `verified.json` under
`build/regular-amr-initial-donor-level1-20261029/`. This is same-input native
kernel evidence, not a newly captured full native CGS donor trajectory.

This establishes the observed initial normal-velocity flattening as shared
limiter behavior for these donors, not a candidate arithmetic mismatch.
Replacing shared scales with componentwise or analytic velocity interpolation
would change the fixed-reference numerical policy and therefore is not adopted.
The analytic-ghost pressure improvement remains diagnostic, not a completed
fix or convergence result. Broader maximum-norm convergence remains open even
though fixed-reference trajectories and these isolated kernel checks pass.
The two failed observer configurations remain preserved and are not counted
as successful captures. Python syntax, project contract and diff whitespace
checks pass.


### 2026-10-29: bounded finer interface-resolution comparison

The observed initial sound speeds are 1129--1139 m/s. Over300ns this is only
0.54 fine-cell widths at n8 and1.09 at n16 (fine dx0.625/0.3125mm). This scale
estimate motivates one n32,2050K recursive comparison to examine a more
resolved interface disturbance; it is not proof of asymptotic convergence or
a reason to disregard the existing failed maxima.

The checker and probe now allow `--cells32` in addition to4/8/16, without
changing default4, acceptance limits, physical domain, refinement topology or
the pinned33-point interpolated profile. The clean GNU/OpenMPI probe rebuild,
provenance and link audit pass. A default n4 periodic control passes original
gates and reproduces all eight earlier CSVs byte for byte; cells64 is rejected
before work-directory creation. Python syntax, project contract and whitespace
checks pass. The n32 serial trajectory is running under
`build/regular-amr-grid32-interface-20261029/`; exact commands and source/binary
identities are in `record.json`, and default/invalid controls are in
`extension_controls.json`. `analyze.py` compares n8/n16/n32 state and evolution
increments. The fixed33-point profile remains piecewise interpolated, so this
is not a smooth analytic-profile convergence test. MPI32-grid equivalence,
original500K spatial failures and other C2 gates remain separate obligations.


### 2026-10-29: n32 completes but fails unchanged reference gates

The n32,2050K serial comparison completes both native and candidate runs
(native942.67s, candidate2715.94s; checker3810.65s including exports). It is
not a timeout. Four times are present, native AMReX finalization is recorded,
and source hashes are unchanged. Nevertheless, the original field checker
fails18 criteria: density, all momentum components and total/internal energy
at each evolved time. At100ns x momentum differs by7.98455e-6 kg m^-2 s^-1
against1e-6, and total energy scaled error is6.47612e-8 against1e-8, both at
fine cell(32,6,61). No inventory gate fails. These failures are preserved in
`build/regular-amr-grid32-interface-20261029/verified_failure.json`; they are
not superseded by output completion or convergence trends. MPI qualification
for n32 has not been attempted.

The n8/n16/n32 spatial analysis also remains incomplete: at300ns the candidate
x/y/z momentum increment Linf differences decrease from0.035537/0.016471/
0.006124 to0.025509/0.006833/0.005891 kg m^-2 s^-1, but z-momentum Linf increases
at100/200ns in both solvers. Eight solver/scope/time rows remain nondecreasing.
The fixed native shows closely corresponding spatial behavior, without
excusing the candidate's failed direct reference gates. State and increment
norms, exact definitions and input hashes remain in `spatial_convergence.json`.

The next run, `build/regular-amr-grid32-donor-20261029/`, records the coarse
donor at zero-based(16,2,30), associated with the first-time maximum. It runs
only the first coarse step and must reproduce the four corresponding candidate
CSV files from the failed comparison. This diagnoses initial donor arithmetic;
it does not itself reproduce or explain the full native/candidate discrepancy.
No acceptance limits, profile or production numerical policy changed.

### 2026-10-29: n32 donor arithmetic and native capture boundary

The isolated candidate donor observer completed the first coarse step and
reproduced all four step0/step1 state and budget CSV files exactly. Its120
interpolation events contain5760 directional component slopes. Replaying
those saved SI inputs through the actual pinned AMReX limiter reproduces
all5760 slopes exactly. The initial shared x coefficient is zero because of
a y-momentum extremum; the x-momentum component alone would admit coefficient1.

Converting those same saved inputs to CGS before the native-kernel replay
changes the shared coefficient by at most4.05329828e-6. The worst case is
event21, direction3:0.197041000164 in SI versus0.197045053462 after conversion.
O2 density differences of order1e-12 determine that coefficient. This measures
rounding sensitivity of nearly constant donor fields, not the actual native
donor trajectory, and does not explain away the18 failed reference gates.
Commands, source hashes and exact-output checks are retained under
`build/regular-amr-grid32-donor-20261029/` in `record.json`, `analysis.json`
and `unit-replay.json`.

A separate observation-only native build under
`build/regular-amr-native-grid32-donor-20261029/` records the same zero-based
coarse point(16,2,30) and its six neighbors in native CGS component order.
It uses the current pinned TGReact link inputs with only the observed AMReX
interpolater object substituted. Its first-step exports must equal the saved
unobserved native exports before the donor trace can support comparison.
The run record owns its completion status; initiation is not a passed check.

The native observation run completed normally (AMReX finalization; native
runtime369.60s). Both initial and100-ns exported tables match the unobserved
native tables byte-for-byte. The trace contains six interpolation calls,
each with seven points and16 components. `compare_donors.py` explicitly maps
native `[rho,momentum3,E,e,T,species9]` to candidate
`[rho,momentum3,E,species9,T,e]` and converts CGS to SI. The first native
stencil matches candidate events1--20 exactly in all112 values. Native calls
2--4 are value-nearest to candidate events21--80; calls5--6 to81--120.
These are value matches, not independently observed stage identities.

Replaying actual native CGS inputs through the pinned kernel gives a
z-direction shared coefficient0.174476749972 for call2, versus candidate
SI0.197041000164 for event21. Both coefficients are controlled by O2 density.
The native O2 center/minus/plus, converted to SI, are approximately
0.05190851594548345/0.05190851594526368/0.05190851595030201 kg m^-3;
the candidate values are0.05190851594698098/0.05190851594673054/
0.05190851595181453. Tiny differences in an almost constant species field
therefore change a common coefficient acting on all momentum components.
The resulting y-momentum z-slope differs by-6.21992305e-4 kg m^-2 s^-1
(native minus candidate). This is a measured interpolation-output difference,
not a demonstrated final-cell error budget or a justified limiter change.

`donor_comparison.json`, `actual_replay.json`, and `limiting_components.json`
retain the comparison, pinned-kernel results, and controlling components.
The native input replay preserves the original CGS bytes, avoiding a CGS-to-SI
round trip. Candidate conversion by division versus multiplication itself
changes last bits, so its converted-input coefficient is only a sensitivity
control. The next discriminator is explicit interpolation time/stage context
and coarse donor update accounting, followed by propagated face-flux and
momentum closure. Production numerics and all acceptance gates remain unchanged.

The saved middle donor stencil is exactly the arithmetic midpoint of the
saved initial and final donor stencils in each solver (all112 values, native
calls1/2/5 and candidate events1/21/81). This algebraic check is retained in
`temporal_inputs.json`; it does not by itself identify physical times.
An explicit time/stage observer and root MOL border/RHS/final-state observer
are now isolated in `build/regular-amr-native-donor-context-20261029/`.
Its output identity and completion remain subject to its run record.

The explicit native context observer completed normally in367.48s and its
two exported field tables match the unobserved native outputs. The six donor
calls occur at fine-level times0,50,50,50,100,100ns with50-ns fine steps;
stages are old,predictor,corrected,old,predictor,corrected. Every donor value
is reproduced exactly from the observed root old/final states using linear
time interpolation (six calls x112 values). See
`build/regular-amr-native-donor-context-20261029/analysis.json`.

The native root O2 field starts with zero z-neighbor differences. Its
predictor differences become-2.30677e-12 and-1.17321e-12 g cm^-3;
the corrected differences are9.82815e-15 and1.01649e-13, and final differences
4.39549e-16 and9.63711e-15. Large cancellation across MOL updates leaves the
small final gradient that controls the shared interpolation limiter. This
locates the native history but does not yet isolate the candidate/native
update mismatch. An isolated candidate root-stage observer is running under
`build/regular-amr-candidate-root-stages-20261029/`; its comparison requires
four unchanged first-step candidate CSVs before interpretation.

The candidate root-stage observer completed successfully and reproduced all
four first-step candidate CSVs. At the seven selected root points, the initial
state matches native in every component. The initial hydro-plus-transport RHS
already differs: maximum O2 rate difference8.20769792e-6 kg m^-3 s^-1.
The predictor O2 state differs by8.20767065e-13 kg m^-3. Thus the mismatch
precedes coarse-to-fine interpolation. Initial O2 z-rate differences are
-0.0230667763/-0.0117315971 in the candidate versus
-0.0230676668/-0.0117320791 kg m^-3 s^-1 in native. The final small gradients
then differ after cancellation across successive updates. This identifies an
earlier discriminator, not the unique cause of all downstream field failures.
Full component and stage records are retained in
`build/regular-amr-candidate-root-stages-20261029/analysis.json`.

A native face-flux observer in
`build/regular-amr-native-root-flux-20261029/` now records extensive diffusion
fluxes before adding hydro, and total fluxes afterward, on all six faces of
these seven cells at each root MOL stage. Its planned divergence comparison
must first pass unchanged native output checks and close against the saved
native RHS. No production algorithm or acceptance tolerance was changed.

The split native face-flux observer completed normally in364.88s with unchanged
initial/first-step native exports. At the initial root stage, the maximum O2
candidate-minus-native rate differences over the seven points are
8.20775524e-6 kg m^-3 s^-1 for transport and1.72873771e-10 for hydro.
Total extensive-face divergence closes against the saved O2 RHS within
1.94511074e-13 kg m^-3 s^-1. Thus the initial O2 RHS mismatch is localized to
transport in this stencil. Later stage hydro differences grow after the states
have diverged; they do not establish an independent initial hydro defect.
Component arrays and closure residuals are in
`build/regular-amr-native-root-flux-20261029/analysis.json`.

The next isolated native observer records left/right primitive states and
harmonic face transport coefficients alongside the diffusion fluxes under
`build/regular-amr-native-transport-inputs-20261029/`. This will distinguish
primitive-state preparation from coefficient/flux arithmetic using identical
inputs. Its run record must establish output identity before replay results
are interpreted. Production transport and reference gates remain unchanged.

The primitive/coefficient observer completed normally with unchanged native
exports. Same-native-primitive replay through the current Fortran
`face_transport_data` and `species_face_flux` gives maximum coefficient
relative difference1.05333e-15 and O2 flux difference4.14790e-17 kg m^-2 s^-1
over126 faces/stages. The replay explicitly uses native QPRES index5 (QGAME
is index4), QTEMP7 and species8--16; an initial mapping error was corrected
before interpreting results. `face_replay.json` retains the corrected results
under `build/regular-amr-native-transport-inputs-20261029/`.

Replaying the seven candidate root states with total-energy EOS gives initial
maximum pressure/temperature differences against native of1.01863e-10 Pa and
1.36424e-12 K. Replaying with the stored internal energy, as the current
regular-AMR transport path does, instead gives0.00546129482 Pa and
0.000105128528 K. At subsequent predictor/corrected stages both conversions
share the recorded state discrepancies. The same-input transport arithmetic
and initial EOS comparison point to old-state primitive preparation as the
next correction target; they are not a full propagated field-error proof.
`eos_comparison.json` and `held_eos_comparison.json` retain the values.

An isolated control under
`build/regular-amr-transport-total-energy-control-20261029/` removes the held
energy argument only from the regular-AMR physical transport flux call and
runs the complete original n32 serial field/inventory checker against the
saved native reference. It does not alter production source, tolerances,
input profile or held-energy evolution. This is a diagnostic nonrotating
control; adoption requires a scoped fix, regression tests and MPI evidence.

### 2026-10-29: total-energy transport control and scoped implementation

The isolated n32 control completes through300ns. With the identical input
hash and unchanged contract, violations decrease from18 to1. All100/200-ns
field/inventory gates pass. At100ns maximum x-momentum error decreases from
7.98455e-6 to9.30516e-9 kg m^-2 s^-1; scaled total-energy error decreases from
6.47612e-8 to2.33513e-10. The remaining300-ns y-momentum error is
1.43111670e-6 against1e-6 at fine cell(2,25,32), identical in value and
location to the previous run. This near-wall discrepancy is separate from
the diagnosed coarse/fine-interface maximum and remains a failed gate.
The complete control is therefore NOT passed.

The production implementation now selects total-energy transport primitives
only in nonrotating MOL sources, in both serial and MPI regular-AMR paths.
An optional internal API flag defaults to the existing held-energy behavior;
CTU and frame-present source calls retain that default. Stored energy remains
in donor interpolation and evolution. MPI collectively rejects inconsistent
primitive-energy settings before computing or changing outputs. Focused tests
exercise distinct stored/total energies, default held behavior, serial/MPI
rates and fluxes, and mismatch rejection. Their build/check status is recorded
under `build/regular-amr-transport-production-20261029/`; test source addition
is not itself a passing test. This change does not qualify n32, C2 convergence,
EB transport, rotating-frame parity, or physical validation.

The production n16 serial and MPI3 comparisons pass all four times with
unchanged reference gates; all eight state/budget CSVs are byte-identical
between serial and MPI. Five focused transport CTests pass5/5, including
stored/total-energy separation, serial/MPI face rates/fluxes, and rank-setting
mismatch rollback. Source hashes remained unchanged through the comparisons.
The first production n32 MPI3 run hits the3600-s command limit after writing
initial/100-ns outputs. Its four saved tables exactly match the isolated
serial control, but the run is incomplete, not a passing qualification. The
probe process was confirmed absent before starting a separate14400-s-limit
rerun in `build/regular-amr-transport-n32-mpi-extended-20261029/`. Numerical
inputs and acceptance limits remain unchanged; only the command timeout changes.

The stored-energy spray-runtime test targets are unavailable in this
configuration and are not counted as passed. The non-spray MOL conservation
and rotating-MOL checks are recorded separately in `mol-tests.log` under
`build/regular-amr-transport-production-20261029/`.

For the remaining y-momentum discrepancy, the200/300-ns top144 errors are all
in the two fine-cell columns closest to the wall. At200ns all144 candidate
values are unchanged by the transport correction. Paired observation-only
native/candidate runs now record fine cell(2,25,32) and its six neighbors
through all six fine intervals, under `regular-amr-wall-native-stages-20261029`
and `regular-amr-wall-candidate-stages-20261029`. They must reproduce the
original four native exports and eight corrected candidate state/budget CSVs
before stage comparison is accepted.

The native wall-near observer completes with all four exports unchanged.
At the target and six neighboring cells, accepted y-momentum increments
close against0.5*dt*(old RHS + corrected RHS) within6.93347e-17 kg m^-2 s^-1
across all six fine intervals. Directional face divergence reproduces the
saved y-momentum RHS with maximum residual3.71169e-8 kg m^-2 s^-2; this
rate residual must be retained when comparing small differences, rather than
calling the reconstruction bit-exact. The native analysis explicitly discards
42 pre-advance rows also labelled `final` by the generic `set_body_state`
observation hook, accepting final states only after the same interval/point's
corrected RHS. This observation-label ambiguity did not change solver outputs.
The ordering rule and raw rows are retained in
`build/regular-amr-wall-native-stages-20261029/native_analysis.json`.
The candidate and extended MPI runs remain live; their results are not yet
used to infer the wall discrepancy's cause.

The corrected candidate wall-stage observer completes with all eight baseline
state/budget CSVs unchanged. At fine cell(2,25,32), y-momentum differences
remain of order1e-13 through the first three fine intervals. In interval4
(150--200ns), the predictor-state hydro RHS difference jumps to
-4.72842135 kg m^-2 s^-2, while the transport RHS difference is
-1.32424e-6 in the same units. The predictor-state y momentum itself differs
by only2.22489e-13 kg m^-2 s^-1. Across all six intervals, accepted old/final
hydro RHS differences integrate to-1.43893346e-6 kg m^-2 s^-1; transport
contributes+7.81676500e-9. Their sum,-1.43111669978e-6, reproduces the final
failed-cell y-momentum error to approximately2.54e-16. This accounting
localizes the observed error growth to the hydro response at interval4,
without yet identifying the reconstruction/Riemann branch responsible.
`comparison.json`, `state_comparison.json`, and `accepted_error.json` are
retained under `build/regular-amr-wall-candidate-stages-20261029/`.

Paired200-ns stencil observations are now running under
`regular-amr-wall-candidate-stencil-20261029` and
`regular-amr-wall-native-stencil-20261029`. They record the four primitive
states and sound speeds around each face, plus directional face fluxes.
They shorten only the diagnostic endpoint to200ns and must match the
corresponding saved prefix outputs; full300-ns reference gates are unchanged.

The native200-ns stencil observer completes with three unchanged exports.
Replaying its actual four-state primitive/sound stencils through the current
Fortran MOL reconstruction/Riemann kernel gives maximum y-momentum face-flux
difference2.91038305e-11 Pa over72 faces (four fine intervals x three stages
x six faces). Native hydro flux is obtained by subtracting recorded diffusion
from total extensive flux and converting units/area. At the interval4
predictor, all six differences are at this scale or smaller; this cannot
account for the observed4.73 kg m^-2 s^-2 hydro-rate discrepancy. The replay
therefore directs the next comparison toward the solvers' actual stencil
inputs and their preparation, rather than assuming a different Riemann
formula. Unit mapping explicitly uses native QPRES5, QTEMP7, QFS8 and QC1.
`native_input_replay.json` and the replay input/source are retained under
`build/regular-amr-wall-native-stencil-20261029/`. The candidate observer is
still running; no conclusion about its stencil values is drawn yet.

The candidate200-ns stencil observer reproduces all six saved prefix CSVs;
all72 captured face reconstructions/fluxes replay exactly through the current
Fortran kernel. Direct native/candidate stencil comparison identifies a
physical wall ghost difference: at interval4 predictor, the x-low face's
outer stencil has candidate tangential velocity+6.62990141458 m/s versus
native-6.62990141459 m/s. Normal velocity already has the reflected sign.
The regular-AMR MOL hydro ghost fill currently omits `viscous_wall`, so the
shared sampler reflects only normal momentum even for no-slip/isothermal
walls; native no-slip reflects tangential momentum too.

A frozen-stencil control negating only the two tangential ghost velocities
reduces the event17 y-momentum RHS discrepancy from-4.72842133 to
-1.17874394e-5 kg m^-2 s^-2. This is a local intervention, not a completed
trajectory fix. `reflect_control.json` and `stencil_comparison.json` retain
all inputs/results under `regular-amr-wall-candidate-stencil-20261029`.
An isolated nonrotating-MOL `viscous_wall` control now runs the full original
n32 checker under `regular-amr-mol-wall-reflection-corrected-20261029`.
The preceding `regular-amr-mol-wall-reflection-control-20261029` attempt
failed compilation because of a missing Fortran continuation marker; it ran
no simulation and remains separate. Production hydro is not yet changed.

The extended pre-wall-fix n32 MPI3 run completes through300ns with all eight
CSV files byte-identical to the isolated serial total-energy-transport control.
The original checker still fails exactly the same300-ns y-momentum criterion
(1.43111670e-6 versus1e-6). Four times and unchanged source hashes are recorded
in `build/regular-amr-transport-n32-mpi-extended-20261029/verification.json`.
This establishes serial/MPI agreement for that incomplete trajectory, not
fixed-reference success. The earlier timeout remains an incomplete attempt.

A focused single-root MOL regression now creates a nonzero wall-normal
variation of tangential velocity and compares tree rates to a halo explicitly
filled with no-slip reflection. Normal-only reflection is a distinct negative
control. The current production hydro fails the intended exact-rate assertion;
the isolated corrected hydro object passes the same test. Existing auxiliary
MOL and rotating-MOL test modes also pass with that isolated object.
`build/regular-amr-mol-wall-regression-20261029/` preserves the failing log,
passing log, exact override-object hash and commands. The test source's new
`wall_mol` mode is not yet registered in CTest pending production adoption.
Full n32 wall-reflection reference verification is still running.


The isolated n32 wall-reflection trajectory completed through300ns and passed
all original field/inventory criteria at all four exported times
(`regular-amr-mol-wall-reflection-corrected-20261029/comparison/*/checks.json`,
status `passed`, empty violations). The corrected nonrotating MOL reflection
is now adopted in the serial and MPI hydro ghost fills; CTU and frame-present
selection remains unchanged. The wall regression is registered as
`unit_amr_patch_tree_no_slip_mol_3d`. MPI MOL spatial tests now include no-slip
and isothermal walls in addition to the existing boundary cases.

Under `build/regular-amr-mol-wall-production-20261029/`, the rebuilt targets
pass five focused tests (wall oracle, serial rotating MOL, MPI MOL spatial
2/4 ranks, MPI rotating MOL spatial2 ranks). Fixed serial/MPI applications
also build; the MPI runtime audit identifies ABI40. Provenance, project
contract, Python syntax and diff checks pass; exact commands and returncodes
are in `static-checks.json`. Production fixed-reference runs for n16 serial,
n16 MPI3 (500K wall) and n32 MPI3 (2050K wall) remain in progress. These
must be checked individually, including source hashes and CSV agreement;
completion of their wrapper is not itself a passing gate.

An additional wall-plus-transport MOL application restart check fails before
restart: the first checkpoint write is rejected. Both common-timestep and
recursive serial attempts fail. An artifact-only diagnostic checkpoint object
identifies the retained internal-energy average-down invariant as the rejecting
condition (`checkpoint-diagnostic/run.log`). Linking the pre-wall-fix hydro
object reproduces the same rejection (`old-hydro-run.log`), so this is not
introduced by the new tangential reflection. Eight diagnostic full-run controls
show periodic cases complete with transport on/off and common/recursive steps;
all four no-slip cases fail the same invariant, including transport disabled
(`checkpoint-controls/results.json`). This is a checkpoint consistency defect
requiring further investigation of covered-cell energy reset versus serialized
state. No invariant or tolerance has been relaxed. Neither MPI restart nor
configuration-change rejection is reached in these failed attempts; those
checks must not be reported as passed. C2 convergence and C3-C5 remain open.


Production wall-fix reference verification now completes: n16 serial and MPI3
at500K wall, and n32 MPI3 at2050K wall, all pass every original criterion at
four export times through300ns. Source hashes remain unchanged over the matrix.
The n16 serial/MPI outputs and n32 MPI/isolated corrected-serial outputs are
each byte-identical for all eight CSVs (`csv-verification.json` in
`regular-amr-mol-wall-production-20261029`). This is bounded trajectory
agreement, not a spatial-convergence or physical-validation completion.

Checkpoint diagnostics at the first5-ns wall step show a covered root-cell
held-energy restriction discrepancy of2.41598173e-5 J/m3 (common timestep) and
2.41598318e-5 J/m3 (recursive). The retained field instead matches
rho_E minus kinetic-energy density exactly. MOL's covered-cell energy reset
is applied after restriction: these two operations do not commute. Changing
the evolved or saved field to its fine-grid average would change continuation.

`regular-amr-mol-checkpoint-reset-20261029` records a separate checkpoint-only
change. Serial/MPI read and write interfaces now accept an explicit optional
`mol_energy_reset` mode, default false. Nonrotating MOL application callers
select it; MPI requires rank agreement and rejects a frame-present request.
The unchanged64-epsilon check accepts a held value equal to restriction or,
in that explicit mode, equal to total energy minus kinetic energy. Conserved
state restriction and held-energy EOS-temperature consistency remain required;
payload bytes and evolution are untouched. Unit coverage includes strict-mode
rejection, indexed/unindexed exact roundtrip, arbitrary held/conserved changes,
failed-write preservation and mismatched MPI flags.

The wall-plus-transport application test completes all common/recursive and
regrid/no-regrid variants in serial and MPI2-to1 restart. Gas CSVs and complete
checkpoints are byte-identical to uninterrupted runs; changed growth/LF
configuration is rejected. The test now reads both retained-energy V1/V2
headers and has a registered `application_amr_tree_mol_wall_restart` variant.
Its original negative missing-history fixture required an envelope-offset
correction; the first SDC history run failed in Python fixture construction,
not in the solver, and that failed log is preserved. Five other regression
tests passed (three serial/MPI checkpoint tests plus MOL default/history).
SDC history is being rerun after correcting that fixture.

A frozen copy of the wall-corrected serial reference probe now reruns n8/n16
at2050K wall in `regular-amr-corrected-spatial-20261029`; its analysis combines
those outputs with the verified n32 MPI outputs. It retains original33-point
profile sampling and dt100ns. State and increment differences are analyzed
separately using conservative factor-two volume restriction and SI L1/Linf;
this study is not an asymptotic space-time proof.

The corrected SDC missing-history fixture passes its full application rerun
(42.05s, `sdc-history-retest.log`). Final checkpoint evidence is six distinct
passing CTests across the initial/retest runs, plus the independent complete
wall-transport application test. The initial fixture-construction failure is
not counted as a solver pass. No full configuration matrix or installed-binary
qualification was performed for this change.


The corrected2050-K n8/n16 reruns complete and pass the original criteria.
Combining them with the verified n32 result, all three momentum increment L1
norms decrease at100/200/300ns in both solvers. Maximum z-momentum differences
remain nonmonotone: for PeleF,8/16 versus16/32 differences are8.84259289e-4
versus1.32766042e-3 kg m^-2 s^-1 at100ns, and3.10403792e-3 versus
3.58127245e-3 at200ns. Corresponding observed maximum-norm orders are
-0.58634481 and-0.20632609; fixed PeleC gives-0.58634479 and-0.20632593.
At300ns the order is only0.05573174. Thus fixing direct reference disagreement
does not establish spatial convergence. Input hashes, SI definitions and
state/increment results are in `regular-amr-corrected-spatial-20261029/`.
A matched dt50-ns n8/n16/n32 matrix with the same original profile, pinned
native binary and wall-corrected frozen candidate now runs under
`regular-amr-corrected-time50-20261029/` to measure temporal contamination.

`zmom_localization.json` and `localize.py` now locate every solver/time/pair
maximum of the corrected increment comparison. All12 maxima are in the first
fine comparison cell adjacent to x=5mm, outside the two-cell physical-wall
mask. For16/32, the comparison-grid fine cell has i=16 and x=4.84375mm;
y/z positions differ by symmetry/roundoff between solvers. At100/200/300ns,
the two-cell interface mask contains38.34/42.89/46.11% of the candidate
volume-weighted absolute increment difference; the wall mask contains
14.35/13.58/12.87%. These masks are defined using each comparison cell width
and are not a causal decomposition. The saved spatial field does not identify
the responsible interpolation, flux or synchronization stage.


The dt50-ns n8/n16/n32 matrix completes through300ns and passes every original
reference criterion at seven times for each grid. The unchanged frozen probe
and inputs are recorded under `regular-amr-corrected-time50-20261029/`.
The z-momentum increment Linf orders remain-0.51072521 and-0.18730351 at
100/200ns (PeleF), with fixed PeleC-0.51072416 and-0.18730361. For16/32,
halving dt changes the maximum differences by-0.4096%,-1.9554%,+1.2966%
at100/200/300ns. For8/16 the changes are+4.9497%,-0.6541%,-1.5693%.
Temporal contamination exists but this halving does not remove the shared
nonmonotonicity. Full field/scope/norm comparisons are in `time_sensitivity.json`;
this is not a completed temporal-order study.

The next diagnostic is a paired100-ns stage observer at the original dt100ns:
n16 fine cell(16,31,27), and its corresponding eight n32 fine cells
(31:32,61:62,53:54). It also observes the adjacent coarse donor cells.
Artifact-local split/MOL objects record old/predictor/corrected/accepted
states and hydro/transport RHS, all six directional hydro face fluxes,
reconstructed left/right pressure and four-cell primitive/sound stencils.
The production source is untouched. Both observers must reproduce their
saved initial/final state and budget CSVs exactly before interpreting traces.
`regular-amr-interface-stage-pair-20261029/` contains the code, source hashes,
commands, notification and prepared conservation/replay analysis.


The paired n16/n32 interface observers complete with all eight saved prefix
CSV files byte-identical. `analysis.json` closes each accepted z-momentum
increment to1.0842e-19 kg m^-2 s^-1; six-face hydro divergence matches the
recorded rate to5.4570e-12 kg m^-2 s^-2. All396 recorded reconstructions and
fluxes replay exactly through the same-input MOL kernel (`replay-result.json`).
At100ns, n16 minus restricted n32 z-momentum increment is-1.32766041884e-3:
hydro contributes-1.21873022178e-3 and transport-1.08930197064e-4.
The directional hydro contributions are x=-8.26672908e-9,
y=+1.70960230e-7,z=-1.21889291528e-3 kg m^-2 s^-1.
This is an exact update decomposition, not a numerical-policy intervention.

Initially the recorded hydro z-rate is zero and reconstructed z pressures are
uniform at101325Pa. The first fine predictor develops nonuniform pressure;
the mean hydro z-rates are-18852.1092 and-4732.2602 kg m^-2 s^-2 for n16 and
the eight n32 cells. A frozen-stencil first-order reconstruction replay reduces
the accumulated z-direction difference to-7.34914556e-4, while second order
reproduces-1.21889292e-3. It does not remove the existing pressure variation,
and cannot predict an evolved first-order trajectory or justify adoption.

`regular-amr-interface-pressure-origin-20261029` now expands the stage trace
to the five n16 z-stencil cells and their forty n32 fine counterparts, adding
all nine species reaction rates and reaction-energy rate. It repeats the
unchanged-prefix gate. Prepared EOS replay evaluates every subset of the
first predictor hydro/transport/reaction update, verifies the full predictor
closure, and assigns nonlinear pressure response symmetrically over the three
components. These frozen-state interventions are diagnostic; no production
numerical source or reference baseline is changed.


The expanded pressure observer reproduces all eight prefix CSVs. Its initial
three-term pressure attribution does NOT close the actual predictor: only
total energy has a material residual (up to16.2513 J/m3), giving up to4.97382Pa
pressure discrepancy. `pressure_origin.json` is retained as an incomplete
three-term diagnostic. The probe explicitly uses `allow_negative_energy=false,
eta2=2`; the resulting retained-energy reset must be included.

For all45 observed cells, replaying the unchanged production reset with the
recorded predictor held energy reproduces total energy and pressure exactly.
At the nine central cells, independently replaying UEINT face fluxes from the
captured primitive/sound stencils reproduces predictor held energy exactly
(`held_flux_closure.json`). Thus reset is not merely assigned from a pressure
residual. Four-term frozen additive attribution and its full-state closure are
in `reset_pressure_origin.json`; reset uses recorded held values outside the
nine independently flux-checked cells.

A stronger local operator control evaluates all32 subsets of directional
hydro, transport and reaction updates at those nine cells. Every subset
recomputes held advection from UEINT fluxes and applies the same reset policy.
The full predictor pressure is reproduced; directional conserved RHS sums
have a maximum absolute residual1.78814e-7 in mixed SI rate components
(large energy rates), with no pressure-closure error. Symmetric mean pressure
contributions are xhydro+23.58568,yhydro-14.33501,transport+0.0263817,
reaction-0.467838Pa at n16; n32 values are+18.27915,-14.25009,+0.0263803,
-0.467839Pa. zhydro is zero initially. These are frozen-input operator
attributions, not evolved source-removal runs.

At the interface x-high face, the two coarse-donor ghost velocities are
identical instead of varying across the fine stencil: n16 u=-0.91083024m/s;
n32 representative u=-0.39526698m/s. Replacing only their normal velocity by
the analytic INITIAL TG value, retaining rho,p,Y,sound and every other face,
reduces the first-predictor n16-minus-mean-n32 pressure difference from
5.22161703 to0.94898027Pa. Replacing tangential velocity too gives the same
pressure result. `ghost_velocity_control.json` records this bounded control;
analytic replacement is not a general interpolation method or adopted fix.

A native interpolation-only observer now records actual16-component coarse
donors for n16/n32 at the corresponding interface in
`regular-amr-native-interface-donors-20261029`. It must preserve four native
prefix exports. Timing callbacks are deliberately absent, so initial events
must be identified by matching saved initial conserved state, never by assuming
a callback index. The next check is which component determines the shared
linear limiter and whether it reproduces the observed flat velocity ghosts.


Native donor observation completes with four parsed prefix exports unchanged.
Initial donor event1/box4 matches all16 saved initial components exactly at
both grids, including temperature and retained energy. Periodic y-neighbor
indices must be wrapped when comparing the coarse export; the first analyzer
attempt omitted that wrap and failed before producing results. The corrected
analyzer records native shared slope factors[0,0,1] at n16 and[0,1,1] at n32.
Several components bind x slopes, not just roundoff in one trace species:
n16 density, transverse momentum, total/retained energy, temperature and major
species; n32 transverse momentum, total energy, temperature and H2. The
resulting two x-ghost velocities reconstruct n16(-0.91083024,-0.91083024,0)
and n32(-0.39526698,-1.45649572,0)m/s, matching candidate stencils to rounding.
This identifies a shared fixed-method interpolation behavior, not a remaining
PeleF/PeleC implementation disagreement.

An isolated conservative momentum-slope control now runs n8/n16/n32 through
100ns in `regular-amr-momentum-limiter-control-20261029`. It retains the shared
factor for all other components, but uses each momentum component's own
limited slope. Linear child offsets still sum to zero; runtime inventories
and positivity remain checks, not assumptions. Production source/defaults
and pinned native inputs are unchanged. The first wrapper completed n8 but
failed on a misspelled budget filename; its failure is retained. `resume.py`
verifies the already-completed n8 output and continues n16/n32 without rerunning
or overwriting it. This intentionally different method is only an attribution
experiment; neither smaller errors nor conservation alone can establish fixed-
reference qualification or justify adopting it.


The isolated independently limited momentum control completes100ns at
n8/n16/n32, preserving all six initial state/budget CSVs exactly. All16
increment L1 and Linf grid-pair differences decrease. z-momentum Linf order
changes from-0.58634 (shared) to1.19916; its16/32 maximum difference is
0.158713 times the shared-policy result. x-momentum is0.220776 times shared,
while y-momentum is1.005542 times shared: this is not uniform improvement.
All closed-domain element, mass/momentum/total-energy ledger and species-sum
checks pass the unchanged1e-10 threshold; largest listed residual is
3.33e-16. Species densities remain nonnegative.

The same control fails the original fixed-reference field limits: seven
fields at n8/n16 and nine at n32 at100ns (`reference_difference.json`).
This tradeoff is explicit. The default shared-policy production path retains
its passing reference evidence; the control is not adopted or counted toward
C2 fixed-reference completion. `regular-amr-momentum-limiter-full-20261029`
now repeats the identical isolated binary through300ns with12 required
initial/100ns prefix matches. Full state/increment norms, all-step conservation
and native-field differences are prepared as separate analyses.


The isolated momentum-limiter300-ns control completes all three grids and
reproduces all12 initial/100-ns prefix CSVs exactly. At100/200/300ns all16
state and increment L1/Linf grid-pair differences decrease; this is a bounded
three-grid observation, not an asymptotic or physical-validation claim.
z-momentum Linf orders are1.19916,1.18701,1.21612. Its16/32 maxima are
15.8713%,23.6707%,29.5819% of the unchanged shared-policy values. The y-momentum
maximum at100ns is0.5542% larger than shared; other observed increment maxima
are unchanged or smaller. All nine grid/time conservation checks pass the
original1e-10 thresholds, with maximum ledger/element residual3.33e-16 and
nonnegative species.

Original fixed-reference field limits still fail at every tested grid/time
(seven to ten observables per comparison). `reference_difference.json`
retains these violations, `policy_comparison.json` all component tradeoffs,
and `conservation.json` every inventory result. The isolated variant remains
unadopted. The shared default keeps its previously verified trajectory gates,
while its maximum-norm convergence and separate EB gates remain open.

The next independent C2 obligation is extrapolation/outflow x boundaries.
`regular-amr-outflow-spatial-20261029` runs the unchanged shared production
method at n4/n8/n16 through300ns with reaction, Simple transport and recursive
clocks. It uses the original checker including open-domain boundary supply,
species/element and energy ledgers. The frozen rebuilt binary and source-file
manifest are retained; wall/EB limitations are not cleared by this family.


The extrapolation/outflow n4/n8/n16 matrix completes through300ns, passing
all original field and candidate open-domain ledger gates. The362 recorded
source files and frozen probe hash still match (`verification.json`); maximum
reported relative ledger/element residual is2.71e-16. Native boundary flux
is not independently exported: the H2-consumption comparison uses candidate
boundary supply for both inventories and is not a native flux audit.

Spatial convergence remains open. Both solvers have nondecreasing increment
Linf for density, total/internal energy, temperature and major species at
100/200/300ns. Candidate total-energy increment Linf orders are-0.146324,
-0.037876,-0.012098. `maximum_locations.json` shows density, total-energy
and temperature maxima in the first cells adjoining physical x boundaries
for both grid pairs at all three times. x-momentum maxima are there at
100/200ns and move to the coarse/fine interface at300ns. This differs from
the corrected thermal-wall z-momentum interface localization.

`regular-amr-outflow-zero-control-20261029` runs a matched zero-initial-velocity
control in both solvers with unchanged transport, reaction, outflow, shared
interpolation and time step. It tests dependence on the initial velocity;
it does not replace the TG baseline or prove a boundary-compatibility cause.
Original field tolerances are retained, while this explicit control requires
exactly zero initial speed rather than the baseline minimum-speed condition.


The zero-initial-velocity outflow control completes all three grids through
300ns; original field/ledger tolerances pass and initial speed is exactly zero.
Source hashes remain unchanged. `baseline_comparison.json` compares both
solvers: at100ns candidate8/16 density, total-energy and temperature increment
Linf are0.0271%,0.2649%,0.7058% of the corresponding TG baseline values.
The maxima of these three fields and x-momentum now lie in the first fine
cell adjoining the x=5mm interface for both pairs at every output.
Transverse momenta remain zero. Density/energy/temperature maximum norms
still do not decrease; total-energy orders are-0.74191,-0.80998,-0.90290.
Thus the large physical-boundary error depends strongly on the initial flow,
but removing that flow does not establish all-field convergence.

`regular-amr-outflow-equilibrium-control-20261029` now compares zero-velocity
nonreacting hydro-only and nonreacting transport at the same three grids and
clocks. These are process-removal controls for the remaining interface error,
not replacements for the original reacting TG qualification case.

The hydro-only n4 control completes evolution but fails original density,
x-momentum and energy field gates, native/candidate open-inventory agreement,
and candidate x-momentum relative ledger gates. The first driver correctly
stops on that failure; its record is retained. A diagnostic continuation
reuses those complete outputs and preserves every failed criterion while
collecting the remaining grids. This is not a passing qualification matrix;
the known separate pinned nonreacting-MOL2 limitation must be considered
before interpreting these process-removal controls.

All six nonreacting process-removal controls now complete their outputs;
none passes all original gates. Hydro-only n8/n16 fail only candidate
x-momentum relative ledger gates. Across these six controls the absolute
x-momentum inventory residual is at most3.20e-22 kg m/s; normalization by
near-zero momentum can produce a relative residual up to0.0167. This
quantifies the failure scale and does not change or waive the criterion
(`momentum_ledger_audit.json`).

At100ns candidate8/16 hydro-only increment maxima are1.64e-13 kg/m3
(density),3.69e-7 J/m3 (total energy), and zero temperature. With transport
they become8.69e-10 kg/m3,0.00195843 J/m3 and0.000382862 K; all three peak
at the first fine interface cell and have nondecreasing grid-pair maxima.
The reacting transport control is substantially larger again. These are
process-removal observations, not additive causal decompositions. An
artifact-only stage/face-pressure capture is being prepared for the
nonreacting transport interface cell and its eight fine counterparts;
its CSV identity must pass before its diagnostics are used.

The corrected artifact stage capture reproduces all eight initial/100ns
state and budget CSVs exactly at n8/n16. The initial capture omitted
`--reaction-off`, failed prefix identity, and remains excluded. Nonreacting
MOL has old/predictor RHS stages only; the analyzer was corrected to that
branch rather than assuming two reacting correctors.

For the fine-interface cell and its eight counterparts, the100ns x-momentum
increment difference is9.2591980e-7 kg/m2/s: hydro contributes9.2881048e-7
and direct momentum transport contributes-2.8906788e-9. All hydro contribution
is x-directed. Accepted-update closure is7.45e-23 and directional flux/RHS
closure3.56e-15. Initial face pressure is uniform to rounding; the first
predictor produces101325.01935–101325.02854 Pa at n8 and
101325.01872–101325.02638 Pa at n16. These observations locate the pressure
response after transport but do not yet separate heat/species flux causes.
Artifacts: `regular-amr-outflow-transport-stages-20261029/analysis.json` and
`corrected_record.json`. A separate directional transport-face capture is
running; its first preparation failed before solver launch and is preserved.

The directional transport capture completes with eight exact prefix CSVs.
Recorded face divergences reproduce all14 transport rates exactly, and those
rates match the split-stage records exactly. At the first fine predictor,
all state components except total energy equal the raw hydro+transport update.
The held energy is unchanged; the reset removes0.680158 J/m3 at the n8
interface cell (n16 paired mean raw transport increment0.555762 J/m3).

A same-state production-reset replay separates species and total-energy
transport increments while retaining observed held energy and hydro. Full
replay reproduces observed predictor energy and pressure exactly. With species
transport, pressure rises0.0285429722 Pa at n8 and0.0233278552 Pa averaged
over n16 counterparts. Without species transport both remain uniform to
3e-11 Pa, regardless of inclusion of transport total-energy increment.
Thus this local pressure response is due to species transport under the
fixed held-energy reset, not direct retained heat input. This is a frozen
first-predictor attribution, not a physical-validation claim. Pinned native
`Diffterm.H` places diffusive heat/species enthalpy flux in UEDEN, not UEINT;
`Advance.cpp` updates and resets the predictor before its second RHS.

Artifacts are `regular-amr-outflow-transport-faces-20261029/analysis.json`,
`pressure_replay.json` and `verification.json`. An isolated first-old-fine
transport-ghost control now mirrors the two x-high ghost layers using the
initial cosine profile symmetry about x=L/2. It changes no hydro ghosts,
later interpolation or production source; it is only a local causal control.

The initial-symmetry transport-ghost control completes with all four initial
state/budget CSVs exact. The first-predictor pressure mismatch changes from
0.00521511694 Pa to-1.66223617e-6 Pa. Old H2 transport rates become
0.160634711 and0.160629482 kg/m3/s at the n8 cell and n16 paired mean,
versus0.221843283 and0.181294380 with original coarse interpolation.
This intervention confirms a large first-predictor contribution from the
initial coarse/fine transport ghosts. At100ns the paired x-momentum
increment difference falls only from9.25919804e-7 to6.64494594e-7 kg/m2/s
(28.2341% reduction). Later ghost fills are unchanged and this does not
explain all subsequent error or establish convergence. The control is
specific to the initial symmetric profile, not a general interpolation
method, and is not adopted. All production numerical defaults and original
acceptance limits remain unchanged. See `baseline_comparison.json` and
`verification.json` under `regular-amr-outflow-transport-mirror-control-20261029`.

The original reacting TG outflow family now has a dt50ns n4/n8/n16
comparison running in `regular-amr-outflow-time50-20261029`. It retains
the original profile, shared interpolation, native binary and tolerances.
`compare_time.py` compares spatial norms at common100/200/300ns with the
completed dt100ns family; two time steps alone are not an order estimate.
The362 recorded source files match the previous default probe manifest.

A separate boundary-capability audit confirms that isothermal slip needs an
explicit velocity/thermal distinction: the current `isothermal_wall` label
means no-slip, while configuration wall_temperature defaults to300K even on
adiabatic slip faces. Treating positive/present wall temperature as a new
activation condition would change existing cases silently and is rejected
as an implementation approach. No boundary or numerical default changed.

The dt50ns outflow n4/n8/n16 matrix completes, passing all original
seven-output field and candidate ledger gates;362 source hashes remain
unchanged. At common100/200/300ns, total-energy increment Linf orders
remain-0.141933,-0.035475,-0.008813. The8/16 maxima change by only
-0.5002%,-0.3822%,-0.4550% relative to dt100ns. Density and temperature
likewise retain nondecreasing maxima at physical-boundary cells. Thus
halving time step does not remove the observed spatial nonmonotonicity.
At200ns the x-momentum worst cell shifts to the coarse/fine interface;
the other stated maxima remain at the physical boundary. Full224 common-
time state/increment comparisons are in `time_sensitivity.json`.

`regular-amr-outflow-n32-time50-20261029` now adds the original n32 family
at50ns, retaining the frozen33-point initial profile and all original
criteria. Its prepared analyzer combines it with completed n8/n16 outputs
for8/16 versus16/32 differences; this tests the finer-grid trend without
changing the problem or retrospectively dropping the coarser failures.

While n32 runs, `initial_profile_regularity.json` records the frozen table's
one-sided derivatives. At x=L/2 the temperature slope jumps from
-3074.355 to+3074.355 K/m and H2 mass-fraction slope from-0.122974
to+0.122974 per metre. Native `prob.H` uses linear table interpolation.
The endpoint slopes are also nonzero, unlike the underlying analytic cosine.
This is a regularity/compatibility limitation of the fixed input, not an
excuse for failed gates or proof of the measured-error cause. The frozen
profile remains untouched.

A fresh MPI3 n16/dt50ns outflow comparison is independently running under
`regular-amr-outflow-mpi3-time50-20261029`, rebuilding the existing clean
GNU/OpenMPI target, using the completed native capture and requiring all14
serial state/budget CSVs to match exactly. It does not substitute for n32
convergence, installed-binary or full configuration-matrix qualification.

The rebuilt MPI3 n16/dt50ns outflow comparison completes successfully,
passing the original seven-time field and candidate ledger gates against
the captured pinned-native trajectory. All14 state/budget CSVs are byte-
identical to serial, and all recorded source hashes remain unchanged.
The independent receipt is `regular-amr-outflow-mpi3-time50-20261029/verification.json`.
This extends bounded CPU/MPI reference evidence; n32 spatial convergence,
other C2 families and C3–C5 remain open.

The original n32/dt50ns reacting outflow serial case completes through300ns,
passing all seven-time field and candidate-ledger gates. The recorded362
source files remain unchanged and the binary matches the checks receipt
(`regular-amr-outflow-n32-time50-20261029/verification.json`). Its combined
8/16/32 spatial analysis is being evaluated separately. An MPI3 run of this
same case now reuses the completed native reference and requires all14
serial state/budget CSVs to match; it is not yet a completed result.

The8/16/32 dt50ns spatial analysis completes416 solver/time/scope/field
rows. All state and increment L1 differences decrease in both solvers at
all six positive times. At250/300ns all16 state and increment Linf
differences also decrease. Earlier times remain incomplete: density,
energy, temperature and major-species increment maxima increase slightly
through200ns, while z-momentum maxima increase at50/100ns near the
coarse/fine interface. Candidate total-energy increment orders at
100/200/300ns are-0.00723356,-0.00553174,+0.0644569. The earlier density/
energy/temperature maxima are physical-boundary-adjacent.

`nondecreasing_norms.json` retains all64 nondecreasing combinations across
both solvers; `convergence_summary.json` prevents selecting only the passing
late endpoints as whole-trajectory convergence. These finer-grid results
show improvement at later times but do not close C2 convergence.

## 2026-10-30: n32 outflow MPI closure

The n32/dt50ns MPI3 comparison completes through300ns. All seven original
field and candidate-ledger gates pass; all14 state/budget CSVs independently
rehash identically to the completed serial run. All362 recorded source
files remain unchanged and the probe hash matches the checks record.
Evidence: `build/regular-amr-outflow-mpi3-n32-time50-20261029/verification.json`.
This is bounded fixed-reference and serial/MPI agreement. The early-time
nonmonotone maxima in the8/16/32 spatial study remain unresolved; neither
physical validation nor full release/configuration qualification follows.

The outflow matrix now has serial n4/n8/n16/n32 dt50ns reference evidence,
MPI3 n16/n32 exact agreement and separate dt100ns sensitivity. Further
identical repetitions are not the next step. Remaining work is the recorded
early-time convergence problem and distinct C2 boundary/EB obligations;
the isothermal-slip implementation must explicitly separate thermal and
velocity conditions while preserving existing adiabatic defaults.

## 2026-10-30: explicit isothermal-slip transport foundation

The physical transport and wall-timestep kernels now accept an optional
`isothermal_slip_faces(6)` mask, independent of existing wall-temperature
values and velocity boundary labels. Only explicitly marked slip faces
receive thermal flux; velocity ghosts retain normal-only reflection.
Unmarked slip faces stay adiabatic even when wall_temperature is supplied.
The new focused test covers six faces, heating/cooling, Fourier flux,
zero tangential shear/species exchange, thermal timestep and invalid masks.

Compilation and three focused CTests are pending under
`build/isothermal-slip-kernel-20261030/`; saved pre-turn source and test
copies plus `change.patch` distinguish this change from existing work.
This is kernel implementation only: AMR/MPI propagation, rank-consistency
checks, application/probe configuration and fixed-native comparison remain
to be completed before reporting isothermal-slip support.

The explicit thermal-slip kernel build completes and all three focused
CTests pass: reactive boundary transport, constant transport and rotating
transport. The first includes six-face hot/cold analytical heat flux,
no tangential shear/species exchange, thermal timestep, unchanged default
adiabatic behavior and rejection of invalid masks/disabled conduction.
`isothermal-slip-kernel-20261030/verification.json` records commands and
source identity; this is focused kernel evidence, not reference parity.

The mask now also reaches serial AMR SDC transport and wall timestep APIs,
and MPI SDC/transport and timestep APIs. Each patch clips the global mask
to actual physical slip faces, excluding coarse/fine and interior faces.
MPI compares all six mask bits collectively before computation and rejects
rank disagreement without changing caller outputs. Three added/extended
AMR tests are building under `isothermal-slip-amr-20261030`; results remain
pending. Split-integrator forwarding, probe/configuration and fixed-native
comparison are still unfinished. Saved pre-change copies isolate this work
from the previously dirty checkout.

## 2026-10-30: thermal-slip AMR checks and integrator verification

The three focused AMR transport tests pass after correcting the test oracle:
serial AMR, MPI4 SDC transport and MPI4 constant SDC transport. The initial
oracle incorrectly required identical slip/no-slip rates for a moving,
nonuniform state, whose ghost states differ by design. The corrected test
uses a stationary uniform state, retains exact rate/flux/timestep comparisons
and invalid-mask rollback assertions, and does not relax tolerances.
The original failure and corrected results are preserved in
`build/isothermal-slip-amr-20261030/{tests.log,corrected-tests.log,corrected-verification.json}`.

The explicit mask now propagates through serial/MPI common and recursive
SDC integrators and the fixed-reference probe. MPI collectively validates
mask agreement and thermal prerequisites before entering the stage.
Verification is running under `build/isothermal-slip-integrator-20261030/`:
five focused tests, n4 500 K slip-wall common/recursive serial/MPI3 reference
comparisons, sixteen expected identical state/budget CSVs and a no-slip
regression. These are planned checks, not yet reported passes.

Application configuration and checkpoint binding remain pending. Existing
`wall_temperature=300` defaults must not enable thermal slip implicitly.
The shared reactive configuration serves additional unsupported applications;
an accepted but ignored thermal-slip option would be incorrect. Explicit
capability validation and restart rejection for changed thermal masks are
required before claiming application support. C2 convergence remains open.

The integrator verification completes: all five focused CTests and all five
n4 reference cases pass original gates. The four thermal-slip cases cover
500 K x walls through300 ns with common/recursive clocks and serial/MPI3;
the fifth preserves the no-slip baseline. Sixteen serial/MPI state/budget
CSVs were independently compared byte-for-byte after completion.
Evidence: `build/isothermal-slip-integrator-20261030/verification.json`,
`tests.log` and each case's `checks.json`. This bounded comparison does not
establish spatial convergence or physical validation.

## 2026-10-30: thermal-slip application and restart wiring

Regular serial/MPI tree gas drivers now read `isothermal_slip_faces(6)` from
`amr_tree_timestepping`, defaulting to allfalse. Enabled bits require slip
velocity boundaries and SDC/MOL with enabled thermal conduction; wall
values remain in `reactive_3d.wall_temperature`. Both common and recursive
automatic timestep paths and SDC/MOL stages receive the explicit mask.
Enabled masks append a checkpoint-context record; absent/allfalse keeps the
previous context representation. Changed masks must reject restart.

The first build failed on a Fortran line-length limit; continuation lines
correct it and both fixed applications build. The first application fixture
used PLM for SDC and was rejected by existing validation; the fixture now
uses PPM for SDC and PLM for MOL. Original failure logs are retained under
`build/isothermal-slip-application-20261030/`. The retry runs fixed/selected
serial/MPI, common/recursive SDC/MOL, thermal/default behavior, rank-changing
restart and changed-mask/temperature refusal. Application checks remain
pending until the retry record completes. No numerical tolerance changed.

Application testing exposes a pre-existing nonrotating SDC checkpoint issue:
the moving-state fixture refuses a checkpoint with thermal slip enabled and
also with the mask disabled. Stationary SDC and moving MOL controls write
successfully. An artifact-only checkpoint diagnostic reports covered-cell
held energy127885.70652404375 J/m3 versus restricted child held energy
127885.70650161128 J/m3, while total energy minus kinetic energy equals the
stored held value exactly (original tolerance1.8173651955641468e-9 J/m3).
The SDC recovery routine calls the same dual-energy reset as MOL, but the
application enabled explicit reset-aware checkpoint validation only for MOL.

Both gas drivers now request that existing identity check for nonrotating
SDC as well. No checkpoint payload, numerical update, tolerance or low-level
strict default changes. Arbitrary held-energy changes and conserved-state
restriction errors remain rejected by the existing validator. This separate
checkpoint correction is saved as `sdc_checkpoint_change.patch`; controls
and diagnostic object/source/log remain under `checkpoint-controls/` and
`checkpoint-diagnostic/`. A fresh notification-managed verification under
`sdc-checkpoint-fix/` includes fixed/selected application matrices and four
transport-subcycle/checkpoint unit tests. Results remain pending.

The SDC checkpoint fix verification passes: fixed and selected applications
each complete68 runs covering common/recursive SDC/MOL, serial/MPI2 full
runs, MPI2-to-serial/MPI1/MPI3 restart, exact serial/MPI checkpoint bytes,
changed mask/temperature refusal, disabled-conduction/wrong-face refusal,
and explicit/omitted adiabatic defaults. Four CTests pass: serial/MPI3
transport-subcycle and checkpoint validation, including corruption refusal.
The full source manifest and eight application matrices were independently
rechecked in `sdc-checkpoint-fix/independent_receipt.json`. This is application
and restart evidence; it is separate from the n4 fixed-reference comparison.

The reusable checker is now registered as fixed/selected CTests. The MPI
subcycle test additionally compares six-face slip thermal timestep limits
against existing isothermal walls and serial results, and requires collective
rollback for rank-inconsistent masks and disabled conduction (including an
empty-owner rank). Verification of these additions and six n8/n16 reference
cases is pending under `build/isothermal-slip-regression-20261030/`.

## 2026-10-30: thermal-slip finer-grid reference and wall-scale diagnosis

All three registered tests pass: fixed/selected thermal-slip application
matrices (68 runs each) and the extended six-face MPI3 recursive timestep
check. All six500 K reference comparisons pass original gates: n8 common
and recursive serial/MPI3 and n16 recursive serial/MPI3 through300ns.
All24 state/budget CSV pairs match exactly;363 recorded source files remain
unchanged in the independent receipt. Evidence is under
`build/isothermal-slip-regression-20261030/`.

The4/8/16 recursive spatial audit retains all224 state/increment rows and
finds40 with at least one negative observed order across both solvers.
Candidate energy-increment Linf orders at100/200/300ns are
-0.7159185007,-0.9274500219,-0.9769229288; native behavior agrees closely.
Energy and temperature increment maxima lie in the first fine cell at the
physical wall. This is not spatial convergence, despite fixed-reference
agreement. No tolerance or boundary condition is changed for acceptance.

An artifact-only probe calls the production wall heat evaluator on the
imported initial states. Initial total wall powers for n4/n8/n16 are
-41.9416789,-84.3311450,-168.8688202 W. The local scale
`sqrt(k_wall/(rho*cv_cell)*300ns)/dx_fine` is0.01155,0.02313,0.04628.
Thus even n16 resolves substantially less than one cell across this local
conduction scale. Actual boundary energy divided by initial Fourier power
times time is0.998551..0.999964 across the nine positive-time samples.
This supports an unresolved wall-temperature-jump contribution; the local
scale is not an exact reacting-flow solution or an exemption from convergence.
`build/isothermal-slip-wall-scale-20261030/` retains evaluator source, build
recipe, per-face data, hashes and the original V3-parser failure plus correction.

A2050 K diagnostic control is running on4/8/16 grids to reduce the initial
wall temperature jump and isolate remaining interface effects. It does not
replace the500 K case or close its convergence requirement.

The2050 K thermal-slip diagnostic completes and passes all original reference
gates on4/8/16 grids. Its224-row spatial audit reduces the number of rows
with a negative order from40 to8. Both energy and temperature increment
L1/Linf differences decrease at all three positive times; their maxima move
from the physical wall to the coarse/fine interface. At300ns the8/16 energy
and temperature increment maxima are37.44% and12.36% of the500 K values.
This intervention supports a wall-temperature-jump contribution to the cold
case; it does not qualify that case or establish global convergence.

The warm case retains x-momentum increment Linf orders-0.1881217 and
-0.0399932 at100/200ns (positive0.0749113 at300ns), peaking at n8 fine
cell(8,10,7). z-momentum increment L1 is also slightly nonmonotone at300ns.
An artifact-only capture under `build/isothermal-slip-interface-stages-20261030/`
records old/predictor/corrector/accepted rates, reconstructed face pressures
and directional hydro fluxes at that cell and its eight n16 children.
It must reproduce all12 original prefix state/budget CSVs through200ns
before the stage decomposition is used. No production numerical change is
made for this diagnostic.

The thermal-slip interface capture completes with all12 prefix state/budget
CSVs unchanged. Accepted momentum updates close within5.14e-17 kg m^-2 s^-1;
directional hydro flux divergence closes within5.83e-11 kg m^-2 s^-2 and
accumulated updates reproduce saved outputs within2.17e-19 kg m^-2 s^-1.
At100ns the n8-minus-restricted-n16 x-momentum increment is-0.01636082938,
comprising hydro-0.01606272595 and transport-0.00029810342 kg m^-2 s^-1.
The hydro difference is dominated by the x direction (-0.01606480259).
Pressure/flux and stencil records retain every observed stage and direction.

All792 captured faces replay with exactly matching reconstructed pressures
and fluxes. Replacing only first-old fine ghost primitives beyond the interface
by their initial midplane reflection reduces the initial x-momentum rate
mismatch from-192110.384835 to-588.752695 kg m^-2 s^-2 (99.6935%).
This frozen-input intervention does not integrate a new trajectory and worsens
the y-momentum rate mismatch; it is not adopted as a general interpolation.

A same-input call to the actual pinned AMReX slope kernel reproduces all96
slopes for the two current coarse donors exactly. Reconstructed initial ghost
density/velocity values match all19 captured entries exactly. Both x-direction
density extrema force the shared scale to zero, although x-momentum's own
slope ratio is1. This confirms shared reference-policy behavior for the
thermal-slip diagnostic, rather than a candidate interpolation arithmetic bug.
Evidence: `analysis.json`, `replay-result.json`, `mirror-replay.json`,
`thermal-native-verified.json` and `independent_receipt.json` under the capture
directory. This evidence does not close convergence.

A frozen-probe n32 warm-wall serial/MPI3 comparison is now running under
`build/isothermal-slip-warm-n32-20261030/` to examine the8/16/32 sequence
with unchanged reference gates. The500 K wall convergence requirement remains
separate and unresolved.

## 2026-10-30: n32 warm thermal-slip comparison

The n32/2050 K recursive serial and MPI3 comparisons complete through300ns
with original fixed-reference field/inventory gates passing. All eight
state/budget CSVs independently compare exactly and all363 source-manifest
files remain unchanged. The serial/MPI wrapper elapsed times are3734/6768s;
these include native reference/export work on a shared host and are not a
solver parallel-speedup measurement. Evidence:
`build/isothermal-slip-warm-n32-20261030/verification.json` and
`independent_receipt.json`.

The8/16/32 audit contains224 state/increment rows. x/y-momentum,
energy and temperature L1/Linf differences decrease at all positive times.
z-momentum Linf remains nonmonotone at100/200ns for state and increment in
both solvers (eight rows). Candidate increment orders are-0.5863495318,
-0.2063378449,+0.0563575504 at100/200/300ns. L1 decreases. Its maximum
16/32 difference occurs at n16 fine cell(16,31,27), adjacent to the interface.
The passing300ns endpoint does not qualify the full trajectory.

An artifact-only observer now captures that cell and its eight n32 children
through100ns, with all three momentum components, stage pressure and six
face fluxes. `build/isothermal-slip-transverse-stages-20261030/` requires eight
original prefix CSVs to remain exact before its analysis is accepted. This
separates the remaining transverse maximum from the previously diagnosed
4/8/16 x-momentum maximum. C2 and500 K wall convergence remain open.

## 2026-10-30: transverse stage pressure and normal-velocity control

The transverse observer completes with all eight original prefix CSVs exact.
An independent receipt rechecks those pairs, all four captured-source hashes,
and 396 reconstructed face pressures/fluxes with zero replay error. Accepted
momentum updates close within 1.25e-17 kg m^-2 s^-1, directional flux divergence
within 2.92e-11 kg m^-2 s^-2, and saved outputs within 1.09e-19 kg m^-2 s^-1.
At 100 ns the n16-minus-restricted-n32 z-momentum increment difference is
-0.001327660419 kg m^-2 s^-1: hydro contributes -0.001218730222 and transport
-0.0001089301971. The z-direction hydro contribution is -0.001218892915;
the x/y-direction terms are much smaller. This is an observed update budget,
not a claim that transverse transport or synchronization is generally absent.

First-predictor pressure reconstruction reproduces the captured baseline
exactly, including the production internal-energy reset. Algebraic subset
accounting attributes mean pressure changes of 9.250672477 Pa (n16) and
4.029057321 Pa (n32 children) to hydro; transport contributes approximately
0.02638 Pa and imported old reaction approximately -0.46784 Pa on both grids.
These mean-pressure contributions are not a pressure-gradient decomposition.

A frozen first-predictor control changes only the initial x-ghost normal
velocity to its midplane-reflected value, retaining captured transport and
old-reaction rates. Conservative fluxes, internal-energy fluxes and the reset
are recomputed with production routines. The unmodified replay has exactly
zero hydro-rate, held-energy and pressure error. The mean n16-minus-n32
predictor pressure difference changes from 5.221617031 Pa to 0.0001399300963 Pa
(99.9973% reduction). This identifies a normal-velocity-to-pressure pathway
at the sampled interface intersection. It is a frozen diagnostic, not an
evolved solution, general interpolation fix or replacement reference policy.

Evidence is under `build/isothermal-slip-transverse-stages-20261030/`:
`analysis.json`, `replay-result.json`, `pressure-partition.json`,
`normal-velocity-control.json` and `independent_receipt.json`. The separate
`flat-pressure-control.json` is only a frozen face-flux sensitivity: fixing
face pressures while retaining captured density/sound speed is not an
EOS-consistent evolution and is not adopted.

Time-step sensitivity is running under
`build/isothermal-slip-time-refinement-20261030/`: n8 at 50/25 ns, n16 at 50 ns
in serial/MPI3, and n32 at 50 ns, all through 300 ns with unchanged reference
gates. Spatial and temporal analysis scripts retain physical-time matching,
SI conversion, composite-cell masking and volume weights. Pending comparisons
are not counted as passes. This 2050 K diagnostic does not replace the
unresolved 500 K wall requirement or close C2.

The n8 50/25 ns cases have now passed the original reference gates. Comparing
100/50/25 ns at common 100/200/300 ns times gives 96 solver/time/field rows;
all L1/Linf self-differences decrease. Candidate z-momentum Linf temporal
orders are 1.03752, 1.04568 and 1.05204. This approximately first-order
behavior on one mesh is not a second-order or spatial-convergence claim.
`temporal_convergence.json` and `temporal_receipt.json` retain the analysis,
input hashes and independent recheck of all 363 source-manifest files.
The n16/n32 50 ns spatial comparison remains pending.

## 2026-10-30: completed warm-wall time refinement

All five planned cases finish with exit code zero and unchanged original
reference gates passing. Independent verification rehashes all 363 source
files and confirms 14 exact n16 serial/MPI3 state/budget CSV pairs.
The dt50 ns spatial audit contains 416 rows. Twenty have negative Linf order:
only z-momentum state/increment, both solvers, at 50/100/150/200/250 ns.
Their L1 orders remain positive and every n16 maximum is on the fine-side
coarse/fine interface (fine i=16). Candidate increment Linf orders at
100/200/300 ns are -0.510726725, -0.187310169 and +0.014301655.

Halving dt changes the 16/32 z-increment maximum by -0.410%, -1.955% and
+1.296% at those common times. Thus temporal refinement does not remove
the early-time spatial nonmonotonicity. The maximizing cell can move along
the interface; at 100 ns it is (16,18,11), rather than the dt100 ns location
(16,31,27). Evidence: `independent_receipt.json`, `spatial_convergence.json`,
`spatial_time_sensitivity.json` and `spatial_summary.json` under
`build/isothermal-slip-time-refinement-20261030/`. These are self-difference
and reference-compatibility results, not physical validation.

A separate n8 common-timestep 100/50/25 ns comparison is running under
`build/isothermal-slip-common-time-20261030/` to distinguish dependence on
recursive subcycling from behavior shared by both timestep modes. It retains
the original reference acceptance and the 2050 K diagnostic boundary.

The three common-timestep cases now pass their original reference gates;
all 363 source hashes remain unchanged. All 96 temporal comparison rows
have decreasing L1/Linf differences. Candidate z-momentum Linf orders are
1.04113, 1.04736 and 1.05299 at 100/200/300 ns, compared with recursive
1.03752, 1.04568 and 1.05204. The six common/recursive ratios of the two
temporal Linf differences lie between 0.99450 and 1.00154. Thus removing
subcycling does not restore second-order behavior in this control.

The region audit contains 216 rows for density, momentum, energy and
temperature. Regions partition the n8 composite grid into interface cells
(root i=5 or fine i=8), physical-wall cells (root i=8 or fine i=1), and
interior cells. L1 is normalized by each region's volume. z-momentum temporal
maxima occur on the coarse side of the interface in both timestep modes;
this spatial association alone does not identify a defective operator.
Evidence: `temporal_convergence.json`, `mode_comparison.json`,
`temporal_regions.json` and `independent_receipt.json` under
`build/isothermal-slip-common-time-20261030/`.

An artifact-only synchronization observer under
`build/isothermal-slip-sync-stages-20261030/` now records root states before
local stepping, before/after flux synchronization and after state recovery.
It covers both timestep modes at all three dt values and must reproduce all
96 original state/budget CSVs before its decomposition is used. Production
numerics and reference acceptance remain unchanged.

## 2026-10-30: coarse-interface synchronization budget

The synchronization capture finishes with all 96 original CSV pairs exact.
Recorded interval-start states and post-recovery states match saved outputs
exactly; the accumulated update budget closes within 8.68e-19 across all
conserved fields. At the 12 z-momentum temporal maxima (two modes, three
times, two dt pairs), the direct synchronization contribution is at most
0.04508% of the total temporal difference. Momentum recovery contributes
zero. The local pre-synchronization update dominates this signed budget.
For example, common stepping at 100 ns with dt100-minus-dt50 gives total
-0.0002496664660, local -0.0002497078742 and synchronization
+4.140813659e-8 kg m^-2 s^-1 at root cell (5,4,3).

This does not exclude indirect effects of earlier synchronization on later
local updates. It locates the observed difference and motivates recording
the local hydro/transport rates and directional face pressures/fluxes.
Evidence: `build/isothermal-slip-sync-stages-20261030/record.json` and
`analysis.json`. The initial analyzer incorrectly assumed that common MOL
used a separate stepping engine; source inspection confirmed both modes
dispatch to the recursive engine, with subcycling disabled for common mode.
The analyzer now verifies the actual mode from the input and the observed
hook separately; the original failure is preserved in
`analysis_initial_failure.txt`. No numerical capture was changed.

The next artifact-only capture is
`build/isothermal-slip-root-stages-20261030/`, covering all 64 coarse-side
interface cells, both modes and all three dt values. Its acceptance requires
the same 96 exact original CSVs and agreement with the independently recorded
local-update budget before interpreting hydro/transport contributions.

## 2026-10-30: local pressure forcing and delayed synchronization feedback

The root-stage capture reproduces all 96 original CSV pairs. z-momentum
local updates close within 5.83e-19 kg m^-2 s^-1, directional face divergence
within 7.28e-12 kg m^-2 s^-2, and the independent synchronization-capture
local budget matches exactly. At common dt100-minus-dt50, 100 ns, root
(5,4,3), local difference -0.0002497078742 comprises hydro -0.0002484618599
and transport -0.000001246014276 kg m^-2 s^-1. The z-direction hydro term
is -0.0002484872706; x/y terms are much smaller.

All 4536 selected face stencils replay exactly in pressure and flux. A frozen
control replacing only z-face pressures by 101325 Pa reduces the corresponding
12 temporal hydro differences by 96.33--99.34%. This retains captured density,
velocity and sound speed and is not an EOS-consistent evolved intervention.

A separate production-EOS replay evaluates 10752 recorded conservative
states around synchronization. Its interval-start pressures and post-recovery
pressures match the actual current/next old-stage pressures exactly. At the
same representative common-step cell and time, pressure changes from
101324.0077480 Pa before synchronization to 101316.2571756 Pa after conservative
reflux, then 101318.1193909 Pa after recovery. Thus small direct z-momentum
reflux does not imply small pressure feedback: the net pressure change is
-5.888357142 Pa and is seen by the next hydro evaluation.

A second frozen control adds this recorded post-recovery-minus-pre-sync
pressure change to each corrector z-face stencil cell, leaving old-stage
pressures and all other Riemann inputs unchanged. The 12 temporal hydro
differences shrink by 82.97--98.15%; at common 100 ns, dt100-minus-dt50,
the difference changes from -0.0002484618599 to +0.000005497314850 kg m^-2 s^-1.
This identifies sensitivity to delayed synchronization pressure feedback.
It uses future recorded corrections and is not an evolved algorithm, a
general correction or a passing fixed-reference implementation. A numerical
change must still preserve conservation/restart and pass the original
reference gates; spatial convergence and the 500 K boundary remain unresolved.

Evidence under `build/isothermal-slip-root-stages-20261030/`:
`analysis.json`, `pressure-control.json`, `sync-pressure.json`,
`sync-pressure-control.json` and `independent_receipt.json`.
The initial standalone EOS helper omitted the probe's explicit EOS policy;
its 1.822e-5 Pa residual failed the original check. Configuring the same
policy gives exact pressure agreement without relaxing that check. The
failure is preserved in `sync_eos_initial_failure.txt`. Replaying the two
earlier first-predictor reset diagnostics with the explicit policy produces
byte-identical outputs (81 and 36 rows); see `prior-reset-policy-audit.json`.

## 2026-10-30: isolated composite-stage synchronization trial

An artifact-only common-MOL prototype under
`build/isothermal-slip-stage-sync-prototype-20261030/` evaluates all levels
at each stage, synchronizes hydro/transport and internal-energy spatial rates
with the existing conservative reflux routine, and uses those rates in the
predictor/corrector. The final restriction/reset retains its energy ledger.
Unlike the preceding frozen sensitivity, it does not use future recorded
states. `prototype.patch` and `design.md` describe the exact change and its
limited serial common-step scope. Production source remains unchanged.

The first n8/dt100 ns case completes but fails 21 original field gates across
three positive times: density, three momenta, total/internal energy and
temperature. At 300 ns the z-momentum difference is 0.002253862418 kg m^-2 s^-1
against the unchanged 1e-6 limit, and temperature differs by 0.2297186403 K
against 0.001 K. Gas ledger residuals are below 3.40e-16 and energy residuals
below 7.41e-17 in this bounded case. Conservation does not compensate for
reference failure; the prototype is not adopted. dt50/25 ns runs remain
pending to measure temporal self-convergence separately. This experiment
changes stage/level ordering as well as synchronization and does not isolate
a single operator or qualify MPI, recursive stepping, restart or physics.

The three prototype runs now finish. Original reference gate failures number
21/42/84 at dt100/50/25 ns, always the same seven fields. Maximum gas and
energy ledger residuals are 3.40e-16 and 8.19e-17; H/O/N relative inventory
drifts remain below 2.79e-16. Candidate z-momentum temporal Linf orders are
1.58898, 1.32538 and 1.39238 at 100/200/300 ns, with maxima now on the fine
side of the interface. This improves on the production common-step orders
but does not establish second-order convergence or fixed-reference agreement.
At dt25 ns the endpoint z error is still 0.001441851927 kg m^-2 s^-1 and
temperature error is 0.2312734879 K. `independent_receipt.json` rechecks all
363 original source files and records the failed gates and separate budgets.

Inspection found an additional difference introduced by bypassing the
production dispatch: the whole-tree path reaverages the covered old energy,
whereas production MOL preserves its already-reset old state. Accordingly,
these results cannot be attributed solely to stage flux synchronization.
Three companion trials complete a two-by-two comparison: stagewise/end
synchronization crossed with preserved/reaveraged covered old state. They are
under `isothermal-slip-whole-stage-control-20261030/`,
`isothermal-slip-stage-sync-old-state-20261030/` and
`isothermal-slip-whole-stage-old-state-20261030/` in `build/`.
End-sync controls also synchronize held energy at the final restriction.
`compare_variants.py` in the original prototype directory verifies terminal
records, source hashes, temporal inputs and original gates before producing
the comparison. All variants remain isolated and unadopted while pending.

The four-way comparison now completes. Every variant finishes all three
trajectories and every original reference comparison fails; none is adopted.
The following z-momentum temporal Linf orders use common physical times:

| Isolated whole-tree variant | 100 ns | 200 ns | 300 ns |
| --- | ---: | ---: | ---: |
| Stage synchronization, reaveraged old state | 1.58898 | 1.32538 | 1.39238 |
| End synchronization, reaveraged old state | 1.04177 | 1.04804 | 1.05334 |
| Stage synchronization, preserved old state | 1.58937 | 1.32652 | 1.39464 |
| End synchronization, preserved old state | 1.04171 | 1.04803 | 1.05334 |

Preserving the production old-state invariant therefore does not remove the
observed contrast between synchronization timings. This supports the temporal
feedback diagnosis, but does not demonstrate second-order convergence or
repair spatial maxima. Whole-tree ordering remains distinct from the
production level-ordered path. `variant_comparison.json` records all 12
original gate results, field errors, conservation/element budgets, 363-file
source checks for each variant, binary/patch hashes and hashed temporal
inputs. `compare_variants.py` requires terminal successful experiment wrappers;
individual reference-checker failures remain explicit. Production solver
code is unchanged and C2 remains open.

## 2026-10-30: planned coupled 16/32/64 refinement

The reference checker and serial/MPI probe admission lists now also accept
64 root cells per direction; numerical algorithms and acceptance limits are
unchanged. The new run under `build/isothermal-slip-coupled64-20261030/`
first requires n8 serial/MPI3 reference passes and 16 exact comparisons with
the archived state/budget CSVs before starting n64/dt25 ns through 300 ns.

The reason for this refinement is recorded before the result: captured
initial n32 fine-interface sound speeds are 1128.9155--1130.2025 m/s,
corresponding to a fine-cell crossing time of approximately 138.25 ns.
The unresolved 100 ns maximum therefore occurs before one local cell
crossing. n64 halves that crossing-time scale. This is a local scale estimate,
not a domainwide stability bound or proof that refinement will converge.

The frozen comparison uses n16/dt100 ns, n32/dt50 ns and n64/dt25 ns, with
the same recursive 2050 K thermal-slip case. Analyze all 16 fields, both
solvers, state and evolution-increment norms at 100/200/300 ns, retaining
all original intermediate reference gates. This measures coupled mesh/time
refinement, not fixed-dt spatial order. It does not replace the failed
earlier spatial maxima, the 500 K boundary or EB obligations. `design.md`,
`acoustic_scale.json`, source manifests and notification records retain
the rationale, resource preflight and pending verification. No n64 result
is counted as passed while the job is running.

## 2026-10-31: completed n64 trajectory, original reference failure

Both native and candidate trajectories reach 300 ns at n64/dt25 ns, with
13 compared output times. The comparison exits 1 because 17 original field
gates fail, not because the trajectories are incomplete. x momentum fails
at steps 2/6/7/8/9/10/11/12; y momentum at 10/12; density at 11; total and
internal energy at 10/11/12. The largest x-momentum difference is
7.508831494e-6 kg m^-2 s^-1 at 250 ns, versus the unchanged 1e-6 limit.
All eight failing x-momentum maxima lie on fine i=64, the coarse/fine
interface. This is a regular-grid case without cut cells.

The maximum gas-ledger relative residual is 1.16e-16. All 363 source hashes
remain unchanged, and the n8 serial/MPI3 regression comparisons pass with
16 exact prior state/budget CSV pairs. `independent_receipt.json` retains
the completed times, original violations and exact identities. Conservation
and regression passes do not qualify the failed n64 reference comparison.
The n64 wrapper took 78170.4 s including native calculation, export, candidate
calculation and comparison; this is not a solver speed or scaling result.

The original success-gated coupled analyzer remains unchanged. A separately
named `analyze_coupled_diagnostic.py` records the n64 reference failure while
analyzing the completed trajectories' self-convergence. Its result remains
qualification=false and cannot replace the original field gates. A separate
saved-endpoint donor-stencil audit targets coarse (33,53,61), associated with
the first failing fine cell (64,105,121) at 50 ns. It is a screening diagnostic;
actual in-step interpolation and causality are not inferred from endpoints.

### Saved-endpoint donor conditioning and next stage capture

At coarse donor (33,53,61), the saved 25 ns endpoint has a shared y-direction
slope factor of 0.0121511432724471 in the fixed reference and
0.00176709814744512 in the candidate. Both are controlled by z momentum.
The corresponding left differences are 2.3049690124e-10 and
3.3445254876e-11 kg m^-2 s^-1. The inferred fine ghost (65,105,121)
y-velocity difference is 1.0233212332e-3 m/s. Substituting only the reference
z-momentum stencil into the candidate screening input restores the reference
factors and reduces that velocity difference to 3.73e-14 m/s. This is an
algebraic intervention on saved endpoints, not an evolved solver correction.

`donor-native-replay.json` independently records 288 exact slopes from the
actual pinned AMReX kernel: both endpoint stencils at 0/25/50 ns, 16 components
and three directions. Thus the screening arithmetic matches that kernel;
it does not prove the actual intermediate halo values or the origin of the
small z-momentum difference. Input and native binary hashes are retained.
Neither limiter changes nor tolerance changes are adopted.

The next observer-only run is
`build/isothermal-slip-donor-stages64-20261031/`. It repeats only the first
two 25 ns steps, recording root (32:34,52:54,60:62), fine
(63:64,104:106,120:122), local stage states, reconstructed pressures,
directional face fluxes, and root pre/post-reflux and recovery states.
All six state/budget CSVs must match the completed trajectory byte for byte.
The covered root column is included for donor tracing, not composite norms.
Native in-step comparison remains required. This capture is pending and is
not counted as a passed check or a causal result.

### Completed coupled refinement diagnostic

The separately gated diagnostic analysis finishes successfully with 224 rows
(all 16 fields, both solvers, state at 0/100/200/300 ns and increments at
100/200/300 ns). No row has a negative observed L1 or Linf order for coupled
16/32/64 root grids and 100/50/25 ns timesteps. Candidate z-momentum increment
Linf orders are 0.205993, 0.430282 and 0.610083 at 100/200/300 ns; the fixed
reference gives 0.205993, 0.430282 and 0.610084. Corresponding candidate L1
orders are 1.518842, 1.538494 and 1.577112. All these z-momentum maximum
locations remain on the fine coarse/fine interface.

This is evidence of decreasing differences along this coupled refinement
sequence, with low maximum-norm order. It is not fixed-dt spatial convergence,
second-order convergence, physical validation, or a replacement for the 17
failed n64 original reference gates. Earlier nonmonotonic spatial sequences
remain part of the evidence. The full definitions, maxima, original gate
statuses and input hashes are in `coupled_diagnostic_convergence.json`.

### Completed candidate n64 stage capture

The observer-only first-two-step capture completes in approximately 3 h 4 min
and all six original state/budget CSVs match byte for byte. Independent
`analyze.py` checks 378 local momentum budgets and 162 root synchronization
component records. Maximum accepted-update closure residual is 5.55e-17
kg m^-2 s^-1; maximum directional face-divergence residual is 5.83e-11
kg m^-2 s^-2. These checks establish observation consistency, not agreement
with native intermediate states.

At root donor (33,53,61), the first25 ns z-momentum change is
8.2416973945e-8 kg m^-2 s^-1: local evolution contributes
2.7869516532e-9 and reflux contributes7.9630022292e-8, with zero momentum
change from recovery. The local contribution separates into hydro
7.4843490047e-10 and transport2.0385167528e-9. Thus endpoint donor momentum
is dominated by reflux in this candidate run; the origin of the much smaller
native/candidate difference is not yet established.

An isolated native observer under
`build/isothermal-slip-native-donor64-20261031/` is prepared to capture the
same first50 ns, including actual filled states, primitives, sound speeds,
reconstructed face pressures, hydro fluxes and pre-hydro transport fluxes.
It overrides only local copies of Diffusion.cpp and MOL.cpp objects and
retains the original executable and source tree. All three exported native
state CSVs must match the original run before interpreting the observations.
This native capture is pending, not counted as passed.

The first native observer attempt builds but aborts before stepping because
the diagnostic runner omitted required `profile.csv`. Its failed record and
log are retained. The corrected runner in
`build/isothermal-slip-native-donor64-r2-20261031/` copies the unchanged
profile along with the other original inputs; no numerical change is made.
The native face comparison script is prepared but remains unexecuted until
all three output identity checks pass.

### Paired native stage observations and first ghost amplification

The corrected native observer completes and all three exported state CSVs
match the original native trajectory exactly. `face_comparison.json` compares
2268 candidate/native face records; 21210 repeated native tile-face records
are exact before deduplication. Within each level, hydro calls are matched
chronologically. Native CGS pressure and momentum flux are converted to SI;
pre-hydro transport flux is divided by its native face area before conversion.

Captured root and fine initial states are identical. At the first fine
predictor, interior momentum differences remain at approximately1e-14
kg m^-2 s^-1, but the adjacent ghost (65,104,121) has candidate-minus-native
velocity differences (+6.82023122e-6,+3.87807513e-4,-8.15507186e-12) m/s.
This is an actual in-step filled ghost, not the previous saved-endpoint
screening. At the corresponding x-high face of fine(64,104,121), x-momentum
hydro flux differs by-2.99773747e-4 Pa while reconstructed left/right
pressures differ by only5.01e-8 and2.01e-7 Pa. The target(64,105,121)
has smaller differences at this first predictor; later stages propagate the
neighboring discrepancy. Root corrector chemistry/species and pressure also
differ slightly, so this observation alone does not attribute all error to
one upstream mechanism.

`state_comparison.json`, `stencil_comparison.json` and
`neighbor_stencil_comparison.json` retain the paired values. The neighbor
stencil audit is explicitly x-direction only: its initially attempted full
stencil includes y cells outside the captured native region and is not
claimed as a complete neighborhood audit.

Two observer-only follow-ups capture the actual coarse interpolation inputs:
`build/isothermal-slip-actual-donors-20261031/` stops deliberately after the
first fine predictor hydro evaluation and requires byte-identical prefixes
of the earlier stage captures plus both initial state/budget CSVs;
`build/isothermal-slip-native-actual-donors-20261031/` observes the fixed
AMReX MF interpolater inputs and requires three exact original native outputs.
Neither changes the interpolation arithmetic. These follow-ups are pending.

### Actual interpolation input replay (native trajectory still pending)

The candidate early-stop observer completes as intended. Its two initial
CSV outputs are exact and its four stage/face/synchronization files are
byte-identical prefixes of the completed50 ns capture. It stops after the
first fine predictor hydro call; this is not a completed50 ns trajectory.

At the actual12.5 ns interpolation for donor(33,52,61), ghost(65,104,121),
the z-direction shared factor is0.889407506501707 in PeleF versus
0.884411824720706 in the native CGS input. Both are controlled by H2O2.
The native H2O2 centered slope is3.9587943232e-24 g cm^-3, versus
3.9708159071e-21 kg m^-3 in the candidate; the left differences are
1.7506022555e-24 g cm^-3 and1.7658367373e-21 kg m^-3. This actual first
amplification differs from the saved25 ns endpoint screen, where z momentum
controlled a different donor's y limiter. Endpoint screening alone therefore
would not have identified the first in-step controller.

Replay in each solver's original units reproduces the previously verified
actual ghost velocities (native within4e-24 m/s, candidate exactly). The
actual pinned AMReX slope kernel gives96 exact same-input slopes. Thus the
observed shared-limiter amplification is explained by its inputs in this
sample, not a different implementation of the slope arithmetic.

Converting the native stencil to SI *before* interpolation changes the
predicted y velocity by-1.7454e-7 m/s because very small neighbor differences
are sensitive to conversion rounding. This diagnostic artifact is retained,
not silently treated as the native answer. Substituting only native H2O2
into the candidate SI stencil reduces its y-velocity difference from
3.8781e-4 to-1.7454e-7 m/s, essentially the conversion-rounding residual.
This is an algebraic input substitution, not an evolved fix or a reason to
change acceptance limits or suppress trace species.

`build/isothermal-slip-actual-donors-20261031/actual-conditioning.json`
records the selected actual input row, original-unit kernel replay and the
unit-conversion sensitivity. Native input capture was still running when
this analysis was written; its final three-output identity checks remain
pending. The earlier native stage observer already passed those identity
checks and supplies the ghost values used in this replay. Next trace the
small H2O2 donor-input discrepancy through the root reaction/source update,
keeping it separate from the now-observed interpolation amplification.

### Native actual-input capture finalized; upstream reactor replay

The native MF interpolation observer finishes with exit0 and three exact
original state CSVs. An independent receipt rehashes each captured/original
pair, and rerunning `analyze_actual.py` records the terminal
`captured_outputs_identical` status. The96 exact same-input slope comparisons
and the H2O2-controlled actual-ghost amplification remain unchanged.

The next isolated capture is
`build/isothermal-slip-donor-reactor64-20261031/`: the first25 ns native step,
recording root donor(33,52,61) and six neighbors before/after initialization
and both root reaction calls. Two native output CSVs must remain exact.
The prepared replay links the current regular-AMR reference-probe libraries,
uses its reduced-coordinate `advance_forced_reactor_cvode` and the same
configured NASA7 EOS, and applies the production density-scaled absolute
tolerance. The scaling is explicit in `reactive_3d_mod.F90`: the reactor
receives input atol multiplied by old density. It must not be inferred from
comparing raw input tolerance numbers across SI and CGS.

This same-input replay will distinguish reactor-generated H2O2 differences
from preceding transport/source differences. It is pending. No tolerance,
kinetic convention, limiter, source arithmetic or original field gate is
changed, and the earlier EB full-coordinate/tolerance experiments are not
substituted for this regular-probe configuration.

### Identical-input root reactor comparison

The native first25 ns reactor capture finishes with both original state CSVs
exact. The current regular-probe reactor replay completes21 samples (seven
cells, initialization and both root reaction calls). With the same native
rhoY/T/energy/source/duration, the reduced-coordinate production integrator
has maximum H2O2 endpoint differences of1.3973e-16,1.4166e-16 and1.3964e-16
kg m^-3 for the three events. The replayed half-time H2O2 donor values differ
from actual candidate interpolation inputs by at most1.26e-21 kg m^-3.
Thus same-input reactor behavior explains the leading endpoint discrepancy;
remaining source/rounding differences still matter to much smaller spatial
increments.

Replacing only H2O2 in the candidate stencil by these same-input replayed
half-time endpoints gives a z common factor0.8896678090, compared with
actual candidate0.8894075065 and native0.8844118247. The resulting ghost
y-velocity difference is4.08014e-4 m/s, versus actual3.87808e-4 m/s.
This is an algebraic diagnostic, not an evolved correction.

A second21-sample replay changes only to the existing full-coordinate CVODE
entry, retaining current linked mechanism, configured EOS and effective
absolute tolerance. H2O2 endpoint differences increase to approximately
1.1103e-15 kg m^-3, but the inferred spatial factor improves to0.8850923011
and H2O2-only ghost velocity difference falls to5.28244e-5 m/s. Endpoint
absolute error and limiter-sensitive spatial differences therefore give
different rankings; neither alone justifies adoption.

The isolated evolved discriminator
`build/isothermal-slip-full-coordinate-prefix64-20261031/` switches only the
reference probe's callback to the existing full-coordinate integrator and
stops after the first fine predictor. Initial state/budget CSVs must remain
exact; changed subsequent stages are compared against the original native
observer. This is pending and does not qualify the original50/300 ns gates,
MPI/restart or convergence. Production sources and acceptance limits remain
unchanged. Reactor records and conditioning algebra are in
`build/isothermal-slip-donor-reactor64-20261031/`.

### Evolved full-coordinate prefix and bounded gate experiment

The full-coordinate prefix completes with two exact initial state/budget
CSVs. The actual ghost y-velocity discrepancy falls from3.87808e-4 to
4.41391e-5 m/s (approximately88.6% reduction); its shared z factor becomes
0.8849804178, versus native0.8844118247. Across the captured first fine
predictor faces, maximum momentum-flux discrepancy falls from2.99774e-4
to3.41485e-5 Pa. Maximum reconstructed-pressure discrepancy there falls
from1.26865e-6 to1.90979e-7 Pa. In contrast, root corrector maximum pressure
discrepancy grows from3.07860e-7 to3.73067e-7 Pa. All five captured hydro
stages are retained in `stage_summary.json`; no universal improvement or
full-trajectory qualification is inferred.

The next experiment, `build/isothermal-slip-full-coordinate-trial64-20261031/`,
first requires original n8 serial/MPI3 gates through300 ns and eight exact
state/budget CSV pairs. Only if those pass does it run n64/dt25 ns through
50 ns, the first original x-momentum gate failure time. Its isolated checker
retains every original acceptance threshold and check, but explicitly sets
the prefix to two steps and qualification=false; the original300 ns checker
and contracts remain untouched. Both serial and MPI experimental probes
change only the callback to the existing full-coordinate integrator. All363
production-source manifest entries are unchanged. No production adoption,
MPI64/restart, convergence or physical-validation result is claimed.

### Full-coordinate small-grid gates passed; n64 launcher recovery

The original n8/300 ns serial and MPI3 reference comparisons both pass;
eight state/budget CSV pairs match exactly, and all recorded comparison steps
and budgets are identical across ranks. `small-grid-comparison.json` retains
all64 field/time comparisons against the previous reduced-coordinate run,
including improvements and regressions. This is bounded reference/MPI
evidence, not full-coordinate production adoption or C2 closure.

The subsequent n64 diagnostic checker aborts before any calculation because
its artifact directory does not contain the tools module import path
(`pelephysics_transport_import`). The failed wrapper and log are retained.
`resume_n64.py` revalidates both n8 check records, all eight output pairs,
experimental binaries and363 source hashes, and restarts only n64 with an
explicit `PYTHONPATH` pointing to the repository tools. Checker source,
acceptance values, binary and numerical inputs are unchanged. The retry is
recorded separately in `retry_record.json` and `retry_notification/` and is
pending. No interrupted or unstarted n64 trial is counted as passed.

### Full-coordinate n64/50 ns trial fails original field limits

The corrected launcher completes both25 ns steps and all three output times.
The bounded reference comparison fails at50 ns in x momentum, total energy
and internal energy. Maximum x-momentum difference is1.6181571944e-6
kg m^-2 s^-1 at fine(64,8,97), versus the unchanged1e-6 limit and the
reduced-coordinate baseline1.2931953390e-6 at fine(64,105,121).
Scaled total/internal energy errors are1.1549807101e-8 and1.1550128677e-8,
respectively, above the unchanged1e-8 limits. The favorable first-predictor
local observation therefore does not persist as a global field improvement.
This full-coordinate candidate is not adopted.

The experiment wrapper exits0 because the bounded trial and its comparisons
completed; its recorded checker exits1 with all three failures retained.
All363 source hashes remain unchanged. `bounded-final-receipt.json` records
source identity, field failures and conservation residuals. The retry wrapper
took approximately6 h32 min including reference verification/export,
candidate calculation and comparison; this is not a solver scaling result.

A saved-data audit now compares all16384 finei64 interface cells at50 ns,
including mean, median,95th percentile, maximum and momentum-limit counts
for both reduced and full-coordinate variants. It also labels saved25 ns
old/new maximum donor stencils explicitly as endpoint screening. This
requires no new CFD trajectory and is pending; a single improving local
sample is not a basis for another solver variant.

### Whole-interface audit and fixed-state kinetics comparison

The saved-data audit completes for all16384 finei64 cells at50 ns. For
x momentum, reduced/full absolute-error mean is1.22878e-8/1.31486e-8,
median6.76275e-10/3.69506e-10,95th percentile4.08773e-8/6.74929e-8,
and maximum1.29320e-6/1.61816e-6 kg m^-2 s^-1. The number above1e-6 rises
from8 to12. y/z momentum means,95th percentiles and maxima also worsen,
although medians decrease. This supports rejection of the full-coordinate
candidate beyond the single reported maximum.

Saved25 ns endpoint donor factors at the new maximum's parent(33,4,49)
are[0,1,0] for native, reduced and full variants. At the old parent
(33,53,61), full improves the y factor relative to reduced. These are
endpoint screens only and cannot explain all intervening stage behavior.
`native_endpoint_screen_original_units.json` independently reads original
CGS data without a SI roundtrip and is the authority for native endpoint
factors; `interface_audit.json` records the whole-plane errors and hashes.

A42-state pre/post native-reactor audit under
`build/isothermal-slip-donor-rates64-20261031/` compares chemical production
and the9x9 density Jacobian at identical T/rhoY using current regular-probe
libraries. Maximum norm-scaled differences are3.11358e-6 for production and
2.42975e-5 for the Jacobian; maximum nonzero H2O2 production relative
difference is5.74652e-6. This excludes forced-energy recovery and the
independent temperature Jacobian column.

The existing pinned native activation-basis diagnostic is verified to change
only23 activation-energy entries relative to the canonical bundle, with
species/thermo/transport and all other reaction parameters unchanged. Its
Fortran module is freshly generated from the verified imported bundle and
linked against the same current core. The42-state production/Jacobian
norm-scaled differences fall to3.29061e-10 and1.18914e-10. No canonical
mechanism, generated production source or golden hash is overwritten.

On the21 identical native reactor inputs, this activation-basis diagnostic
with full coordinates reduces maximum H2O2 endpoint difference to4.54303e-17
kg m^-3 (reduced coordinates give1.20471e-15). Yet the H2O2-only half-time
stencil substitution still gives z factor0.8850182632 and ghost y-velocity
difference4.70769e-5 m/s; the native factor is0.8844118247. Smaller absolute
chemical error therefore still does not establish the required spatial-ratio
agreement. These algebraic diagnostics do not justify another full trajectory
or adoption. Reaction-constant conversion, integrator differences and
ill-conditioned shared limiting remain distinct evidence classes.

### Matched forced RHS and CVODE policy discrimination

`build/isothermal-slip-forced-rhs64-20261031/` directly invokes pinned
PelePhysics `fKernelSpec` with the e-reactor type on42 native pre/post
states. The candidate uses the current full RHS with the same affine energy,
external species/energy sources, elapsed time and temperature seed. Maximum
inverse-EOS temperature difference is9.09495e-13 K. With canonical kinetics,
maximum temperature-rate relative difference is2.55998e-6 and species-RHS
norm-scaled difference is3.19975e-6. With the previously verified activation
basis, these fall to4.10875e-10 and3.37944e-10 respectively, without changing
the recovered-temperature bound. This checks the actual forced RHS, beyond
the earlier fixed-temperature rate oracle, but not internal CVODE trial
sequences or the complete approximate Jacobian.

The reduced/full integrator comparison also changes solver policies:
reduced coordinates use CVODE's default maximum BDF order with explicit
maximum step and stop time, whereas the existing full-coordinate entry uses
the pinned maximum order2 and ordinary dense output. To distinguish these
from coordinate effects, `build/isothermal-slip-cvode-policy64-20261031/`
compares three source-derived policies on the same21 native reactor inputs,
retaining canonical kinetics, EOS, effective tolerance and reduced variables.
No production source changes or arbitrary tolerance sweeps are made.

| Reduced-coordinate policy | Steps per captured case | Maximum H2O2 endpoint difference (kg m^-3) | H2O2-only ghost y residual (m/s) |
| --- | ---: | ---: | ---: |
| Order2, existing stop policy | 122 | 9.79621e-16 | 1.27901e-4 |
| Default order, native dense-output policy | 30 | 1.35843e-16 | 4.00383e-4 |
| Order2, native dense-output policy | 122 | 9.79824e-16 | 1.38113e-4 |

The original reduced replay uses30 steps, while the full-coordinate replay
uses116. These are local integrator counts, not a performance benchmark.
The table's ghost residuals are algebraic H2O2-only substitutions into the
saved donor stencil, not evolved trajectory results. Matching the order
improves that local spatial ratio but does not remove its discrepancy;
changing dense-output policy alone has little effect. These controls do
not justify another long trajectory or production adoption. The rejected
full-coordinate n64 gate result remains authoritative.

A further *fixed-reference-policy* check uses the native input's constant
species density atol1e-18 g cm^-3, converted to1e-15 kg m^-3, on the21
activation-compatible full-coordinate replay inputs. This is one
source-derived policy, not a tolerance sweep or a relaxed acceptance gate.
The native typical-value option is disabled by default and absent from the
case input. The candidate scalar temperature atol remains1e-15 K instead
of native1e-18 K, so full vector-weight equivalence is not claimed.
Maximum H2O2 endpoint difference falls to1.49469e-20 kg m^-3, yet the
H2O2-only spatial factor becomes0.8860188522 and ghost y residual
1.24751e-4 m/s. This stronger absolute agreement still fails to recover
the tiny neighbor-difference ratio. The result is retained in
`native-density-atol-comparison.json`; no tolerance or integrator policy is
adopted and no additional CFD trajectory is launched.

### 2026-10-31 forced-RHS density roundtrip discrimination

The artifact `build/isothermal-slip-density-roundtrip-20261031/` removes
arithmetic roundtrips in isolated copies of `constant_volume_reactor_mod.F90`.
Both copies retain the original RHS call for its input validation and temperature
recovery. `rates_only` recomputes molar production from the original
`density*y/molecular_weight` concentrations, without dividing the resulting
mass rate by density and multiplying it back. `direct_density` also computes
concentrations directly from rhoY. These are diagnostic overrides, not a
production implementation or an efficiency measurement.

Both variants completed 42 matched native pre/post RHS evaluations and 21
same-input reduced-coordinate integrations using the canonical mechanism,
existing EOS and unchanged tolerances. Total accepted-step counts were 630
for each 21-case set. The maximum species-RHS absolute changes relative to
the baseline were 1.3877787807814457e-17 and 5.551115123125783e-17 kg/m3/s;
both retain the native norm-scaled discrepancy 3.1997511278338604e-6.
Maximum H2O2 endpoint discrepancies remain 1.4166370363885772e-16 and
1.4166370350961075e-16 kg/m3. In the previously defined H2O2-only algebraic
donor substitution, z limiter factors are 0.8896829929913804 and
0.8896714226202874, with transverse ghost-velocity differences
4.0919312885456804e-4 and 4.0829493776861625e-4 m/s. This does not resolve
the amplification; no production change or long CFD run is justified by
this control. It is not evidence about all internal trial states, evolved
trajectory agreement, MPI/restart, or physical validation.

`run.py` records exact compile/link commands, input/source/executable/output
SHA-256 values and per-case outputs in `comparison.json`. The production
reactor source is checked unchanged before/after the experiment. The next
remaining discriminator is the reaction/EOS arithmetic and internal solver
trajectory at the seven actual donors, rather than this density roundtrip.

### 2026-10-31 native reactor replay and internal trial-state audit

`build/isothermal-slip-native-reactor-replay-20261031/` adds a standalone
single-cell CPU box-API replay using the actual fixed native ReactorCvode
objects and native CGS inputs from the seven donor cells at three reaction
events. All 21 species-density vectors, temperatures and energy densities
match the previously observed CFD post-reaction states exactly. The standalone
wrapper initializes both AMReX and its SUNDIALS context, as the native main
does. An initial missing-object link failure and a missing-context startup
failure are retained separately; neither is counted as a completed test.

An observer-only ReactorCvode object records every RHS evaluation. All 21
endpoint output records remain byte-identical with and without the observer.
There are 2,667 RHS evaluations, 2,478 accepted internal steps and zero error-test
failures; 67 evaluated trial states have at least one negative species density.
The observer and replay preserve native order2, dense analytical Jacobian,
rtol1e-10, species atol1e-18 g/cm3 and temperature atol1e-18 K.

Current Fortran full-RHS evaluation at each of those native trial states
completes successfully. Against the same-input native RHS oracle, recovered
temperature differences are zero for all 2,667 states. The canonical mechanism
retains maximum species norm-scaled discrepancy3.206005680735004e-6 and
thermal-rate relative discrepancy2.560520181694672e-6. The existing verified
activation-basis diagnostic reduces these to3.390610416370696e-10 and
4.1119480394103294e-10, respectively. Maximum absolute H2O2 rate differences
are5.656924499143623e-8 and8.430756093247282e-16 kg/m3/s. Native CGS-to-SI-to-CGS
recomputation has species norm-scaled noise at most8.137368729523851e-16.
This extends the earlier endpoint RHS audit to the actual native internal
trials; it is not an assertion that different integrators visit identical states.

An isolated full-coordinate candidate with the existing activation-compatible
mechanism now uses native-derived vector absolute tolerances:1e-15 kg/m3 for
species and1e-18 K for temperature. It completes all21 inputs with2,478 accepted
steps. An observer-only candidate trace preserves endpoint bytes and records
2,667 RHS calls, with all21 per-case evaluation counts matching native.
However, even the first positive trial times differ in floating-point arithmetic;
for case0 they are3.6327455412951863e-19 s native versus
3.6327455412951593e-19 s candidate. Initial species states agree after unit
conversion; initial temperature-coordinate differences reach9.094947017729282e-13 K.
Ordinal trace comparisons are descriptive, not same-time state comparisons.

Matching vector weights does not resolve the donor sensitivity: maximum
H2O2 endpoint difference remains1.494692102428593e-20 kg/m3, while the
H2O2-only algebraic donor substitution gives z factor0.8859985247175095 and
transverse ghost-velocity difference1.2317321377253165e-4 m/s. This is not an
evolved CFD result, and no solver tolerance, mechanism or production algorithm
is changed. C2 remains open. Further work should discriminate the residual
reaction arithmetic and its spatial differences, rather than assume missing
temperature weights or grossly different adaptive step counts explain them.

Artifacts: `comparison.json`, `trace/receipt.json`, `trace-comparison.json`,
`vector_atol/comparison.json`, and `candidate_trace/comparison.json` retain
per-case results, commands and SHA-256 identities. Native output reproduction,
observer invariance, internal RHS comparisons and the candidate vector-weight
control are separate evidence classes; none qualifies MPI/restart or physics.

### 2026-10-31 reaction-resolved attribution and foreign-kernel controls

`build/isothermal-slip-reaction-progress-20261031/` compares forward and
reverse progress rates at all2,667 native internal trial states, supplying
identical molar concentrations and recovered temperatures. Native GET_RMAP
ordering is applied and all21 reaction stoichiometric vectors are checked
against the Fortran mechanism. Reconstructing species production from the
recorded progress rates agrees with each implementation's production kernel
within6.54e-16 in the density-rate norm used here.

For the existing activation-compatible mechanism, the largest forward-rate
relative difference is2.4025199084966834e-12; the reverse maximum is
3.148914885728517e-7. Reaction8 (zero-based), H + O2 (+M) <=> HO2 (+M),
has the largest absolute species-rate contribution difference,
1.466310141339574e-9 kg/m3/s. Its reverse-rate discrepancy is attributable
to the already documented rounded equilibrium reference gas constant:
native8314.46 versus canonical8314.46261815324 in SI kmol units.
Analytically adjusting each reverse rate by (Rnative/Rcanonical)^delta_nu
reduces the species-production norm discrepancy from3.29282169977043e-10
to5.6323490946680257e-14, and the maximum reverse relative discrepancy to
2.435805111845459e-12. This is reaction-resolved attribution on the current
trial states, not a new production constant or an adopted mechanism change.
It does not supersede the earlier EB evidence that reducing absolute rate
error need not resolve spatial limiter sensitivity.

`build/isothermal-slip-native-kernel-control-20261031/` directly separates
reaction arithmetic from the candidate integrator by calling native kernels
only in isolated diagnostic objects. Both controls retain the existing full
Fortran CVODE path, native-derived vector weights, activation-compatible
approximate Jacobian, and Fortran initial/final thermodynamic validation.
The species-only control uses native RTY2WDOT at the candidate recovered
T, then retains candidate thermal-rate assembly. The full-RHS control calls
actual native fKernelSpec, including affine-energy recovery and thermal rate.
The C binding reproduces all numerical outputs of the standalone native RHS
oracle exactly on2,667 saved trial inputs (`binding/receipt.json`). This adds
no production foreign-solver dependency.

Both controls complete21 integrations with2,478 accepted steps. The maximum
H2O2 endpoint differences are2.6211285660274064e-23 kg/m3 for native species
rates and2.626298444855863e-23 for the complete native RHS, compared with
1.494692102428593e-20 for the preceding activation/vector-weight candidate.
Nevertheless, H2O2-only algebraic donor substitution gives z factors
0.8847252557108484 and0.8849416538668086, with transverse ghost-velocity
differences2.4331192049587003e-5 and4.11298662665871e-5 m/s. These are neither
full-state substitutions nor evolved CFD trajectories. The fixed-reference
factor remains0.884411824720706. Matching reaction kernels therefore reduces
absolute endpoint error but does not establish the required spatial ratio or
trajectory agreement. Native and candidate arithmetic units, approximate
Jacobians, initial/final recovery and solver builds remain distinct in these
controls; the residual cannot be uniquely assigned to any one of them.

No numerical production change is adopted, no long CFD is launched, and C2
remains open. Scripts, exact compile/link commands, per-reaction and per-case
results, inputs and binary/source/output hashes are retained in the two
artifact directories. Provenance, project contract, new Python syntax and
`git diff --check` pass; no full configuration or physical validation claim
is made.

### 2026-10-31 full native CPU Jacobian control

`build/isothermal-slip-native-jacobian-control-20261031/` extends the preceding
foreign-kernel controls with the complete native CPU dense Jacobian assembly.
It uses native RTY2JAC, the actual molecular-weight scaling and the native
unwritten T,T entry (zero in a freshly zeroed matrix), followed by the exact
CGS-to-SI row/column transformation. A separate executable invokes the actual
linked native `cvode::cJac` and compares every entry at all2,667 recorded
native trial states. All100 entries match exactly in every state. This verifies
the extracted diagnostic matrix against the native callback, not the production
Fortran approximate Jacobian against a mathematical derivative. An initial
oracle-wrapper compile failure attempted to copy the noncopyable SUNDIALS
context; changing the wrapper to bind a reference fixed it. The failed build
log is retained and is not counted as a passed check.

Both same-input21-case controls complete with2,478 accepted internal steps.
With native species rates and the complete native Jacobian, maximum H2O2
endpoint difference is2.654732778412373e-23 kg/m3 and the H2O2-only algebraic
substitution gives z factor0.8847086784346628, transverse ghost-velocity
difference2.304432220157082e-5 m/s. With native full forced RHS and the complete
native Jacobian, the corresponding values are2.6340532630985475e-23 kg/m3,
0.8845881997665744 and1.3691738367072048e-5 m/s. The preceding full-native-RHS
control with the candidate approximate Jacobian had4.11298662665871e-5 m/s.
The intervention therefore changes the sensitive spatial ratio but does not
eliminate its discrepancy. It does not justify adopting a source change or
claiming evolved trajectory agreement.

The diagnostic still integrates SI density coordinates with the candidate
SUNDIALS build and retains Fortran initial/final thermodynamic validation;
those distinctions remain to be separated. All original production tolerances,
mechanism files, reference limits and numerical sources are unchanged. No
long CFD run or production numerical fix is introduced. Commands, complete
per-case outputs and hashes are in `comparison.json`; the independent actual
native callback verification is in `jacobian-receipt.json`. C2 remains open.

### 2026-10-31 Jacobian refresh-frequency mismatch and implementation

The combined CGS/native-kernel control in
`build/isothermal-slip-cgs-native-core-20261031/` uses raw native inputs,
native initial/final inverse EOS, native forced RHS and native full Jacobian
through the candidate Fortran SUNDIALS API. Before matching the refresh policy,
all21 cases have matching step/RHS counts (2,478/2,667), temperature and energy,
but species differ by at most4.756288522179909e-24 g/cm3, including H2O2 at
2.592512478918333e-26 g/cm3. Relinking the same program against the actual
native SUNDIALS shared libraries, while retaining the Fortran interface
archives, produces byte-identical output to the candidate static SUNDIALS
build. The library paths are verified with `ldd` in `shared-ldd.txt`.

Inspection of the native CPU `ReactorCvode::initCvode` identifies an explicit
`CVodeSetJacEvalFrequency(...,100)` call. The candidate full-coordinate path
had left the CVODE default (`CVLS_MSBJ=51`) in effect. Adding only this
source-derived setting to the CGS control makes all21 species vectors,
temperatures and energies exactly equal to the native CFD reactor observations,
while the total step/RHS counts remain2,478/2,667. Thus aggregate solver counts
had concealed a distinct Jacobian update schedule. `frequency100-comparison.json`
retains this one-setting intervention separately from kernel and unit changes.

The production `advance_forced_reactor_cvode_full` now explicitly sets
`FCVodeSetJacEvalFrequency(memory,100_c_long)` and checks the return status
inside the existing transactional solve. The reduced-coordinate path is
unchanged. This matches the pinned CPU policy without changing RHS/Jacobian
formulas, mechanism data, coordinate units, tolerances or acceptance gates.
The independent CGS/native-kernel control is a diagnostic and is not part of
the production implementation.

The rebuilt canonical production full-coordinate replay completes all21 saved
inputs with2,436 accepted steps. Its maximum H2O2 endpoint discrepancy remains
1.1103327962228815e-15 kg/m3; maximum species change from the previous full
implementation is5.3336605897419224e-21 kg/m3. The exact result from the native
kernel control must not be attributed to production PeleF. This implementation
fixes one proven CPU policy mismatch, not C2 trajectory or convergence.
Three focused tests pass: `unit_sundials_forced_reactor_full`,
`unit_sundials_forced_reactor` and `unit_forced_signed_trials`, including
analytic mass/energy, signed intermediates, admissible endpoints, context reuse
and failure rollback. MPI/restart verification for this increment is recorded
separately below after completion.

The refresh-policy implementation's application verification subsequently
completed successfully. `application/reacting_restart_vo97eecy/checks.json`
under the CGS-core artifact reports15 expected outcomes (7 successful runs,
8 deliberate refusals). An independent receipt rehashes every run input,
both rebuilt executables and the report; serial, MPI2, MPI2-to3 restart and
MPI-to-serial restart final-state CSVs are byte-identical. Changed-reactor
restart identities and the existing malformed-input refusals remain enforced.
This is bounded reacting EB-AMR lifecycle and conservation evidence, separate
from fixed-reference trajectory agreement. `application-receipt.json` records
the verification; no pending application job remains for that artifact.

A separate current-source n8/300ns serial/MPI3 fixed-reference comparison is
recorded under `build/isothermal-slip-full-jac-refresh-20261031/`. It uses the
original TGReact checker and limits with the explicit full-coordinate probe,
reusing the unchanged fixed-native reference. Its status must be read from
that artifact's `record.json`; the completed lifecycle checks above do not
imply a pass for this comparison or the outstanding n64/C2 gates.

The post-fix n8/300ns fixed-reference comparisons now pass original gates in
serial and MPI3. All8 state/budget CSV pairs are byte-identical. Independent
`verified.json` rehashes both check reports, both executables, the checker/probe
inputs and all363 current numerical sources. Relative to the previous full
candidate, initial CSVs are unchanged and evolved outputs differ slightly.
At matched times the maximum absolute momentum discrepancies are
4.6371417816715166e-10 (x),4.4428641321481166e-11 (y) and
3.981735947899967e-11 (z); the previous full candidate had
4.6370118855776354e-10,4.442644863100753e-11 and3.979095525297183e-11.
`n8-field-comparison.json` preserves64 matched field/time records, maxima,
locations, original statuses and input hashes for post-fix full, previous full
and baseline reduced cases. This confirms the bounded regression passes;
it does not indicate a material improvement of the existing n8 error.

The same artifact's `run_n64.py` launches the previously failing n64/50ns
comparison with the post-fix full-coordinate executable. The fixed reference,
prior bounded checker and original thresholds are unchanged. Its independent
`n64-record.json` and `n64_notification/` own that run's pending/terminal status;
no pass is inferred from the successful n8 cases. Once complete,
`compare_field_gates.py --cells 64` can compare all three candidates at the
matched0/25/50ns times without recomputing CFD or changing gates.

### 2026-11-01 post-refresh n64 result: improved, still not qualified

The post-fix n64/50ns full-coordinate comparison completed in23,075.889 s
(about6h25m). The wrapper returned0 because the comparison completed; the
original checker returned1 and correctly reports failure. `n64-verified.json`
rehashes the checker, executable, report and all363 source files and confirms
all three0/25/50ns output times and the unchanged original thresholds.

At50ns, only x-momentum fails: maximum absolute error1.3056688216970424e-6
versus limit1e-6, at fine(level1,i64,j121,k32). This is19.31% below the
previous full candidate's1.618157194350456e-6, but slightly above the baseline
reduced candidate's1.2931953389840045e-6. The previous full candidate also
failed total/internal energy; those fields now pass, with scaled maxima
9.33208650326718e-9 and9.332345812877956e-9 against1e-8. Maximum y/z momentum
errors are2.498293370095661e-7 and1.986263631103863e-7, also below the old full
candidate but above the reduced baseline. Gas/energy ledger residuals remain
at most4.002664524582677e-17. The improved pressure/energy agreement does not
close the remaining momentum gate or establish spatial convergence.

`n64-field-comparison.json` preserves matched field/time maxima and locations
for the three candidates. The maximum location moved from old full
(i64,j8,k97) and reduced(i64,j105,k121) to the point above; no causal claim
about that new point is made from saved endpoints. A separate saved-CSV audit
of all16,384 finei64 cells is tracked in `interface_notification/` and will
quantify error distribution and exceedance counts without another CFD run.

The saved-plane audit completed successfully (`postfix-interface-audit.json`).
For all16,384 finei64 cells at50ns, post-fix x-momentum mean absolute error is
1.281680887297433e-8, median3.6752061376849054e-10 and95th percentile
6.680426991711313e-8. The old full values were1.3148581962149008e-8,
3.695057237615429e-10 and6.749291702693588e-8. Thus both aggregate errors and
the maximum improve, but12 cells still exceed1e-6, unchanged from old full
and worse than the reduced baseline's8 cells. The reduced baseline also has
lower mean and95th-percentile errors (1.2287787635372785e-8 and
4.0877289419691265e-8), despite its larger median.

The12 post-fix exceedances form six adjacent-k pairs at (j,k pairs)
(121,31:32),(8,97:98),(57,33:34),(57,95:96),(8,31:32),(72,33:34), all at
finei64. The current maximum is not the location of the prior captured first
fine-predictor discrepancy. The earlier local H2O2/limiter causal evidence
therefore cannot be relabelled as a direct observation of this new maximum.
Further stage-local donor/pressure/face-flux capture should cover the new
cluster before attributing the remaining trajectory error. Endpoint screening
alone does not prove that the same component controls the in-step limiter.
No additional numerical change is made and no further CFD job is launched
as part of this plane audit. The refresh-policy correction remains a verified
implementation correction; switching the C2 trajectory to the full-coordinate
path remains unqualified.

### 2026-11-01 paired observations for the new maximum (pending)

Observer-only runs in `build/isothermal-slip-newmax-stages64-20261101/` and
`build/isothermal-slip-native-newmax64-20261101/` target fine(64,121,32).
Candidate root cells32:34,59:63,14:18 and fine cells63:64,119:123,29:35 record
stage states, reconstructed pressures, all directional face fluxes and root
synchronization phases. The native observation region includes the matching
faces and central reconstruction stencils. The candidate links the post-fix
full-coordinate probe and runtime; neither production source nor numerical
settings are changed. Six candidate and three native output CSVs must match
their respective completed baselines exactly before stage comparisons are
interpreted. Prepared analysis scripts pass Python syntax checks; actual stage
analysis and observation invariance remain pending in these run records.

The native new-maximum observer has completed. Independent
`native-capture-receipt.json` in its artifact verifies the three exported CSVs
against the unchanged reference, source/executable identities and capture
coverage. There are17,712 unique state records with41,760 exact duplicates,
and47,088 unique face records with63,900 exact duplicates. All requested
7,740 face observations,1,290 selected cell-state observations and432 central
reconstruction-stencil observations are present across six root and twelve
fine hydro calls. The candidate observer is still running; no paired pressure,
flux or causal conclusion is inferred from the native-only integrity result.

### 2026-11-01 new-maximum paired stage budgets and donor reconstruction

The candidate observation has now completed. Its six state/budget CSVs are
byte-identical to the post-Jacobian-refresh full-coordinate n64/50ns run;
the native observation's three exports match the unchanged pinned reference.
Artifacts are `build/isothermal-slip-newmax-stages64-20261101/` and
`build/isothermal-slip-native-newmax64-20261101/`. Re-running `analyze.py`
and `momentum_attribution.py` verifies1,290 candidate momentum budgets,
450 synchronization component records and paired momentum closure at
fine cells(64,121,31:32). Candidate update closure is at most1.152e-17;
paired candidate-minus-native closure is at most2.105e-15.

At fine(64,121,32), the x-momentum difference is only -4.442e-11 at25ns,
then reaches5.798e-7 at37.5ns and1.305669e-6 at50ns. The last two fine
intervals contribute respectively3.827973e-7 and5.232540e-7 through the
transport-rate difference, and1.970190e-7 and2.026523e-7 through x-direction
hydrodynamic face fluxes. The y/z hydrodynamic contributions are each below
6.1e-12 in magnitude. This is an update decomposition, not proof of a unique
transport or pressure defect: the candidate observation contains aggregate
transport rates, while native directional transport face fluxes are recorded.
The adjacent k=31 cell develops the opposite-sign error.

`reconstruct_candidate.py` blends the seven root donor states at interval2
stages0 and3 with the source's exact half-time weighting, then replays the
shared limiter. Reconstructed ghost(65,121,31:32) velocities match the actual
event14 observation within2.776e-17 m/s. At donor(33,61,16), the candidate
shared factors are0 in x,0.3605697617880513 in y (controlled by temperature),
and0.6575792553414581 in z (controlled by z momentum). This is distinct from
the earlier maximum's H2O2-controlled stencil. No native in-step donor match
or causal correction is claimed yet.

A separate native-only observer in
`build/isothermal-slip-native-newmax-donors-20261101/` is collecting those
actual interpolation inputs. Its completion must verify all three reference
exports before `analyze.py` accepts the paired donor comparison. No new
production numerical changes or acceptance-threshold changes were made.
Python syntax, mechanism provenance, project contract and `git diff --check`
pass for this diagnostic increment. C2 trajectory and convergence remain open.

### 2026-11-01 actual new-maximum donor comparison

The native donor observer completed in29min22s including build/export. Its
three reference exports are byte-identical and its source hashes are unchanged;
`build/isothermal-slip-native-newmax-donors-20261101/analyze.py` verifies these
conditions before comparing actual preceding-event13 interpolation inputs.
Native CGS replay reconstructs event14 ghost(65,121,31:32) velocities within
1e-12 m/s. The independently pinned AMReX slope executable matches all96
same-input slopes exactly in `kernel-check.json`, using each solver's original
units. This does not establish identical evolved interpolation inputs.

Both solvers select temperature for the y-direction shared limiter and
z momentum for z. Native factors are0.3599453517855147 and
0.6516726950702059, compared with candidate0.3605697617880513 and
0.6575792553414581. In `factor_attribution.json`, substituting only the native
z factor into the candidate stencil reduces ghost(65,121,31) x-velocity
difference from -6.092137e-5 to9.441900e-8 m/s, and y-velocity difference
from -4.934457e-4 to -7.059794e-6 m/s. Replacing both y/z factors leaves
x/y differences -1.699431e-8 /3.378792e-11 m/s. These are algebraic
substitutions, not a CFD correction or a qualification result.

`root_budget.py` traces the donor z-momentum difference through observed root
face fluxes and transport rates, using actual full-time interpolation inputs
before fine synchronization as the endpoint. At root(33,61,16), interval2
starts at3.198155e-12 and ends at3.565103e-10 kg/(m2 s). The z-hydrodynamic
increment difference is2.680998e-10 and the transport increment difference
8.522074e-11; x/y hydro increments are below7.2e-15. Across three donor
cells and two root intervals, paired momentum closure is at most1.371e-15.
Thus the z-limiter factor difference includes prior evolved pressure/transport
feedback; changing the limiter alone is not justified by this evidence.

A backward trace identifies an earlier nearby ghost(65,122,33) at event5.
Replaying its candidate root donor(33,61,17) half-time stencil reproduces
all three observed velocities exactly. Its y factor0.7706878639251975 is
controlled by HO2, while x/z factors are zero. The current native donor
capture does not include its complete seven-cell stencil, so another isolated
native observer, `build/isothermal-slip-native-newmax-seed-20261101/`, is
collecting that stencil. Its analysis requires exact reference exports and
same-input slope verification. This earlier seed is not yet assigned a
chemical or interpolation cause. No production equations, thresholds or
default solver selection changed in these diagnostics. Python syntax checks
pass; C2 remains incomplete.

### 2026-11-01 first-interval seed identified in the native donor stencil

The additional native seed capture completed with exit0 and two-step exports
(all three files) exactly matching the original reference. In
`build/isothermal-slip-native-newmax-seed-20261101/`, `analyze.py` verifies
those hashes and reconstructs ghost(65,122,33) at event5 from the actual
preceding-event4 CGS stencil. Velocity residuals are0,0,5.294e-23 m/s;
all96 same-input AMReX slopes match exactly.

Both solvers' y limiter is controlled by HO2 at root donor(33,61,17).
Native factor is0.7665521678057533 versus candidate0.7706878639251975.
Native half-time HO2 densities along y-minus/center/y-plus are
5.249572569496193e-8,5.2495725695404485e-8,5.2495725695509404e-8 kg/m3.
Candidate-minus-native differences are approximately -6.273503e-14 kg/m3,
but the variation of that difference across the stencil is only of order
1e-22 kg/m3. The near-flat species profile makes its ratio-based shared
limiter sensitive to these tiny differential discrepancies.

`factor_attribution.py` substitutes the native dimensionless factors into
unchanged candidate SI donor values: ghost x/y velocity differences fall
from8.219749e-7/-4.676276e-5 to -1.379903e-13/2.181588e-14 m/s.
This identifies the local algebraic amplification, not its chemical origin
or an evolved remedy. The interpolation implementation is not changed.

A first-root-interval reaction observer is now collecting native pre/post
states at these seven donors in
`build/isothermal-slip-newmax-seed-reactor-20261101/`. Prepared controls
replay21 actual native inputs through production full-coordinate chemistry,
actual native ReactorCvode and the previously verified native-CGS Fortran
control. Two additional isolated variants use the already generated pinned
activation-basis source, with/without the native LiDryer equilibrium reference
concentration constant. The latter changes only a local copy of the generic
kinetics module, verifies the pinned native header SHA and does not change
canonical mechanism data or production source. Both variants compile/link;
no new replay outcomes or CFD improvement are claimed before data arrive.

The seed-reactor diagnostic also compiles two native-weight controls
(species absolute tolerance1e-15 kg/m3 and temperature1e-18 K, with/without
the equilibrium-constant change), preserving the100-step Jacobian policy.
These are isolated module copies, not production tolerance changes. A smoke
replay of the older donor's21 inputs completes for the activation and
activation-plus-equilibrium controls (2,436 steps each); HO2 maximum endpoint
differences are1.642267e-17 and1.904064e-18 kg/m3, respectively. These older
inputs do not establish the new seed's behavior or a corrected CFD trajectory.
The new21-input comparison remains pending capture completion. Python syntax,
mechanism provenance, project contract and diff-whitespace checks pass.

### 2026-11-01 new-seed identical-input reactor controls and bounded CFD trial

The seven-donor first25ns reactor capture completed with exit0; both exports
match the original reference. `replay_analysis.py` verifies capture/source
hashes and replays21 actual inputs. Actual native ReactorCvode and the
native-CGS Fortran diagnostic with Jacobian refresh100 reproduce all21
observed species/temperature states exactly. The full-coordinate production
path and source-derived controls have these maximum HO2 endpoint differences:

| Control | max HO2 difference, kg/m3 |
| --- | ---: |
| Production full-coordinate canonical mechanism | 1.258338337686366e-13 |
| Native activation basis | 1.642267182806040e-17 |
| Activation plus native equilibrium reference concentration | 1.905175621721934e-18 |
| Activation plus native component absolute weights | 1.830235705202590e-17 |
| Activation, equilibrium constant and native weights | 3.269017780809565e-21 |

These are identical native input replays, not evolved CFD agreement. In the
HO2-only algebraic substitution, the combined control leaves y-velocity
difference -7.400842e-6 m/s at ghost(65,122,33). Even substituting exact
native endpoint HO2 after CGS-to-SI conversion leaves -2.195017e-7 m/s:
unit conversion can perturb the ratio of nearly cancelling differences.
Reduced endpoint norms therefore do not prove matching limiter factors or
convergent trajectories.

An isolated full-coordinate CFD trial is running in
`build/isothermal-slip-native-compat-trial64-20261101/`. It links the
source-derived activation, equilibrium and native-weight controls while
keeping the fixed reference and all original field/budget limits. Its campaign
requires n8 serial and MPI3 reference gates and eight byte-identical MPI CSV
pairs before launching n64/50ns. No outcome is claimed until completion.
The trial is diagnostic and is not adopted into the production applications.

Independent invariant checks found a limitation of the uncapped weight
control: it replaces a caller's stricter1e-19 species tolerance with1e-15,
and fails the existing signed-intermediate analytic endpoint test. Reduced
reactor and signed-trial kinetics tests pass (two passes, one failure overall).
An initial test-link orchestration error, assuming that the kinetics-only
signed test links CVODE, is retained separately and not counted as a test.
Neither failure is hidden by changing an expected value or acceptance limit.

A separate `capped_weights/` control uses the smaller of the caller's absolute
tolerance and each native component tolerance. All three existing focused
tests then pass, including the full-coordinate signed endpoint, mass/energy,
rollback and context-reuse checks. Its21 new-donor replay outputs are byte-exact
to the uncapped control for both mechanism variants. This establishes only
those input cases; the running CFD trial retains its original uncapped binary
and source hashes and must not be relabeled as a capped-policy trajectory.
No production source/default/tolerance changed. Python syntax, provenance,
project contract and diff checks are separate from C2, which remains open.

### 2026-11-01 native-compatibility CFD result and rejected arithmetic shortcuts

The diagnostic campaign completed. n8 serial/MPI3 passes the original field
and budget gates with eight byte-identical state/budget CSV pairs. The
n64/50ns calculation takes23,269.5s and still fails the unchanged x-momentum
limit:1.0048058732793283e-6 versus1e-6, at fine(64,26,109). All other original
field gates pass. Maximum y/z differences are2.0724378857472914e-7 and
6.023157189225957e-8. Independent `verify.py` checks the363-source manifest,
linked diagnostic source/object/binary identities, checker records and MPI
pairs; gas/energy budget residual maxima are1.3262e-16/4.4110e-17. A receipt
script initially treated the species-budget dictionary as a scalar; its
aggregation was corrected without modifying solver results or gates.

At the entire16,384-cell fine i=64 plane, x-momentum mean error is
1.2542036996e-8,95th percentile5.2058287354e-8 and maximum1.0048058733e-6.
Three cells exceed1e-6, compared with12 for the production full-coordinate
refresh100 baseline. The maximum is about23% smaller, but remains a failure;
the maximum has moved. This uncapped-weight diagnostic remains unsuitable
for adoption because of the previously recorded tighter-caller-tolerance
failure. No capped-policy CFD trajectory is inferred from it.

Additional bounded controls under
`build/isothermal-slip-newmax-seed-reactor-20261101/` do not justify another
large trajectory on their own:

- `cgs_coordinates/` implements a covariant density-coordinate change inside
  CVODE: SI API/thermochemistry retained, CGS ODE densities, transformed RHS,
  Jacobian, tolerances and endpoint admission. Three focused tests pass.
  On21 actual inputs, the activation/equilibrium control has HO2 maximum
  difference2.95138e-21, but its HO2-only ghost y-velocity difference is
  +8.75903e-6, versus -7.40084e-6 for the capped SI-coordinate control.
- `direct_arrhenius/` uses native-style prefactor multiplication and inverse-T
  evaluation only within the safe exponential range, preserving the logarithmic
  fallback. On42 identical pre/post states, production/Jacobian norm differences
  remain approximately5.16e-14/1.98e-13, essentially unchanged from5.10e-14/
  1.99e-13. The21-input HO2 maximum falls slightly to2.73962e-21, but local
  ghost y-velocity difference increases to1.13185e-5. It is not adopted.
- `native_jacobian_policy/` tests an independent Fortran approximation to the
  native preconditioner policy: chemical-only thermal density row and zero
  temperature-temperature entry, with physical forced RHS unchanged. The
  existing finite-difference species temperature column remains. An initial
  missing empty-reaction branch fails the inert test; that diagnostic-only
  defect and receipt are retained. With an explicit zero-production branch,
  all three focused tests pass. The new donor's y-velocity substitution improves
  to3.13383e-6, but the earlier H2O2-controlled donor worsens from1.21527e-4
  to1.38482e-4. This is not consistent evidence of a trajectory correction.

These comparisons reinforce the distinction between tiny species norm errors
and ratios of nearly cancelling neighboring values. All controls remain
isolated diagnostics; no production equation, mechanism, default selection,
acceptance limit or golden hash is changed. C2 remains open, and neither the
50ns trial nor the focused tests qualify C3--C5 or physical validation.

### 2026-11-01 remaining-tail paired observation launched

The native-compatibility trial's three fine-plane x-momentum exceedances are
(64,26,109),1.0048058732793283e-6;
(64,26,110),1.0041781991006367e-6; and
(64,103,109),1.0000977055855086e-6 kg/(m2 s).
A paired observation now targets both y neighborhoods, without changing the
trial's numerical settings. Candidate artifacts are in
`build/isothermal-slip-compat-tail-stages64-20261101/`; native artifacts are
in `build/isothermal-slip-native-compat-tail64-20261101/`.

The candidate records72 root and72 fine cells, pressures, directional hydro
face fluxes, primitive reconstruction stencils, transport rates and coarse/fine
synchronization stages. Its driver requires six CSVs to match the completed
native-compatibility trial exactly. This is the same uncapped diagnostic,
not a capped-policy qualification run. The native observer records matching
states/faces and an expanded neighborhood of actual seven-cell interpolation
stencils; donor rows now have keys(event,level,i1,j1,k1) followed by112 rawCGS
values. Three reference exports must match before native analysis is accepted.
Both native processes have been verified live. Analysis scripts are prepared,
but no new stage result is claimed while captures are running.

The native remaining-tail observer completed with exit0. Independent
`verify_native.py` confirms all three reference CSVs, source/binary hashes,
38,520 unique state records (115,560 exact tile duplicates),101,520 unique
face records (192,600 exact duplicates), and1,344 unique actual donor stencils.
`analyze_native_donors.py` reconstructs432 native ghost velocities across
12 fine hydro events and36 cells. Maximum reconstruction residual is
8.881784197001252e-16 m/s. The first fine call of each root interval uses
preceding-event keys(3,level0)/(12,level0), not level1; an initial analysis key
lookup was corrected to reflect the observed call ordering.

For native ghosts(65,26,109) and(65,103,109), shared factors are[0,0,0]
at event5, with z momentum among the controlling components; they are
[0,1,1] at event14. These native-only observations do not establish the
candidate's factors or explain the remaining errors. Candidate PID4031259
is verified live on the original capture. Its output identity and paired
stage comparison remain pending; only that remaining job is returned to
the completion watcher.

### 2026-11-02 remaining-tail stage attribution and admitted-composition change

The candidate tail observation completed in22,270.2s with all six CSVs exactly
matching the unchanged native-compatibility trial. Its1,296 momentum budgets
and432 synchronization component records close within5.543e-17; all7,776
paired face records are covered. At the three failed cells the paired momentum
update closes within2.153e-15. At(64,26,109), x-momentum difference rises from
8.119e-11 at25ns to1.265175e-6 at37.5ns; the transport increment difference is
1.175244e-6, versus9.044930e-8 from x-hydrodynamic faces. It decreases to the
previously reported1.004806e-6 at50ns. These measurements are in
`build/isothermal-slip-native-compat-tail64-20261101/`.

`paired_donors.py` reconstructs192 paired ghost states with maximum velocity
residual2.665e-15 m/s. At ghost(65,26,109), event13 has candidate/native z
factors0.4634583794073112/0.4484284557144933, both controlled by z momentum.
The velocity difference is(+4.44096e-5,-2.78169e-3,-1.62816e-10) m/s.
Replacing only the shared factors algebraically leaves x/y differences
-1.19141e-9/-3.90799e-13. At event14 both factors are[0,1,1], with velocity
differences below2e-9. This is a start-of-interval amplification, not a persistent
mismatch at every predictor/corrector stage.

The first25ns synchronization changes the root donor(33,13,55) z-momentum
difference from2.94417e-14 to2.96949e-11. Candidate/native hydro reflux
contributes only -1.19475e-19 of that difference; the inferred transport reflux
contribution is2.96655e-11. Native reflux reconstructed from the measured coarse
face and four fine faces agrees with the observed donor update. For this donor
and its y counterpart(33,52,55), using the observed velocity gradients and
native viscosity inferred from the off-diagonal stress reconstructs the
transport reflux difference within2.632e-19. The dominant difference is the
transverse du/dz term in the x-face z-momentum viscous flux. Candidate transport
coefficients were held at native values in this algebraic reconstruction, so
this is not a general same-input transport-kernel qualification.

An earlier nearby ghost(65,26,108) already has z-factor0.03157503644517107
versus native0.03232456921788181 at event5, again controlled by z momentum.
Before any fine-grid evolution, the precursor root(33,13,54) z-momentum
endpoint differs by -1.22111e-14; its update discrepancy is predominantly the
z-direction hydrodynamic pressure flux. Pressure-face discrepancies at the
last root corrector are of order1e-10 Pa. Thus the measured chain is a small
root pressure/velocity discrepancy, shared-limiter velocity amplification,
viscous interface reflux, and renewed limiter amplification. It does not
establish a defect in the synchronization weights themselves.

`build/isothermal-slip-tail-primitive-20261102/` compiles the actual native
`pc_ctoprim` and the production Fortran primitive converter. Both exactly
reproduce their432 observed root states. Cross-evaluating both sets,864 inputs
per solver with explicit SI/CGS conversion, gives maximum pressure difference
3.637978807091713e-10 Pa and maximum native-input temperature difference
1.1368683772161603e-12 K. An initial reporting index selected QREINT instead
of QTEMP; it was corrected using IndexDefines.H before interpretation.

A concrete source difference is the additional normalization of admitted
rhoY/rho in Fortran; native pc_ctoprim preserves those ratios. An isolated
control removes only that second normalization, retaining positivity, closure
and finite checks. Across both input sets its maximum pressure difference
falls to7.275957614183426e-11 Pa; mean differences fall from3.702e-11/4.261e-11
to1.499e-11/1.452e-11. Cross-solver comparisons still include unit-conversion
roundoff and do not prove an evolved trajectory correction.

A focused production change now adds an optional `preserve_species_ratios`
primitive choice and enables it for nonrotating regular MOL through the shared
MOL hydro adapter. Defaults elsewhere retain the prior normalization. Conserved
states, species-admission tolerances, mechanisms and reference gates are
unchanged. Regression coverage tests the retained ratios, unchanged default,
closure refusal, immutable input, and uniform MOL partial-density face flux.
Preimages and the isolated four-file patch are recorded under
`build/regular-mol-admitted-composition-20261102/`. The focused build completed and
`unit_mol_riemann_3d` plus `unit_reactive_mol_3d` passed (2/2,1.46 s).
Project contract, diff whitespace and six campaign Python syntax checks passed.
Byte provenance passed for h2o2_full and lidryer; applying the same checker
to h2o2_elementary was rejected because that bundle has no source object,
so it is not counted as a passed provenance check. Application/restart
checks and numerical qualification remain pending. The full and reduced
canonical-reactor campaigns under `build/regular-mol-composition-{full,reduced}-20261102/`
first require original n8 serial/MPI3 gates and eight exact CSV pairs before
n64/50ns. This is an implementation, not a passed C2 result.

The admitted-ratio implementation's rebuilt fixed and selected applications
subsequently pass68 thermal-slip outcomes each (136 total), covering common
and recursive SDC/MOL, exact serial/MPI output, rank-changing restart and
intentional refusals. `thermal-verification.json` and
`verification-receipt.json` in the implementation artifact retain command,
binary and report hashes; all136 per-run logs are present.

The first three generic CTU restart regressions failed in their checkers:
serial checks assumed three total steps although both uninterrupted and
restarted runs now take four, and the corruption helper treated the optional
timestep-history tag as the schema/context length. The separate checker fix
preserves exact CSV comparison and now requires restarted step count equal to
the uninterrupted run. It locates the byte-counted context and ledger through
the optional headers, so finite-ledger and index corruption still target the
intended records. No solver tolerance or checkpoint reader was changed.
Initial failures remain in `restart-tests.log`; after the checker correction,
three serial/MPI subcycle regressions pass in102.95s, and the indexed MPI
restart regression passes in29.33s. The checker changes and preimages are
recorded separately in `checker-fix.json`, `checker-fix.patch`, and
`checker-before/`. Python syntax, diff whitespace and project contract pass.
These lifecycle results do not establish fixed-PeleC trajectory convergence.

The canonical reduced-coordinate admitted-ratio campaign completes n64/50ns
in12425.82s. Original n8 serial/MPI3 field gates pass with eight byte-identical
CSV pairs. The n64 comparison still fails one original gate: x-momentum at50ns
is1.309126998955308e-6 against1e-6, at fine(64,105,121). The old reduced run
was1.2931953389840045e-6 at the same cell: this intervention does not improve
the evolved maximum despite improving same-input pressure reconstruction.
At25ns x-momentum maximum is6.784679513504623e-7 at fine(64,41,104).
Maximum gas/energy budget residuals remain approximately1.16e-16. The
`verification.json` in `build/regular-mol-composition-reduced-20261102/`
rechecks363 source hashes, binary/report hashes and both small-grid receipts.
The campaign driver's successful exit means the bounded comparison completed;
the reference comparison itself failed. The full-coordinate counterpart is
still running. Neither trajectory nor C2 convergence is qualified.
The whole fine i=64 interface plane (16384 cells) confirms eight x-momentum
exceedances, unchanged in count. Mean absolute x-momentum difference grows
from1.2287787635372785e-8 to1.2829232506357674e-8; y/z momentum means also
increase slightly. `interface-audit.json` binds candidate and unchanged native
CSV hashes and retains prior reduced/full summaries separately. This rejects
an interpretation that the unchanged maximum merely concealed a general
plane-wide error reduction from this intervention.

The admitted-ratio full-coordinate campaign also completes, in23763.95s,
with original n8 serial/MPI3 gates and eight exact CSV pairs. At n64/50ns
it fails four gates: density scaled1.0312724698608185e-8 (limit1e-8),
x-momentum absolute1.983209212836594e-6 (limit1e-6), total-energy
scaled1.4142906191232826e-8 and internal-energy scaled1.4143300507729766e-8
(both limits1e-8). The x-momentum maximum moves to fine(64,72,33);
density/energy maxima occur at(64,71,34). The interface-plane audit gives
16 x-momentum exceedances and mean1.3308067288963277e-8.
`build/regular-mol-composition-full-20261102/verification.json` independently
rechecks the363-source manifest, binary and report hashes. Thus the primitive
policy correction is not an established trajectory improvement in either
integrator. No tolerances or reference hashes were changed.

A new observational pair is running under
`build/regular-mol-composition-stages64-20261102/` and
`build/regular-mol-composition-native64-20261102/`. It targets root
(32:34,35:37,15:18) and fine(63:64,70:74,29:36), including the new maximum
and its interpolation donors. Candidate observation retains current canonical
mechanism/full-reactor source and adds direct direction-specific transport
face flux recording, alongside hydro pressure/flux, stage momentum and sync
observations. Native observation records actual donor stencils and hydro/
transport faces. Acceptance requires candidate six-CSV and native three-CSV
byte identity with their already completed trajectories before causal analysis.
These captures are pending and do not constitute additional passed checks.

The fixed-native new-maximum observer completes and reproduces all three
reference CSVs exactly. Independent source/binary/output rehash and duplicate
checks are in `regular-mol-composition-native64-20261102/native-receipt.json`:
23952 unique state records,63288 unique face records,972 unique donor records;
all tile duplicates match exactly. An initial copied verifier still selected
the previous region and rejected missing old coordinates; the retained
`verify_current_native.py` checks the actual(33,36,16:17) donor region.
Native-only interpolation replay reconstructs480 ghost velocities with maximum
residual1.1102230246251565e-16 m/s. At event14 the ghost(65,72,33) has y/z
shared factors0.3596660526167657/0.6499767313648196, controlled respectively
by temperature/z momentum. These are native observations only; a candidate
cause attribution awaits its exact-output capture. No extrapolation from the
old-region candidate data is used.

The candidate new-maximum capture completes with six exact original CSVs.
The retained analysis was adjusted for this region's36 root/80 fine cells
and36 synchronized root cells (the copied old-region counts initially failed).
Its1176 component update budgets close within6.6881721515199555e-18;
hydro face divergence agrees within1.4551915228366852e-11, with216 sync rows.
Paired native/candidate7056 hydro faces show maximum pressure reconstruction
difference0.003409425451536663 Pa. At fine(64,72,33), x-momentum differences
are -6.782689061552057e-11 at25ns,8.809584738007258e-7 at37.5ns,
and1.983209212836594e-6 at50ns. The latter two fine intervals respectively
receive x-direction hydro differences2.9860954149626197e-7/3.0619031400419773e-7
and total transport differences5.824209959081246e-7/7.960703160334313e-7.
Paired momentum closure over the three target cells is2.0429006719171024e-15.

Actual native and reconstructed candidate donor replay matches96 observed
ghosts within2.7755575615628914e-17 m/s. At event14 ghost(65,72,33), the
candidate/native y factors are0.3607909244021695/0.3596660526167657
(temperature-controlled), and z factors0.6589305611017544/0.6499767313648196
(z-momentum-controlled). Velocity differences are9.272145912937041e-5 and
-7.245984433702224e-4 m/s in x/y; substituting native factors leaves
2.6326738520113246e-8/-2.4934776465812547e-12. This algebraic diagnostic
identifies shared-limiter amplification, not a production limiter replacement.
At root(33,36,17), first synchronization changes the z-momentum discrepancy
from1.5123246345940446e-14 to-4.897589393982327e-12; the measured difference
in reflux increments is-4.91271264040106e-12, with zero candidate recovery
increment. Sources are paired-donors, sync-attribution and momentum_attribution
JSONs in `regular-mol-composition-native64-20261102/`.

The newly added candidate transport face observer had a one-face indexing
error: AMR transport arrays use1-based normal face indices, unlike the MOL
hydro diagnostic arrays. The source solver and six trajectory CSVs are
unaffected. `transport-divergence.json` retains the failed direct reconstruction
(max residual11006.49); it must not be treated as valid face attribution.
Reindexing only available faces to their actual positions reconstructs1224
interior component rates within9.094947017729282e-13, confirming the diagnosis;
768 cell-stages lack required faces and are excluded, not inferred. A corrected
observation under `regular-mol-composition-flux64-20261102/` is running to
record the missing interface faces directly. The original bad observation is
preserved, and no full-field rerun result is counted before completion.
`root-budget.json` further separates coarse evolution from fine synchronization.
At root(33,36,17), the first pre-sync z-momentum difference1.5123246345940446e-14
is accounted for mainly by z hydro1.5133991837501525e-14 (closure4.03e-17).
The second coarse interval starts at-4.897589393982327e-12 and ends before
fine synchronization at-5.454557781603463e-10, with z hydro-4.104373510926962e-10
and transport-1.3013324747568333e-10 (closure7.61e-17). Thus the second-interval
shared z-factor discrepancy is already present in the provisional coarse
endpoint; it cannot be assigned solely to that interval's later reflux.
The root predictor/corrector temperature differences evolve from about7.23e-9 K
at the first corrector to1.45e-8/1.37e-8 K in the second interval. These measured
state differences identify the next thermochemical/pressure input comparison;
they do not by themselves identify a reaction-rate or EOS implementation bug.

The corrected transport-face observation completes and again reproduces all
six candidate CSVs exactly. Its3528 momentum transport component rates agree
with direct face divergence within2.7284841053187847e-12;1176 stage updates
retain6.6881721515199555e-18 closure. Direct candidate/native directional
transport fluxes now reproduce the previously inferred transport-increment
differences within4.3807485240808346e-21 over the three target cells.
At fine(64,72,33), x-normal viscous-flux contributions to x-momentum difference
are5.671284483465766e-7 during25--37.5ns and7.776259150794261e-7 during
37.5--50ns; y-normal contributions are1.7245099359525644e-8/2.3639029532083883e-8,
and z-normal contributions-1.9525517979790406e-9/-5.194628578079324e-9.
Thus the later transport difference is primarily across the coarse/fine
normal direction. `transport-budget.json` and `transport-comparison.json`
in the native observation artifact bind the corrected flux capture.

A native first-root-step reactor-input observation is now running under
`build/regular-mol-composition-reactor-20261102/`, targeting root(33,36,17)
and its six immediate neighbors. Current canonical full-reactor replay was
rebuilt and hashed in `replay-build.json`; the prepared replay requires exact
native reproduction of observed endpoints before interpreting same-input
species/temperature differences. This diagnostic separates thermochemical
input differences from evolved pressure/transport feedback and is pending.

The native reactor-input observation completes with two exact reference CSVs.
The actual native CPU reactor replays all21 observed endpoints exactly.
Current canonical full-coordinate replay on those identical inputs gives
maximum species-density difference2.398775622580729e-13 kg/m3 and temperature
difference7.256403478095308e-9 K. At root(33,36,17), event3 temperature
difference7.229346010717563e-9 K equals the recorded first coarse-corrector
temperature difference. This is same-input reactor evidence, separate from
full-trajectory proof. Artifacts: `regular-mol-composition-reactor-20261102/`
`replay-comparison.json`, actual pre/post states and rebuilt replay provenance.

Existing hash-verified source-derived controls reduce maximum temperature
error to4.615685611497611e-11 K with native activation conversion alone and
4.524736141320318e-11 K with the native equilibrium concentration as well.
The activation-plus-native-component-weight control gives1.3642420526593924e-12 K;
adding equilibrium concentration gives1.1368683772161603e-12 K. Those weight
controls already include activation conversion: they are not weights-only
experiments. Caller-capped weights reproduce both uncapped controls byte for
byte for these21 inputs; uncapped weights remain inadmissible as a general
policy because they override stricter caller tolerances. `controls.json` and
`capped-controls.json` retain this distinction. Activation conversion is a
pinned-native mechanism convention, not evidence that canonical SI conversion
is physically incorrect. No canonical mechanism or reference gate was edited.

Two current-MOL trajectory diagnostics are running under
`regular-mol-composition-activation64-20261102/` and
`regular-mol-composition-compat64-20261102/`: activation-only, and activation
plus native equilibrium concentration plus caller-capped component weights.
Each requires original n8 serial/MPI3 gates and eight exact CSV pairs before
n64/50ns. They test whether the measured local thermochemical improvement
survives coupled evolution. They remain explicit diagnostic mechanism/policy
variants and cannot be reported as canonical C2 qualification.

The activation-only current-MOL diagnostic finishes n64/50ns in24589.80s.
Original n8 serial/MPI3 gates pass, with all eight CSV pairs independently
rehashed exact. The363-source manifest, diagnostic source/object receipts,
probe and checker report hashes are verified in
`regular-mol-composition-activation64-20261102/verification.json`.
The original n64 gates still fail for x-momentum1.763584799381944e-6
(limit1e-6), total-energy scaled1.2590162292185247e-8 and internal-energy
scaled1.2590513198333983e-8 (limits1e-8). Density now passes. The maximum
remains fine(64,72,33); x-momentum is about11.1% lower than the canonical
ratio-preserving full run's1.983209212836594e-6. This is partial evolved
improvement, not elimination of the shared-limiter amplification or canonical
qualification. The caller-capped-weight/equilibrium diagnostic remains live.

The activation-only plane audit lowers x-momentum mean absolute difference
from1.3308067288963277e-8 to1.1748508687543497e-8 and exceedance count16 to12;
its maximum remains above the original gate. The combined caller-capped
native-compatibility diagnostic then completes in24637.93s. Its n8 serial/MPI3
gates pass with eight exact CSV pairs, and363-source plus diagnostic-object
and binary/report receipts are reverified. At n64/50ns it fails only x-momentum,
1.1002593854306415e-6 versus1e-6, now at fine(64,103,109). The whole i=64
plane has four exceedances, x-momentum mean1.2546609053621312e-8,
y/z maxima2.2690314160911385e-7/7.700942384871118e-8. Density and both
energy gates pass. These different mechanism/policy diagnostics do not qualify
the canonical source or C2; neither a smaller failure nor improved local replay
is counted as success. `regular-mol-composition-compat64-20261102/verification.json`
and `interface-audit.json` retain exact failures and input identities.

Before another trajectory intervention, isolated primitive arithmetic controls
were tested on the archived864 same-state inputs, whose inputs/native outputs
and provenance hashes were reverified. Reciprocal-density velocity evaluation,
velocity-form kinetic energy and reciprocal-density species evaluation do not
reduce the maximum pressure difference7.275957614183426e-11 Pa. Native-style
pressure multiplication in SI modestly reduces means to1.290135076959376e-11/
1.408032538300311e-11 Pa (native/candidate inputs), and a CGS arithmetic
control gives8.421247238638225e-12/7.848602426410825e-12, but neither lowers
the maximum. These are archived in `regular-mol-primitive-arithmetic-20261103/`
and are not adopted production changes or evolved-trajectory evidence.

The remaining combined-diagnostic exceedances are fine(64,103,109:110) and
(64,26,109:110). A current-policy observation is running in
`regular-mol-compat-tail64-20261103/`, with corrected direct transport faces.
It must reproduce all six original diagnostic CSVs exactly. The fixed-native
records for this same region already exist; they are reused by symlink in
`regular-mol-compat-native-tail64-20261103/` and independently reverified,
without rerunning or modifying the fixed reference. This comparison isolates
the first new stage difference from the earlier ratio-normalizing tail capture;
it does not presume that any arithmetic control fixes the remaining excess.

The current compatibility tail observer completes with six exact original
CSV outputs. Candidate1296 momentum updates close within5.5500986835890775e-17,
hydro divergence within2.9103830456733704e-11, and3888 transport component
rates within7.275957614183426e-12. Paired momentum updates close within
2.152524981105092e-15, and direct directional transport increments agree
with independent stage-rate attribution within5.153379160175404e-20.
At the worst fine(64,103,109), x-momentum difference is8.846439856699528e-11
at25ns,1.385351442011358e-6 at37.5ns, and1.1002593854306415e-6 at50ns.
The25--37.5ns increase contains1.2868752641726217e-6 transport and
9.904475533403456e-8 x-normal hydro contributions. At event13, ghost(65,103,109)
uses candidate/native z factors0.4652599395758062/0.44880237037091714;
the x/y velocity differences4.86279370950915e-5/0.003045910248220096 m/s
shrink to-1.2830325274437726e-9/4.4941828036826337e-13 under native-factor
substitution. The first root sync changes the center donor z-momentum
discrepancy from2.9512290134867974e-14 to3.2883459614830897e-11.

`prior-comparison.json` compares the earlier normalized/uncapped observer
with the current ratio-preserving/capped one. The72 root states are identical
at initial and predictor stages. Predictor reconstructed face pressures differ
by at most8.731149137020111e-11 Pa; the next z-momentum states differ by at
most4.656612874555465e-15. This establishes a same-state primitive/face
roundoff difference before the later amplification. It is not a whole-run
single-variable test because the weight cap also differs.

A separate fixed-native self-sensitivity experiment now runs under
`native-pressure-association-20261103/`. An isolated Diffusion.cpp copy
reassociates the ideal-gas pressure multiplication from rho*R*T*sum(Y/W) to
(rho*T)*(R*sum(Y/W)), retaining the same sum order and all numerical inputs.
Only the nonrotating primitive-pressure assignment changes; original source,
binary and frozen reference remain intact. It will require exact initial CSV
identity and compare composite-cell50ns differences to the unchanged reference.
This is a roundoff-sensitivity diagnostic, not a new reference, accepted
correction, justification for tolerance changes or evidence of C2 completion.

The native pressure-reassociation diagnostic completes with an identical
initial CSV and unchanged physical equations/input files. Its composite-cell
comparison uses1179648 cells and the same SI conversion as the original
checker. At25ns x-momentum self-difference is4.9813052271571734e-8; at50ns
it is1.0524018055585682e-6 at fine(64,8,32), with two cells above the unchanged
1e-6 momentum limit. Density scaled5.449422827194646e-9 and total/internal
energy scaled7.473216349862826e-9/7.473425488783125e-9 remain below1e-8.
This proves sensitivity of this evolved maximum to an algebraically equivalent
native pressure evaluation; it does not prove candidate correctness, invalidate
all original gates, or authorize their relaxation. `comparison.json` retains
source/binary/reference hashes and all fields.

On38520 archived native primitive samples, `local-pressure.json` reproduces
the original CKPY pressure exactly and measures at most3 ulps, or
3.492459654808045e-11 Pa, from the alternative multiplication association.
This is a bound for those samples, not a global bound over every evolving
cell. Native paired16/32-grid calculations at the same25ns step and50ns
endpoint are running to quantify resolution dependence of this sensitivity.
A canonical reduced-reactor16/32-grid comparator is prepared in
`regular-mol-fixed-dt-spatial-20261103/`; it will reuse the newly generated
original-native references after their completion and combine with the verified
current64-grid case for fixed-dt spatial analysis. It has not been launched
or counted as passed.

The native16/32-grid pairs complete with byte-identical initial CSVs. All
source/binary and exported-output hashes are checked in `grid-comparison.json`.
At fixed25ns steps,50ns x-momentum self-differences are8.043051294426284e-11
(n16),9.487015392023457e-10(n32),1.0524018055585682e-6(n64). Every maximum
is adjacent to the coarse/fine interface; only n64 crosses the original1e-6
threshold. This finite-resolution result demonstrates strong resolution-dependent
sensitivity to pressure arithmetic; it is not a bound on every equivalent
implementation or a reason to change the reference gate.

The first canonical16/32 comparison launch stopped before solver execution
because the newly generated reference directories lacked reference_manifest.json.
The original failure/log is preserved. Both manifests were then constructed
from actual verified original-binary/source identities, exact expected input
and grid files, and hashes of every captured plotfile artifact; no field data
or contract hashes were substituted. A separate retry runs under
`regular-mol-fixed-dt-spatial-20261103/retry/`. Its prepared conservative
restriction analysis separates state and evolution-increment L1/Linf norms
at25/50ns, with the same dt on16/32/64. Reference failures remain visible;
no convergence result is yet claimed for that running candidate comparison.

The canonical reduced fixed-dt retry completes n16 in192.65s and n32 in
1493.57s. The n32 original25/50ns gates pass, with maximum x-momentum
differences1.596807934801936e-7/1.2348026744157803e-7. The n16 case fails
at25ns: density scaled1.3533807431649876e-8, x-momentum1.2507029330602326e-6,
total/internal-energy scaled1.9035044758689633e-8/1.903632031312775e-8.
Its50ns gates pass, with x-momentum5.024272582074274e-7. Both n16 maxima
are fine(16,20,22). The existing n64 case remains failed at50ns. Source363,
probe and report hashes are independently verified in
`regular-mol-fixed-dt-spatial-20261103/verification.json`; maximum gas-budget
residuals are8.16e-17/3.36e-17. These genuine canonical reference differences
are distinct from the smaller native-only n16/n32 association sensitivity.
The matched-time conservative state/increment spatial analysis subsequently
completes. `fixed-dt-convergence.json` contains160 rows, of which128 are at
positive times (25/50ns); all256 positive-time L1/Linf comparisons strictly
decrease from16/32 to32/64. This count includes both solvers, state and
evolution increment, and16 fields. Initial-state rows are excluded from that
strict-decrease claim. The analysis restricts eight fine cells onto each
same-level coarse cell, excludes covered coarse cells, and uses physical
volume weights. Temperature is volume averaged, not reconstructed from energy.

Candidate momentum-increment orders at25ns are L1/Linf1.395/0.0723 (x),
2.172/1.230 (y), and1.679/0.0585 (z); at50ns they are1.239/0.2207,
2.352/1.385, and1.632/0.2886. The fixed native solver gives nearly identical
orders. In particular,25ns z-momentum maximum spatial differences decrease
only from1.347640070534156e-4 to1.294065525079966e-4. These finite-resolution
self-differences are distinct from candidate/reference errors and do not
establish an asymptotic rate or complete C2. Native pressure-roundoff
sensitivity does not explain away the slow maximum spatial convergence.

`convergence-summary.json` retains the momentum rows and result hash;
`convergence-verification.json` independently verifies that hash, all18 input
CSV hashes,128 positive-time rows and256 strict norm decreases. Source and
binary identities remain those in `verification.json` and the individual
run records. MPI3 checks of the canonical n16/n32 cases are running under
`mpi-notification/`, with commands and outcomes in `mpi-record.json`.
They require six exact serial/MPI CSV pairs per grid and preservation of
the original reference outcomes, including the n16 failure; no MPI pass
is claimed before those checks finish.

The MPI3 driver subsequently completes (exit0), independently rechecked in
`mpi-verification.json`: all12 serial/MPI state/budget CSV pairs match exactly.
The underlying n16 checker exits1 with the same four25ns reference failures;
n32 exits0 with all original gates passed. Elapsed checker times are422.75s
and3271.73s. Driver success means successful comparison of these expected
outcomes, not12 fixed-reference passes. C2 remains open.

The n16/25ns canonical reduced maximum is traced with isolated observers in
`regular-mol-reduced-stages16-20261103/` and
`regular-mol-native-stages16-20261103/`. All six candidate state/budget CSVs
and all three native reference CSVs are byte-identical to their original runs.
The native observer driver used an incorrect exported-CSV directory after
successful CFD execution; the failure is preserved in `record.json` and
`notification/stderr.log`. `verify_exports.py` verifies the successful native
run against the correct exports without rerunning CFD, recording
`verification-record.json`. Native duplicate tile observations agree exactly.

Candidate1464 momentum updates close within5.53e-17; directional hydro
rate residuals are at most1.17e-10, and4392 direct transport rate components
close within3.64e-12. Paired candidate/native momentum increments close
within3.69e-16, and direct directional transport attribution within2.21e-20.
At fine(16,20,22), x-momentum discrepancy progresses from0 to-4.2364647273e-7
at12.5ns and-1.2507029331e-6 at25ns. The two intervals' x-direction hydro
contributions are-3.6364301079e-7 and-7.0891599229e-7; total transport
contributions are-5.9999545677e-8 and-1.1812193162e-7. Thus this early n16
maximum is predominantly hydro, unlike the previously observed n64 tail.

Paired donor replay reconstructs96 ghost observations within8.89e-16 velocity
residual. At ghost(17,20,22), event5, the z-direction common limiter is
controlled by H2O2: candidate/native factors0.21214476254/0.21461219970.
The x/y velocity differences are-3.1840560730e-4/1.0693793068e-3 m/s.
Substituting only native limiter factors in the candidate donor reconstruction
reduces those differences to-6.38e-13/-2.14e-14 m/s. This algebraic
counterfactual identifies local amplification, not an adopted correction or
proof that a species equation should be excluded from the native limiter.
An actual native reactor-input capture at donor(9,10,11) and its six neighbors
is now running to separate the upstream chemistry difference from interpolation.

The donor(9,10,11) native reactor capture completes in
`regular-mol-reduced-reactor16-20261103/`, preserving both original native
CSV exports. Actual native CPU replay reproduces all21 species/temperature
outputs exactly. On identical native inputs, canonical reduced/full maximum
temperature differences are6.3933e-9/7.6661e-9 K and maximum species
differences2.5326e-13/2.5244e-13 kg/m3. At the donor corrector, reduced
H2O2 error is1.38814109649e-16 kg/m3 (native event4, final corrector), compared with the observed half-time
donor difference6.94070501457e-17 kg/m3.

Existing isolated full-reactor controls reduce maximum temperature error to
4.8431e-11 K with native activation convention, and1.1369e-12 K with activation,
native component weights and equilibrium convention. The latter weights here
are uncapped diagnostic binaries, not an adoptable policy; the ongoing CFD
comparison uses the previously checked caller-capped variant. Activation-only
and weight controls must not be described as a weights-only intervention.
`controls.json` verifies all reused binary/source hashes.

`seed-counterfactual.json` is an algebraic sensitivity calculation, not a
conservative coupled update: replacing estimated half-step species/temperature
reaction errors in the candidate donor with exact-native replay errors changes
the z limiter from0.21214476254 to0.21464928125, near native0.21461219970.
This corrected result uses native event4 (final corrector); the earlier
first-iteration substitution was incorrectly described as the corrector and
is retained as `seed-counterfactual-first-iteration.json`, not evidence for
an unexplained large remainder. Identical-native-input chemistry replay explains
most of the local limiter difference, but this algebraic test is not a
conservative trajectory correction.

A bounded n16 trajectory comparison of existing canonical full, activation-only,
and caller-capped compatibility probes runs in
`regular-mol-reactor-controls16-20261103/`. Its initial driver failed before CFD
because the standalone checker lacked its tools PYTHONPATH; original failure
records remain intact. `retry_driver.py` supplies that path and writes separate
retry outputs, `retry-record.json` and `retry-notification/`. No result is counted
before the original field gates are inspected; compatibility variants remain
diagnostic even if they pass.

All three n16 reactor-control trajectories subsequently pass the unchanged
25/50ns field gates. `verified-summary.json` rechecks report and binary hashes.
The canonical full probe has x-momentum maxima4.7791399865e-9/1.9121484604e-9;
activation-only diagnostic4.5732630863e-8/1.8372792957e-8; caller-capped
compatibility diagnostic7.6357672962e-10/4.1406092444e-10. Elapsed times are
368.95/368.87/370.76s. These small-grid successes do not override any n64
failures or promote diagnostic mechanism conventions to production. In
particular, activation-only improves identical-input reactor agreement but
has larger evolved n16 momentum maxima than canonical full integration.

An isolated candidate reaction observer is running under
`regular-mol-candidate-reactor16-20261103/`. It records the actual reduced
reactor inputs and outputs at donor(9,10,11) and six neighbors for the first
three root reaction calls, and requires six exact original state/budget CSVs.
The prepared replay separates reactor differences on candidate inputs from
the native reactor's response to differing inputs; no decomposition result
is claimed until exact candidate replay and reference alignment are checked.

Actual candidate capture completes with all six original CSVs exact.
`regular-mol-candidate-reactor16-20261103/replay-comparison.json` replays all21
captured reduced-reactor outputs exactly. Fourteen first-step cases pair
candidate calls1/2 with native calls3/4: native calls1/2 initialize root/fine
history, whereas the candidate observer starts in the stepping routine.
Candidate call3 is the next root timestep and has no paired native observation
in this one-step native capture. The initial event-number alignment was rejected
and retained in `rejected-event-alignment.json`; it must not be used as evidence.
The corrected pairing is checked against common old composition, temperature,
duration and the native initialization flag/time.

At donor(9,10,11), final corrector H2O2 difference is1.38814100291e-16 kg/m3.
Its decomposition is1.38814109933e-16 from reduced versus native integration
on actual candidate inputs, and-9.64182401507e-24 from the native reactor's
response to differing inputs. Native SI/CGS conversion can add roundoff;
this is a measured local decomposition, not a universal error bound. It
supports investigating component error weights/coordinate choice before
changing the shared AMReX interpolation rule.

`regular-mol-full-capped-weights-20261103/` now isolates caller-capped native
component absolute tolerances in the full reactor while keeping the canonical
mechanism, activation and equilibrium constants unchanged. It reuses the
hash-verified isolated implementation, caps species/temperature absolute
weights at min(caller_atol,1e-15)/min(caller_atol,1e-18) in SI, and retains
rollback and endpoint admission. The campaign requires original n8 serial/MPI3
gates and exact outputs, then n16 gates before n64/50ns. This is a candidate
numerical policy under evaluation, not a production change or completed C2.

The canonical full-reactor component-weight control completes, but is not
adopted. `regular-mol-full-capped-weights-20261103/verification.json` rechecks
363 source identities, all recorded override-source hashes, n8 original gates
and eight exact serial/MPI3 CSV pairs. The n16 original gates also pass. At n64,
25ns x-momentum maximum is1.40451116471e-7 at fine(64,90,72);50ns fails four
original gates: density scaled1.15058474760e-8, x-momentum2.21378633114e-6,
total/internal-energy scaled1.57791226953e-8/1.57795629542e-8. Maximum x error
is now at fine(64,8,32). The bounded n64 checker takes24597.49s and exits1;
the campaign's exit0 only means it completed and recorded that expected
pass/fail classification. It is not reference qualification. Conserved gas
budget residuals remain below4.78e-17 at the two positive endpoints.

This candidate worsens the prior canonical full50ns maximum1.98320921284e-6
and does not justify a production tolerance-policy change. Its maximum cell
coincides with the earlier native pressure-association diagnostic maximum;
a plane-wide signed correlation analysis is being performed before assigning
any shared cause. The frozen native reference and original thresholds remain
unchanged.

The signed50ns x-momentum plane comparison completes in
`regular-mol-full-capped-weights-20261103/pressure-correlation.json`, with
all16384 fine i64 cells and independently rehashed native exports. Correlation
between candidate-minus-original and pressure-associated-native-minus-original
is0.5895711. Both maxima are positive at(64,8,32), but candidate-minus-associated
still has maximum1.8760054762e-6 at(64,121,97), with10 cells above1e-6.
Consequently this coincidence does not establish pressure arithmetic as the
whole cause; the associated native is never used as a replacement reference.

A separate isolated control in `regular-mol-native-pressure-order-20261103/`
changes only admitted-ratio primitive-pressure arithmetic to the native CKPY
CGS operation order, returning SI pressure afterward. Canonical mechanism,
full-reactor scalar tolerances, and every original gate stay unchanged; the
failed component-weight policy is not combined with this experiment. The
original primitive checks and finite-pressure guard remain. This is a
reference-arithmetic control, not a production change or claimed solution.

`primitive-check.json` validates864 archived same-input conversions: the
generic path remains byte-identical, and only pressure changes when admitted
ratios are requested. Among432 native-input cases, exact pressures increase
from140 to295 and mean absolute pressure error decreases from1.499e-11 to
8.421e-12 Pa; for432 candidate-input cases, exact counts increase152 to304
and means1.452e-11 to7.849e-12 Pa. Maximum pressure error stays7.276e-11 Pa
in both input groups. The campaign requires n8 serial/MPI3 exact outputs
and original n16 gates before the n64 comparison; no trajectory improvement
is inferred from these local means.

The native-pressure-order campaign completes in24949.02s for n64, with
original n8 serial/MPI3 and n16 gates passed and eight exact n8 CSV pairs.
`verification.json` rechecks363 source identities and report/binary hashes.
At25ns x-momentum maximum is1.09280339423e-7 at fine(64,26,22); at50ns it
is1.46727344814e-6 at fine(64,72,33). Density now passes, but x-momentum
and total/internal-energy scaled1.04804633984e-8/1.04807549407e-8 still fail.
Budget residuals are below4.02e-17. The pressure-order candidate reduces the
prior full-reactor maximum but does not qualify and is not adopted.

A small isolated reduced-CVODE seed-history control in
`regular-mol-fixed-reactor-seed-20261104/` keeps the admitted initial EOS seed
in RHS/Jacobian callbacks instead of replacing it with the last recovered
value. It changes no mechanism or tolerance. Across21 native n16 inputs,
21 actual candidate n16 inputs and21 native n64 inputs, maximum H2O2 errors
remain essentially unchanged (about1.4--1.6e-16 kg/m3), while maximum temperature
errors increase. Both variants use630 total internal steps on each21-case set.
`comparison.json` records source/binary/input hashes. There is no demonstrated
improvement and no full CFD launch or production adoption for this control.
An absolute-tolerance refinement replay on actual candidate inputs now checks
whether the remaining trace-species difference converges before choosing any
further large trajectory experiment.

The21-actual-input absolute-tolerance sweep completes at factors1,0.1,0.01,
0.001,0.0001, with relative tolerance fixed1e-10. Reduced H2O2 differences
from the tightest run decrease8.763e-17,2.102e-17,1.449e-18,3.125e-19 kg/m3;
but differences from native remain around1.57--1.76e-16 kg/m3. Full-coordinate
H2O2 differences from the tightest run also decrease (1.328e-15 through8.960e-18),
while native differences are nonmonotone and end at1.722e-16. Total internal
steps increase630 to2679 (reduced),2494 to44349 (full). The tightest result
is a numerical comparator, not exact truth, and species-wide differences
are not uniformly monotone. No stricter tolerance is adopted from this sweep.
`verification.json` rechecks the experiment's binaries, inputs and result hashes.

This separates integration-error refinement from fixed-native agreement:
reducing estimated integration error alone does not remove the measured
canonical/native trace-species discrepancy. Earlier references to the
"integrator-side" difference must include reaction-rate conventions, not
just CVODE coordinates or error weights. The existing Cantera importer
`tools/ingest_cantera_mechanism.py` explicitly normalizes activation energy to
J/kmol using pinned Cantera input_data; `mechanisms/lidryer.json` pins that
source and conversion. A native convention cannot silently replace those
canonical coefficients or be counted as canonical qualification. Any supported
native-convention selection would require explicit provenance and a distinct
selection identity, plus original fixed-reference tests; no such production
selection has been added in this diagnostic turn.

The source-selection audit finds an existing explicit diagnostic generation
path, `tools/prepare_lidryer_native_activation.py`; no new mechanism selection
or alternate canonical baseline is needed or introduced. Its activation
validation is strengthened from a sorted nonzero coefficient comparison to
reaction-index, stoichiometric-equation and high-/low-pressure branch checks,
including zero activation energies. The manifest now records23 branch mappings.
Swapping coefficients between reactions or pressure branches is rejected even
when the coefficient multiset remains unchanged. This changes verification,
not rates or source conversion.

`build/native-activation-mapping-20261104/` contains before-images, the focused
increment patch, failures and verification receipt. Seven new Python tests pass;
the actual pinned header validates23 branches. Regeneration with pinned
Cantera3.2.0 (4a8358e) produces YAML/JSON/Fortran byte-identical to the previous
verified diagnostic generation. Canonical LiDryer provenance/regeneration passes.

A separate test-only repair updates stale importer expectations: six failures
in the original12-test run expected currently supported PLOG, Chebyshev, SRI,
Tsang and chemically activated classifications to reject. The corrected tests
check their classification and still reject incomplete parameter records;
custom reaction orders and explicit reverse overrides remain rejected. All13
importer tests pass in the pinned Cantera environment, and both registered
CTest entries pass. Syntax, project contract and diff checks pass. These are
focused tool/provenance results, not full configuration, installed-binary,
fixed-reference trajectory or physical-validation completion. C2 remains open.

### 2026-11-04 explicit equilibrium reference-concentration convention

The elementary-reaction library now exposes
`equilibrium_reference_gas_constant` in J/(kmol K), defaulting to the existing
thermodynamic constant. It changes only the standard concentration factor in
Kc. Arrhenius activation units, the EOS and pressure-dependent forward-rate
pressure remain unchanged. This makes the previously isolated native-R
experiment expressible without rewriting the core source. There is no new
application input, mechanism importer policy or restart-selection contract;
the diagnostic convention is not enabled in canonical applications.

The dirty-main base is `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
`build/equilibrium-reference-policy-20261104/` preserves the kinetics/CMake
before-images, commands, hashes and logs. A clean GNU13.3 Release build with
MPI/SUNDIALS/spray disabled passes four focused tests: elementary kinetics,
the new equilibrium convention test, pressure-dependent kinetics and the
full H2O2 Jacobian. The new test covers net stoichiometry -1/0/+1, unchanged
forward rates, finite-difference Jacobians, invalid constants and unrepresentable
reference concentrations; its separate floating-point-trap run also passes.

The independently rebuilt previous kinetics object and dependent mechanism/test
produce byte-identical basic-test output. Forty-two frozen LiDryer states also
produce byte-identical default rates and density Jacobians against the
hash-verified prior replay binary. The explicit convention produces42 finite
rows. These are fixed-state checks, not complete reactor or trajectory parity.
LiDryer and H2O2 provenance/regeneration checks pass. Receipts are
`default-comparison.json` and `rate-replay.json`.

The first clean MPI/SUNDIALS build completed in500.6 seconds; all five focused
tests and the n8 serial original-reference gates passed. MPI3 failed before
calculation: CMake selected Intel `libmpi.so` despite the OpenMPI compiler
wrapper, through inherited library-search settings. Its failed receipt and
logs remain in the artifact directory. A fresh `clean-openmpi/` build removes
the inherited library/include/prefix search paths and explicitly binds the
system MPI and hwloc libraries. MPI runtime verification remains pending;
the first build is not MPI qualification. `audit_clean_mpi.py` will verify
source/binary hashes, linked dependencies, original gates and eight exact
serial/MPI CSV pairs after completion. C2 trajectory/convergence and physical
validation remain unresolved.


The sanitized clean OpenMPI/SUNDIALS rebuild subsequently completed. Independent
`clean-openmpi/verification.json` checks the frozen source manifest, executable
hashes and dynamic-library resolution (OpenMPI40, no Intel MPI or unresolved
libraries). All five focused tests pass; n8 serial/MPI3 original-reference
gates pass with eight byte-identical state/budget CSV pairs. This supersedes
the pending runtime check above, while preserving the first failed attempt.

`reactor/verification.json` then exercises21 captured candidate reactor inputs
with both full and reduced coordinates. Default outputs, including internal
step counts, are byte-identical to the hash-verified prior binaries. Changing
only the equilibrium reference constant to8314.46 leaves full/reduced step
totals at2494/630 and barely changes the native discrepancies: reduced maximum
H2O2 error moves from1.5701862e-16 to1.5695244e-16 kg/m3 and maximum temperature
error from8.8980414e-9 to8.8955403e-9 K. Full-coordinate maximum temperature
error remains approximately9.53e-9 K. Thus this isolated convention change
does not resolve the captured discrepancy; it is not evidence for adopting
it alone or for rerunning the long64-grid case. Activation-energy conversion
and integration effects remain separate controls. The library option enables
those controls without core-source substitutions; default application behavior
and original reference limits are retained.

### 2026-11-04 activation/equilibrium factorial reactor controls

`build/equilibrium-reference-policy-20261104/activation-reactor/` compiles the
hash-verified regenerated diagnostic mechanism against the current clean core.
The generated source is unchanged; the new public reaction field selects the
equilibrium convention without source substitutions. The21 actual candidate
inputs and scalar caller tolerances are unchanged. With activation conversion
alone, full-coordinate maximum native temperature/species differences fall to
4.8431e-11 K/2.4993e-15 kg/m3. Reduced-coordinate H2O2 difference instead rises
to1.2638e-15 kg/m3. Adding the native equilibrium constant changes these only
slightly. `verification.json` records all four controls, commands and hashes.

An absolute-tolerance sweep (factors1,0.1,0.01,0.001,0.0001; relative tolerance
fixed at1e-10) uses activation plus native equilibrium conversion. Full-coordinate
H2O2 difference from the tightest control falls from1.3285e-15 to8.9077e-18
before the final factor, but difference from the fixed native output grows
from4.4864e-17 to1.2882e-15. Reduced-coordinate difference from its tightest
control falls from8.7620e-17 to3.1248e-19, while its native difference approaches
1.2917e-15. Step totals rise from2494 to44349 (full), and630 to2679 (reduced).
Thus close agreement at the original full tolerance is not evidence of a
converged reference trajectory; tightening the candidate alone does not
monotonically approach the fixed native output. The tightest control is not
assumed exact. At that setting full/reduced H2O2 differs by3.9686e-18 kg/m3,
whereas their maximum temperature difference remains7.6034e-9 K. Temperature
recovery and coordinate coupling therefore still require separate inspection.
The fixed-reference solver/tolerances have not been changed and no long-grid
candidate is adopted from these frozen-input results. Receipts:
`tolerance.json` and `coordinate-comparison.json`.

### 2026-11-04 inverse-temperature stopping and coordinate difference

The preceding tight full/reduced states are replayed through the actual EOS
inverse with their affine final energy and their respective final species.
`build/equilibrium-reference-policy-20261104/eos/verification.json` records252
cases:21 states, two coordinates, three temperature seeds (own, other,
common2050 K), and default/converged pinned inverse policies. With default
stopping, each own seed is returned exactly; exchanging seeds exchanges the
approximately7.6034e-9 K difference. A common seed or the converged inverse
reduces the full/reduced maximum difference to6.2300e-11 K. The source condition
is the pinned default `abs(step)<1e-6` exit before applying the last update;
this is intended reference behavior, not a newly identified solver defect.
Postprocessed composition differences still contribute the smaller residual.

`converged-inverse/verification.json` repeats the entire21-input reactor
integration with only the existing converged-inverse option enabled, using
the activation/native-equilibrium diagnostic mechanism. Full-coordinate native
temperature difference increases to7.5870e-9 K at the original caller tolerance;
its H2O2 difference remains4.4864e-17 kg/m3 with2494 internal steps. Tightening
absolute tolerance by1e-4 gives full/reduced native temperature differences
8.3371e-9/8.3264e-9 K and H2O2 differences1.2882e-15/1.2917e-15 kg/m3.
Thus enforcing a more converged EOS inverse removes seed dependence but does
not improve fixed-default-native trajectory agreement. This control is not
adopted in the canonical path. The remaining native/candidate reactor mismatch
must be assessed with native integrator and inverse histories, not interpreted
as a requirement to tighten the candidate alone. No reference gates, generated
mechanism sources, canonical tolerances or EOS defaults were changed.

### 2026-11-04 paired actual reactor RHS histories

The existing hash-verified native observer binary was run on the21 captured
candidate inputs, producing endpoints byte-identical to their saved native
replay. `build/equilibrium-reference-policy-20261104/native-trace16/` records
2536 internal steps,2753 RHS calls and7 error-test failures; per-case RHS
counts agree with CVODE's reported statistics. This reuses the original
native solver and tolerances and does not modify the fixed CFD reference.

A separate observer build of the current full-coordinate candidate with the
activation/native-equilibrium diagnostic mechanism also preserves all21
unobserved endpoints byte-for-byte. Its2711 RHS calls are two fewer per case
than native; no case has matching call counts. Initial temperature differs
by at most9.0949e-13 K, and initial density coordinates by6.9389e-18 kg/m3
(including CGS/SI conversion). Thus ordinal callback pairing cannot be used
as a same-time trajectory comparison. The candidate retains its scalar
caller absolute tolerance, while native uses its species/temperature vector
weights; the adaptive-history difference is not by itself evidence of a
Jacobian or EOS defect. `candidate-trace16/verification.json` stores the
observer source, commands, hashes and first-positive-time comparisons. Earlier
native-vector-weight controls remain distinct and must not silently relax a
stricter caller tolerance. Next comparisons must evaluate both RHS functions
at identical captured states or explicitly match weight policies.

### 2026-11-04 same-state RHS factorial audit over native trial history

`build/equilibrium-reference-policy-20261104/rhs-trace/verification.json`
evaluates all2753 actual native RHS states through the native kernel and four
candidate controls, with identical time, affine energy, temperature seed,
species and external forcing. Every recovered temperature is exactly equal
to native in all four controls. The maximum norm-scaled species RHS difference
is4.9965e-6 for canonical chemistry,4.9971e-6 with only native equilibrium
reference concentration,6.4483e-10 with activation conversion alone, and
5.9421e-14 with both. Maximum relative temperature RHS difference similarly
falls to8.9261e-14 with both conventions. Native SI/CGS roundtrip reevaluation
has a maximum species norm difference8.3364e-16; candidate roundoff remains
larger than this floor, so bitwise RHS equivalence is not claimed. All controls
use the current clean core and explicit equilibrium field; generated source
and fixed-reference inputs remain unchanged.

This narrows the remaining full-coordinate endpoint mismatch to adaptive
integration/weighting and residual floating-point differences for this
captured set; it does not establish global trajectory parity. The21 scalar
caller absolute tolerances span1.0314544e-15 to1.0377672e-15, so in this set a
caller-capped native species weight of1e-15 would not relax them. A separate
weight intervention can therefore compare adaptive histories without changing
reaction conventions or violating caller bounds. General callers may request
smaller tolerances, which must remain binding. No production weight policy
has been changed by this audit.

### 2026-11-04 caller-capped vector-weight library implementation

The isolated capped-weight observer on the same21 actual inputs matches all
native per-case RHS counts:2753 total calls and2536 internal steps. Maximum
native H2O2 endpoint difference is7.5997e-23 kg/m3 and temperature difference
1.1369e-12 K. Plain and observed binaries produce byte-identical outputs.
This is local reactor evidence, not a new64-grid success; the earlier evolved
native-compatible64-grid candidate still failed the original momentum gate.

The full-coordinate library now accepts optional `absolute_tolerances(:)`
(species density entries followed by temperature), using elementwise
`min(atol,absolute_tolerances)` and retaining the scalar path when omitted.
Invalid sizes/nonfinite/nonpositive values reject before mutation. Four
focused CVODE tests pass after rebuilding the GNU/OpenMPI/SUNDIALS configuration;
the extended full test also passes after adding vector-specific failed-solve
rollback/context-release coverage. A request exceeding the caller cap produces
exactly the same states, temperature and step count as an explicitly capped
request. Invalid input and zero-duration shape checks pass.

`vector-api/verification.json` independently recompiles replay programs against
the changed library:21 default outputs match the pre-change default executable,
and21 explicit vector outputs match the isolated capped implementation exactly.
Before-images are `sundials_before_vector.F90` and
`full_test_before_vector.F90`. Existing earlier build/source manifests describe
the pre-vector revision and must not be treated as the current source identity;
the new receipt hashes the changed source/test and replay executables. No
application selector, checkpoint contract, canonical mechanism or default
integration policy changes here. MPI/application regression for this additional
API change is not included in the earlier pre-vector n8 qualification.


The post-vector fresh GNU/OpenMPI/SUNDIALS build now passes all five focused
tests and n8 serial/MPI3 original-reference gates. Independent
`vector-clean/verification.json` verifies current source/binary identities and
clean MPI dependencies; eight state/budget CSV pairs are byte-identical across
ranks. `vector-clean/default-application-parity.json` additionally records eight
byte-identical comparisons to the pre-vector serial outputs. These replace the
pending MPI/application check for this API revision, without extending the
scope to restart-option selection or64-grid parity.

`vector-multiset/verification.json` exercises the public API on three independent
captured21-state sets (n16 native inputs, n16 candidate inputs, n64 native donor
inputs). With the explicit activation/equilibrium diagnostic convention, maximum
temperature error is1.1369e-12 K in each set; maximum H2O2 error is respectively
2.9210e-24,7.5997e-23,6.2039e-24 kg/m3. Maximum species-wide errors remain at
most6.6614e-16 kg/m3. These63 local replays do not supersede the failed evolved
n64 native-compatible trajectory. Existing paired n64 face records show that
fine ghost velocity differences first grow at hydro event5 while adjacent
interior differences remain around1e-13 m/s; further work returns to the
coarse/fine interpolation inputs and their floating-point conditioning.

### 2026-11-04 remaining compatible64 donor z-momentum attribution

The existing native-compatible64 capture is revisited at the first amplified
ghost (65,103,108), hydro event5. Its root donor center is(33,52,54). The shared
z-direction slope factor is controlled by z-momentum, not H2O2: candidate
0.0315105411093 versus native0.0323245692179. The small right difference is
-8.644881e-13 versus-8.870848e-13 kg/m2/s. This factor change produces observed
velocity differences(2.828909e-6,1.608913e-4,3.3300e-14) m/s. Substituting only
the native factors in the candidate stencil leaves all three velocity errors
below8e-14 m/s. This is an algebraic diagnostic, not an evolved intervention.

`build/compat-root-z-attribution-20261104/verification.json` compares first-root
step face contributions at the center and its two z neighbors. Their actual
half-time z-momentum differences are(-7.840509e-15,-1.789832e-14,1.475615e-14)
kg/m2/s. Integrating directional face differences gives dominant z hydro
contributions; transverse hydro and transport terms are much smaller. Residual
closure is at most3.293e-16 kg/m2/s, bounded explicitly at5e-16, so this is a
dominant-term attribution rather than exact closure. Native divergence is
formed in original CGS before conversion; subtraction order and recorded
face precision matter at these tiny differences. Stage3 z-momentum flux
(candidate minus native) spans approximately-2.47e-10 to4.07e-10 Pa across these
faces. The first-growth mechanism in this remaining compatible case therefore
differs from the earlier canonical H2O2-controlled case. No limiter threshold,
reference tolerance or production flux has been altered from this analysis.

### 2026-11-04 root z-face frozen-input pressure decomposition

`build/compat-face-replay-20261104/verification.json` replays18 root z faces
(three donor cells, two sides, three initial hydro events) through the current
MOL kernel. Candidate-input reconstructed pressures and momentum fluxes match
the recorded evolved calculation exactly. Replaying the corresponding native
primitive stencils in SI separates the input term from the residual kernel/unit
conversion term. The initial event agrees exactly; at the third event the
observed z-momentum flux difference reaches4.0745e-10 Pa, while the same-physical-
input kernel/conversion residual is at most1.4552e-11 Pa. This is not bitwise
native-kernel equivalence because native CGS inputs were converted to SI.

Ninety frozen one-field substitutions are recorded in `substitution.json`.
At event3, replacing only the candidate stencil pressure with the native
pressure reduces the maximum z-flux residual to1.4552e-11 Pa. Replacing only
velocity, density, species or sound speed individually leaves the maximum at
4.0745e-10 Pa. These substitutions are not additive causal shares. They localize
the dominant discrepancy to cell-pressure inputs before face reconstruction,
rather than demonstrating a defective Riemann formula. The next distinction is
between evolved thermodynamic-state differences and primitive-conversion
arithmetic on identical conserved states. No production pressure formula,
limiter or reference gate was changed.

### 2026-11-04 conserved-state versus primitive-arithmetic pressure split

`build/compat-primitive-replay-20261104/` evaluates432 paired root states from
the current compatible64 captures through both native `pc_ctoprim` and the
ratio-preserving candidate converter (864 inputs per solver). Both solvers
exactly reproduce their own recorded primitive fields. Maximum same-input
cross-solver pressure difference is7.2760e-11 Pa and temperature difference
1.1369e-12 K, including SI/CGS conversion roundoff.

At root(33,52,54), third hydro evaluation, the observed pressure difference is
4.0745e-10 Pa. A native-anchored decomposition assigns3.7835e-10 Pa to the
conserved-input difference and2.9104e-11 Pa to same-input primitive arithmetic
and unit conversion. Temperature arithmetic agrees exactly there. Thus a
pressure-formula-only change cannot remove the dominant local difference.

`substitution/verification.json` records1152 frozen native-kernel interventions
at72 third-event root states. At the center, substituting all native species
densities reduces pressure residual to-2.9104e-11 Pa, while substituting density,
momentum, total energy, held energy or temperature seed separately does not.
Among individual species, substituting N2 reduces it to2.9104e-11 Pa; H2 alone
leaves3.3469e-10 Pa and other species leave3.7835e-10 Pa. These inconsistent
one-field controls are sensitivity tests, not additive causes or conservative
solver updates. They identify bulk inert-species accumulation/transport roundoff
as the next local quantity to trace, rather than residual H2O2 reaction-rate
mismatch. All-native substitutions reproduce the reference pressure exactly.

### 2026-11-04 inert N2 update and density-coordinate control

`build/compat-inert-n2-20261104/verification.json` follows N2 through the first
root update at four cells around(33,52,54). At the center, old and predictor
N2 agree exactly, while the first reactor endpoint differs by6.0368e-16 kg/m3.
The affine N2 endpoints computed from recorded transport/hydro sources agree
in double precision. Candidate and native integrated endpoints differ from
that affine value by+2.7756e-16 and-3.2613e-16 kg/m3 respectively. The resulting
endpoint difference closes exactly algebraically. Native source reconstruction
retains divergence-order/recording roundoff, so residuals are not uniquely
assigned to CVODE without the actual reactor input trace. Neighboring source
endpoint differences are at most6.94e-18 kg/m3, much smaller than the endpoint
residuals. This is a bulk inert coordinate effect, not an N2 chemical rate.

An isolated full-reactor control in `build/compat-cgs-coordinate-20261104/`
advances density coordinates in CGS while retaining SI RHS evaluation. Species
absolute tolerances and RHS entries are scaled consistently; analytic density
columns and the finite-difference temperature column use the corresponding
chain rule. The diagnostic requires explicit vector tolerances and does not
modify the standard library source. Across the existing three21-state sets,
step totals stay2508/2536/2478. On the n64 native-input set, maximum species
error falls to4.1633e-17 kg/m3 and N2 error to1.3878e-17 kg/m3. On the two n16
sets, maximum N2 error remains5.8981e-16 kg/m3. Temperature maxima remain about
1e-12 K; H2O2 error stays below7.7e-23 kg/m3. Coordinate scaling therefore
helps one set but does not establish a general correction; no CFD rerun or
production adoption is justified from this control alone. Actual target-cell
reactor inputs and SI/CGS conversion order remain to be distinguished.

### 2026-11-04 target-cell reactor reconstruction and widened control

`build/compat-target-reactor-20261104/verification.json` reconstructs first-
reacting-update inputs from stored root hydro/transport stages and imported
initial reaction history, including the predictor energy reset. The current
public full/vector API with explicit activation/equilibrium convention reproduces
all species and temperatures exactly at(33,52,53..56). On those same inputs,
the isolated CGS-coordinate reactor makes N2 endpoints equal to the actual
native captured state in all four cells, including the6.0368e-16 discrepancy
at the limiter-controlling center. This supports coordinate-roundoff sensitivity
at the observed first amplification; it does not qualify an evolved correction.

The separate widened record `build/compat-target-reactor-wide-20261104/` covers
72 captured root cells. Species replay is exact in all72, but temperature is
exact in48; the remaining comparison includes post-reactor temperature recovery
and differs by at most1.4530e-5 K. The initial whole-state exactness assertion
therefore failed and was replaced by explicit species/temperature counts, not
an enlarged tolerance or a claim of whole-state reproduction. On all72 species-
exact cells, CGS coordinates reduce mean absolute N2 difference from2.4787e-16
to2.5732e-17 kg/m3 and produce31 exact native N2 endpoints. Maximum N2 difference
increases from6.8695e-16 to7.4246e-16 kg/m3. Thus the targeted improvement is not
a uniform maximum-error improvement; exceptional cells and native coordinate
conversion must be checked before any long coupled candidate is promoted.
No production source or reference tolerance changes were made in this study.

### 2026-11-04 exceptional CGS cells: native replay and initial roundtrip

`build/compat-target-reactor-wide-20261104/native/verification.json` runs the
unchanged native reactor on all72 reconstructed candidate inputs converted to
CGS. All72 native/candidate step counts agree. The two largest discrepancies,
(34,12,54) and(34,53,54), are reproduced by native itself: both native and
CGS-coordinate candidate differ from the captured native trajectory by
-7.4246e-16 kg/m3 in N2 and agree exactly with each other in N2. Across the
whole set, same-input candidate/native N2 difference is at most1.3878e-17 kg/m3.
The exceptional error therefore follows input differences rather than a
candidate-only integration defect in this control.

`initial-roundtrip/verification.json` keeps candidate-derived forcing but
substitutes exact native initial fields. Restoring only initial N2 in original
CGS yields70/72 exact captured N2 endpoints; replacing all initial species or
all initial thermodynamic fields gives the same count. The two exceptional
cells retain their full discrepancy. Replacing temperature alone gives46 exact
cells. This distinguishes ordinary SI/CGS initial-state roundtrip effects from
the two remaining forcing-sensitive cases. Native input IDs were corrected to
integer formatting after an initial parser failure; only completed runs are
included. Further diagnosis targets their actual nonreacting source assembly,
not a blanket tolerance reduction or a change to the fixed reference.

### 2026-11-04 exceptional cells: reconstructed source controls

`build/compat-target-reactor-wide-20261104/source-roundtrip/verification.json`
reconstructs native species forcing from recorded direction-specific face fluxes
with exact native initial fields. Separate hydro/transport divergences, combined
face divergence, and substitution of all species sources each retain70/72 exact
N2 endpoints and the7.424616477180734e-16 kg/m3 maximum discrepancy. Substituting
the reconstructed energy source yields no exact N2 endpoints and a maximum
6.661338147750939e-16 kg/m3 difference. These controls do not establish the actual
reactor energy input, so they do not attribute the two exceptional cells to a
specific source term.

A separate observer in `build/compat-exception-native-reactor-20261104/` records
actual reactor inputs/outputs at the two exceptional cells and four root donor
controls. It replaces only a locally linked React.cpp object; the pinned source
and executable remain unchanged. The capture driver requires initial and first-
step CSV SHA256 equality to the original reference before accepting the trace.
The prepared replay separates exact actual input, reconstructed forcing, actual
energy with reconstructed species forcing, and actual species forcing with
reconstructed energy. At this entry's creation the capture is running; no capture
or replay pass is claimed. LiDryer provenance, Python syntax and project contract
checks pass; the initial project-contract invocation omitted the required
`--source` argument and was rerun successfully with `--source .`.

### 2026-11-04 actual first-MOL input and source-assembly attribution

The observer completed with exit0 (959s native calculation, about45s export).
Both exported CSV hashes equal the original reference. `replay.json` in
`build/compat-exception-native-reactor-20261104/` reproduces N2 endpoints exactly
at all six observed cells with actual native inputs. Keeping exact initial
fields but candidate forcing retains the two exceptions. Substitution controls
in `species_interventions.json` show that actual N2 forcing alone removes both
exceptions; replacing any of the other eight species individually does not.
Conversely, retaining only candidate N2 forcing among actual species forcing
restores both discrepancies. Actual energy forcing alone does not remove them.
The CGS endpoint difference converted to SI is7.453889935837843e-16 kg/m3;
converting each endpoint to SI before subtracting gives the earlier
7.424616477180734e-16 value. These are different rounding orders, not different
runs or relaxed criteria.

Pinned `Source/Advance.cpp` constructs the first MOL reacting source by forming
corrected state, applying `LinComb(1/dt, new, -1/dt, old)`, and subtracting the
previous reaction source. Later MOL iterations use an average of source terms.
The earlier face reconstruction incorrectly assumed that later-iteration
average also represented first-iteration floating-point arithmetic. With the
captured native predictor and face fluxes, `source_assembly.json` reproduces
actual N2 forcing exactly at all six cells. N2 initial reaction source is zero.
`wide_n2_assembly/verification.json` expands the same assembly to72 cells: native
replay with exact native initial fields, candidate forcing for other components,
and the reconstructed native N2 source now gives72/72 exact N2 endpoints.
This identifies a source-assembly rounding contribution; it does not establish
an evolved candidate correction or physical validation.

`candidate_assembly/verification.json` applies state-difference arithmetic to
candidate predictor and stage source values. With exact native initial fields,
both SI and CGS assembly retain70/72 exact N2 endpoints and the same maximum
error. With candidate initial roundtrips they give46/72 exact endpoints. Thus
changing this formula alone has not removed the exceptional local discrepancy:
upstream source differences can still select different rounded corrected states.
A separate isolated split-module control in
`build/first-mol-source-assembly-20261104/` tests first-MOL assembly in evolved
flow, using the previously recorded explicit-compatible library artifacts. It
runs n8 first and proceeds to bounded n64 only if the original n8 gates pass.
It is not a test of the new standard CVODE API, an adopted production change,
or C2 completion. At this entry's creation that control is running.

### 2026-11-04 bounded evolved first-MOL assembly control

`build/first-mol-source-assembly-20261104/record.json` completed with exit0.
The isolated split-module change preserves averaged fluxes and changes only
first reacting MOL source assembly; it uses SI arithmetic and the previously
recorded explicit-compatible chemistry artifacts. The n8 reference gates pass
(223s), followed by all original n64 gates at0/25/50ns (25022s). At50ns maximum
absolute x-momentum difference decreases from1.10025938543e-6 to
9.603001674893252e-7 kg/(m2 s), about12.7 percent. The maximum moves to fine
cell(64,26,109). Maximum density difference is2.3558338613227647e-10 kg/m3,
total-energy difference7.408757810480893e-4 J/m3, and temperature difference
1.3032974948146148e-6 K. No field/budget violations remain in this bounded run.
`comparison.json` asserts unchanged source/case contracts, bundle, checker,
reference, timestep mode and chemistry/transport settings against the preceding
compatible control. This pass is close to the momentum limit and does not
prove longer-time agreement, convergence, standard-source API behavior or
physical validation. Production source has not adopted this arithmetic change.
An MPI3 exact-CSV control is running under `mpi-notification/`; it requires n8
original gates and eight exact CSV pairs before bounded n64 and six exact pairs.

### 2026-11-05 MPI equality and current-source integration

**MPI evidence invalidated by the subsequent compile-path audit below.**
`build/first-mol-source-assembly-20261104/mpi-record.json` completed with exit0
and reported `bounded_mpi_exact`, with eight n8 and six n64 matching CSV pairs.
However, its purported MPI object lacks the MPI compile definition and calls
the serial integrator. These were duplicate serial executions under mpiexec,
not MPI verification. Their485s/52249s elapsed times cannot establish MPI
performance or MPI equality. The serial bounded result remains valid.

The same first-MOL source assembly is now applied to the working
`src/amr/amr_patch_tree_split_3d_mod.F90`. Only the first reacting MOL iteration
forms a corrected state, includes the old reaction estimate, recovers source
by reciprocal-timestep state difference, and subtracts that estimate. Later
iterations retain source averaging. The predictor already includes its energy
reset, so the separate first-iteration reset addition is excluded from this
reacting branch. Averaged face fluxes and nonreacting updates are preserved.
The exact previous file is in
`build/first-mol-current-source-20261105/split_before.F90`.

A new GNU/OpenMPI/SUNDIALS build under that directory verifies current source,
standard chemistry settings, focused split/subcycle/checkpoint/history/chemistry
tests and n8 serial/MPI3 original-reference gates. Its first invocation failed
before compilation because extra target names followed CMake's `-j4` option;
that command and logs are preserved in `attempt1/`. The retry puts target names
before the parallel option. At this record's creation the retry is running;
no standard-source test pass is claimed. LiDryer provenance, Python syntax,
project contract and whitespace checks have passed. C2 remains open.

### 2026-11-05 genuine MPI path correction and evidence invalidation

The fresh current-source build completed11 focused tests and both standard n8
reference comparisons successfully, but its final exact serial/MPI assertion
failed. At the final snapshot, maximum x-momentum serial/MPI difference was
8.564364495367016e-13 and temperature difference6.1390892369672664e-12 K.
The separately implemented MPI split module still used source averaging on the
first MOL iteration. It now has the same corrected-state source construction
with owner-local guards as the serial module.

The discrepancy also exposed a flaw in the older isolated compatibility build:
`build/regular-mol-composition-compat64-20261102/mpi.o` was compiled from raw
preprocessed source with `-cpp` but without `-DPELEF_REFERENCE_MPI`.
`nm -u` shows serial advance calls and no MPI advance dispatch. Its downstream
MPI equality claims, including the14 pairs above, are withdrawn. Existing
artifacts and logs are retained. The fresh CMake-built MPI probe has the correct
MPI dispatch; its failed exactness assertion was therefore useful evidence.
See `build/first-mol-current-mpi-source-20261105/mpi_path_audit.json` for object
hashes and symbol evidence. Future isolated MPI probes must retain preprocessor
definitions and verify actual MPI dispatch, not infer it from mpiexec alone.

A new fresh build in that directory is running both corrected source paths,
the same11 focused tests and standard n8 serial/MPI comparisons. No new MPI
pass is claimed until completion. The prior current-source test records are
preserved as failed overall, despite their individual tests/reference gates
passing. No reference tolerance or golden output was changed.

### 2026-11-05 current-source genuine MPI verification

`build/first-mol-current-mpi-source-20261105/audit.py` passes: all11 focused tests,
standard n8 serial and MPI3 reference gates, and eight exact CSV pairs. It also
verifies the unchanged source manifest and the compiled object's MPI dispatch
and MPI initialization symbols. The receipt is `verified.json`. This replaces
the invalid isolated-MPI evidence only for the verified current-source n8 scope;
it does not restore the old n64 MPI claim.

The first-MOL source assembly is implemented in both serial and owner-local MPI
split modules. An additional application run is checking nonreacting thermal
slip behavior and reacting MOL fixed/selected common/recursive restart across
serial/MPI1/MPI2/MPI3. Its copied hotspot fixture explicitly enables MOL, PLM
and PeleC fluxes while retaining chemistry and transport. Separately,
`build/first-mol-public-api-compat-20261105/` builds a current-source compatibility
probe using public full-CVODE/vector-tolerance and equilibrium-constant APIs,
plus the previously generated native-activation mechanism. It preserves CMake's
preprocessor definitions and asserts compiled MPI dispatch. This avoids reusing
the old modified CVODE/kinetics objects. Its first comparison is n8 serial/MPI3;
no completion of these additional checks is claimed here.

### 2026-11-05 reacting restart exposes diagnostic scope mismatch

The additional application build and both nonreacting thermal-slip tests pass.
The common-step reacting MOL restart checker passes continuation-output
comparisons through MPI2-to-MPI3 and MPI-to-serial continuation, but fails exact
checkpoint equality after serial writes its checkpoint. It stops before the
selected/common and recursive cases finish, so the suite is not passed.

A preserved MPI2 rerun in
`build/first-mol-current-mpi-source-20261105/checkpoint-diagnosis/` differs from
the serial checkpoint in exactly one field on one line: `maximum_reflux` is
2.727626953818670669e-4 serial versus3.328653373541728142e-4 MPI. All other
checkpoint lines match. The MPI coupled-reflux path appends held internal energy
to conserved components and previously reported the maximum over that augmented
vector; serial reports the conserved-component maximum.

The MPI reflux module now retains its public interface and forwards to a private
core with an explicit diagnostic-component count. Coupled synchronization still
updates gas and held energy together, but appending the held field cannot change
the public conserved maximum. The reflux-budget regression now supplies a large
auxiliary flux and requires the same maximum as its gas-only call. The job in
`build/mpi-reflux-diagnostic-scope-20261105/` checks MPI2/3, requires this new
regression to fail against the preserved old module, and reruns common/recursive
reacting restart. These results are pending. This change affects diagnostic
metadata, not reference tolerances; no comparison was relaxed.

### 2026-11-05 diagnostic correction and reacting restart verified

The reflux diagnostic regression passes on MPI2 and MPI3. Compiling the same
new regression against the preserved old reflux implementation fails at the
specific conserved-maximum assertion, establishing that the test detects the
bug. The first test attempt changed auxiliary-flux magnitude without adjusting
its analytical transfer/inventory expectations; those fixture quantities were
made consistent (x flux30, mode6 y flux21, transfer ratio3000) without changing
tolerances. Both failed attempts remain under `attempt1/` and `attempt2/`.

After checkpoint byte equality was restored, the existing corruption helper
stopped because it admitted only V1/V2 checkpoint tags. The writer/reader share
the same context/metadata/boundary-ledger prefix for V1 through V6; reaction
sources and heat release follow that prefix. The helper now accepts those six
known schemas. It still rejects unknown tags, and the actual V5 malformed,
truncated, changed-physics and corrupted-ledger guards are exercised below.

`build/mpi-reflux-diagnostic-scope-20261105/record.json` finishes as
`regression_and_restart_passed`: both fixed and selected reacting MOL checks
pass for common stepping and subcycling, including MPI2-to-MPI3, MPI-to-serial,
serial-to-MPI1 continuation, exact checkpoint bytes and rejection guards.
This is distinct from the two earlier passing nonreacting thermal-slip tests.

`build/first-mol-public-api-compat-20261105/record.json` also finishes as
`public_api_n8_mpi_exact`. Current public CVODE/tolerance/equilibrium APIs and
the generated explicit activation convention pass n8 serial/MPI3 original gates
with eight exact CSV pairs and confirmed MPI dispatch. The probe is relinked
against the diagnostic correction before a new n64 serial/MPI3 run. That
bounded50ns comparison is running under `n64-notification/`; no n64 result from
these current executables is claimed yet. Longer trajectories, spatial/time
convergence and the broader C2 release gate remain unqualified.

### 2026-11-06 current public-API n64 serial and genuine MPI3 pass

`build/first-mol-public-api-compat-20261105/n64-record.json` finishes as
`bounded_comparison_complete`. Both current serial and genuine MPI3 pass all
original field/budget gates at0/25/50ns on root64. All six CSV pairs, including
budget files, are byte-identical. Maximum absolute x-momentum difference at50ns
is9.603001674893252e-7 kg/(m2 s) at fine(64,26,109) for both executables.
Source and linked-artifact hashes were independently rechecked. The current
serial outputs also match all six preserved isolated serial files exactly;
`public_api_n64_audit.json` records this migration check. It does not reinstate
the invalid older MPI claim: this new run uses verified MPI compilation and
dispatch. Elapsed times are24590.7s serial and47811.1s MPI3, without a speedup
claim.

The source-assembly correction therefore has bounded evolved evidence with
current public chemistry APIs and real MPI, plus the separately verified
reacting restart checks. The chemistry convention remains explicitly compatible,
not the canonical mechanism or a validated physical model. The50ns prefix is
not a completed300ns trajectory or mesh/time convergence qualification.

Two follow-ups are running: `build/first-mol-compatible-spatial-20261106/` adds
root16/32 at the same25ns step and50ns endpoint, then applies the existing
volume-restriction state/increment L1/Linf analysis to16/32/64 for both solvers.
`build/first-mol-compatible-trajectory-20261106/` runs the original full checker
at root16 through300ns with25ns steps, generating a fresh unchanged pinned-PeleC
reference. Their original field failures, if any, will remain separate from
self-convergence estimates. No result from these running studies is claimed.

### 2026-11-06 compatible n16 full300ns trajectory passes

`build/first-mol-compatible-trajectory-20261106/record.json` completes with no
violations using the original full trajectory checker, a freshly generated
unchanged pinned-PeleC reference, current public-API compatible serial executable,
25ns root steps and recursive refinement. All13 snapshots through300ns pass.
Elapsed time including native/export/comparison is2478.5s. Across the whole
trajectory, maximum absolute differences are5.392805868531525e-9 kg/(m2 s)
in x-momentum,4.482705873165571e-12 kg/m3 in density,
1.4298420865088701e-5 J/m3 in total energy and2.529554876673501e-8 K in temperature.
The final-time x-momentum maximum is3.090824353402244e-9.

`verification.json` independently checks source/binary hashes and byte equality
of the first three candidate and native CSV snapshots against the separate
n16/50ns run. This is a full-duration pass at n16 for the explicit compatible
configuration, not proof of spatial/time convergence or all C2 combinations.
The n32/50ns spatial job remains running. A new n32/300ns matched trajectory job
is running in `build/first-mol-compatible-trajectory32-20261106/`, retaining
25ns steps and generating its own unchanged pinned reference. Its results are
not yet available.

### 2026-11-06 compatible fixed-step spatial audit

`build/first-mol-compatible-spatial-20261106/record.json` completes. The explicit
compatible current-source n16/n32/n64 runs all pass original field/budget gates
through50ns at fixed25ns steps. The analysis has160 solver/time/scope/field rows;
all256 positive-time L1/Linf comparisons strictly decrease from16/32 to32/64.
Input CSV hashes are independently rechecked in `verification.json`.

Decreasing norms do not establish asymptotic second-order behavior. At50ns,
PeleF increment x-momentum orders are1.2393(L1) and0.22070(Linf); z-momentum
orders are1.6317 and0.28864. The native values are effectively the same.
Density and temperature maximum orders are about0.22. Fine-level x-momentum,
z-momentum and temperature maxima remain at the coarse/fine interface in both
mesh pairs (fine i16 in the coarser comparison, i32 in the finer comparison).
The minimum positive-time observed order across the audit is0.05853. These
finite-resolution differences remain distinct from fixed-reference errors and
from a full mesh/time convergence qualification.

The n32/300ns trajectory remains running. In parallel,
`build/first-mol-compatible-temporal-20261106/` adds n16 trajectories with100ns
and50ns root steps to the completed25ns trajectory. Its analysis compares both
solvers at matched100/200/300ns times on the same composite grid, using original
reference gates and volume-weighted L1/pointwise Linf self-differences. These
new temporal results are pending; C2 remains open.

### 2026-11-06 n16 temporal audit at100/50/25ns

`build/first-mol-compatible-temporal-20261106/record.json` completes. All three
root timesteps pass the original full300ns trajectory gates. All192
solver/time/field/norm comparisons decrease as timestep is halved, with minimum
observed order0.54128. Analysis input hashes are independently verified.
At300ns, PeleF maximum-norm orders are1.09315(x-momentum),1.15079(y-momentum),
0.65699(z-momentum) and1.17414(temperature), closely matching native values.
These ratios show temporal reduction in the measured range, not uniformly
second-order or established asymptotic convergence.

A12.5ns n16/300ns run is started in
`build/first-mol-compatible-temporal8-20261106/` to assess the next halving.
Its local checker copy only adds refinement8 and adjusts the relocated script's
repository-root lookup; acceptance limits and duration are unchanged. The first
invocation stopped before calculation because the copied checker looked under
`build/cases`; logs and its original copy are preserved in `attempt1/`.
The corrected retry is running. The prepared analysis combines50/25/12.5ns at
matched100/200/300ns output times. The n32/300ns trajectory is still pending.

### 2026-11-06 n16 temporal extension and region localization

The12.5ns n16 run in `build/first-mol-compatible-temporal8-20261106/` completes
and passes all25 original snapshots through300ns. Source/executable/checker
and analysis-input hashes were rechecked. The50/25/12.5ns comparison also has
192 decreasing L1/Linf norms, with minimum observed order0.36087. At300ns,
PeleF maximum-norm orders are0.95845(x-momentum),0.68597(y-momentum),
1.42830(z-momentum),1.10494(temperature); native ratios are nearly identical.
This does not show uniform second-order asymptotics.

`regions.py` partitions the same temporal differences into disjoint interface,
x-wall and interior regions with1/2/4 local-cell-width bands. At300ns with a
two-cell interface band, global density, momentum, energy and temperature maxima
are in that band. Interface volume-L1 orders are about1.02--1.07; interior
orders for density/x-momentum/energy/temperature are about1.37, and interior
maximum orders are not uniformly second order either. Localization therefore
does not by itself prove an interface-only cause. This analysis changes no
solver, reference or acceptance limit. The n32/300ns trajectory remains running
and will guide the next intervention or qualification study.
