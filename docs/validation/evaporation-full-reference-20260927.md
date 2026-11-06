# Full evaporation reference diagnostic — 2026-09-27

Implementation: `a9cc5569`; candidate binary built from `be25365087344799b8798f8ccd0aadddd9b0bee7`.
The fixture now explicitly selects the existing `pelephysics` drag law.
Abramzon–Sirignano heat/mass transfer was already enabled. No numerical kernel
or acceptance tolerance was changed.

The pinned PeleC Spray-Evaporation executable completed two 1e-7 s steps
with two SDC iterations on the prepared periodic 4³, 400 K gas / 310 K liquid
adaptation. This is not the original upstream 273 K case.
Source: PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`; complete dependency identities
are in `build/pelec-pinned-hit-source/pelef-reference-provenance.json`.
Fresh reference build: GNU 13.3.0 Release, serial, particles ON, tests OFF,
`build/pelec-pinned-spray`; configure/build logs use that path prefix.
Candidate configuration is recorded in the
[SDC integration record](nonreacting-particle-sdc-20260927.md).

## Results

All four parcels have the same values. Changes are measured from the shared
initial state; reporting relative final temperatures would hide this discrepancy.

| Quantity | Reference change | Fortran change | Change-magnitude difference |
| --- | ---: | ---: | ---: |
| Temperature (K) | 0.004727062899 | 0.003210010595 | -32.09% |
| Single-drop mass loss (kg) | 1.133273583e-16 | 1.289211782e-16 | +13.76% |
| Axial speed reduction (m/s) | 0.000697978230 | 0.000615191466 | -11.86% |

**Particle trajectory parity fails to follow from these runs.** Matching the
drag model alone does not resolve the differences. Film transport properties
and integration remain candidates for diagnosis, not established causes.
No gas-field comparison, CTest matrix, multirank or physical validation was run.

## Reproduction and artifacts

Reference executable:
`build/pelec-pinned-spray/Exec/RegTests/Spray-Evaporation/PeleC-Spray-Evaporation`.
It was run with `reference.inp` from `build/spray-evaporation-sdc-run2`, exit 0;
log `reference.log`, output `spray00002.p3d` in that directory.

```sh
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-reference/heptane_3sp.json --output build/spray-evaporation-matched-drag
# From the prepared directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-spray-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
# From the repository root:
python3 tools/compare_evaporation_particles.py --reference build/spray-evaporation-sdc-run2/spray00002.p3d --candidate build/spray-evaporation-matched-drag/fortran.csv.particles.rank0.csv --output build/spray-evaporation-matched-drag/comparison.json
```

One reference run and one matched-drag candidate run completed with exit 0.
The diagnostic records input hashes, CGS-to-SI conversions and per-ID changes;
successful execution does not assert parity. Mechanism provenance, Python syntax,
project contract and `git diff --check` passed. Existing earlier candidate
artifacts retain the default Schiller–Naumann drag setting.

## Current-source refresh — 2026-10-29

The current dirty main at0a0e00ebe9138145fbb0e0e745665d0d768df009 was built
in a new GNU/OpenMPI Release directory with selected pinned heptane Simple,
MPI/spray ON and tests/SUNDIALS OFF. Build elapsed845.61s. The same4^3 periodic,
400K gas/310K liquid, four-parcel, two100ns-step fixture completes. This is
source-bound diagnostic evidence, not a new acceptance gate or physical test.

Artifacts: `build/evaporation-current-reference-20261029/`. `record.json` records
source/binary hashes, configuration and commands; `native-source-audit.json`
checks130 shared native-source hashes and the pinned provenance manifest.
Fresh native final particle and exported gas files are byte-identical to the
historical frozen reference, independently checked in `verify.py`. The native
run completed in0.13s and candidate MPI1 run in1.04s. Source remained unchanged
through the build/run. Existing compiler warnings are retained in build.log.

All four particles have temperature difference3.979e-13K, axial-velocity
difference-3.109e-15m/s and single-drop mass difference8.078e-28kg. The current
candidate's gas maximum differences, normalized by the maximum absolute native
field, are2.670e-6/2.670e-6/2.654e-6 for x/y/z momentum,2.360e-10 for total
energy,2.218e-9 for stored temperature and7.979e-13 for fuel density. The
comparison tool has no parity gate; successful execution is not a parity pass.
The combined candidate energy budget residual is1.897e-19J. Component budgets
are independently reconstructed from the published initial/final and external
terms, rather than accepting only the final log marker.

A separate one-step candidate and native plt00001 export localize onset.
At100ns, maximum energy difference1.235e-6J/m3 is at the parcel cell(2,2,2);
maximum x/y/z momentum differences about6.721e-11kg/m2/s are at its lower
adjacent cells(1,2,2)/(2,1,2)/(2,2,1). At200ns the energy maximum is still at
the parcel cell,2.442e-6J/m3, while momentum maxima thereabout rise to
2.844e-10/2.844e-10/2.854e-10kg/m2/s. Integrated gas-energy differences are
1.255e-15J and2.509e-15J. This is consistent with the previously documented
explicit-native versus conservative-candidate energy-feedback difference;
it does not prove the neighboring momentum discrepancy's stage-local cause.
Do not remove conservative feedback merely to match this residual. The next
discriminator is same-stage deposited energy, pressure and directional face
flux, preserving the original native baseline and candidate conservation.

MPI2 completes the same two-step case and reproduces MPI1 gas and budget CSVs
byte for byte, plus all four particle rows exactly after collecting by ID.
This is a single regular patch with an empty rank, not multi-patch particle
migration or general parallel spray qualification. `extension-commands.json`
and `verified.json` retain commands, regional measurements and hashes.
Python syntax, heptane provenance, project contract and diff whitespace checks
pass. C3 gas/particle trajectory qualification and C2 remain open.

## Measured stage particle-source discrepancy — 2026-10-29

Isolated source copies now record the candidate common particle callback's
input gas and returned conserved rates, and native old/new spray sources before
hydro and after new-source evaluation. Native post-source records also contain
hydro momentum rates. Production source files are unchanged. The instrumented
candidate reproduces ordinary gas, budget and particle CSVs byte for byte;
the instrumented native reproduces final particle/exported gas files exactly.
All compile/link/run/export commands succeed. Artifacts and exact identities:
`build/evaporation-stage-source-20261029/record.json`, `trace-schema.json`.

`analyze.py` validates384 candidate rows and512 native rows,64 unique cells per
sample, clocks and stage ordering. Candidate stages0/2 correspond to old/final
new source values; the intermediate candidate rate retains the old source
exactly. Native sub_iteration0/1 is a different numbering convention. Native
energy-rate CGS values are multiplied by0.1 and momentum rates by10 for SI.
The native phase0 hydro entries are placeholders, not sampled rates.

At the parcel cell, candidate-minus-native old/new gas energy source differences
are12.55009/12.54810J m^-3 s^-1 in step1 and12.54712/12.54507 in step2.
Trapezoidal source integration over the64-cell periodic volume predicts
cumulative energy differences1.2549094e-15 and2.5095186e-15J. Observed gas
energy differences are1.2546880e-15 and2.5093123e-15J; residuals are
-2.215e-19 and-2.064e-19J. This directly closes the global gas-energy discrepancy
from measured particle sources, rather than inferring it solely from endpoint
energy. It does not resolve the local momentum discrepancy or establish
physical validity. Source ordering agrees; the energy-source values differ.

At the same source cell the old/new axial momentum-rate differences are
1.718e-10/4.110e-11kg m^-2 s^-2 in step1 and9.303e-10/1.605e-9 in step2.
Further attribution requires matching actual pressure reconstruction and
hydrodynamic face flux around the neighboring momentum maxima. Native Sborder
at the post-source phase is the transport-stage input and must not be assumed
to be the particle-sampling state without checking its call context.
Conservative candidate feedback remains unchanged. Python compilation,
project contract and diff whitespace checks pass. C3 parity remains open.


## Predictor-only energy intervention — 2026-10-29

A diagnostic copy of the MPI EB hydro module adds the measured first-step
native-minus-candidate old particle energy rate to the assembled CTU predictor
source, after the normal halo/interior consistency checks. The same periodic
mapping is used for interior cells and halo images. The frozen correction is
applied in both predictor evaluations of the one100ns step. Actual conservative
particle exchange and production sources remain unchanged. This is a causal
control of one pathway, not an adopted numerical policy or complete native
source replacement.

The previous halo-callback intervention failed on the first step. Source
inspection identifies its conflict with the exact interior-source consistency
contract; no trajectory or causal result is inferred from that failed run.
Its original record remains at `build/evaporation-predictor-energy-control-20261029/`,
including `failure-analysis.json`.

The corrected placement is recorded at
`build/evaporation-predictor-energy-kernel-control-20261029/`. Compilation,
linking, zero-correction and nonzero-correction runs complete. Zero correction
reproduces the ordinary one-step gas, combined-budget and particle CSVs byte
for byte. `analyze.py` independently checks this identity, completion markers,
64-cell correspondence, both predictor markers, and all eight budget residual
expressions. `record.json` records commands, source/binary hashes and timings.

The transverse momentum maximum difference falls from6.721348e-11 to
5.490127e-11kg m^-2 s^-1 in both x and y: a reduction of18.318%, with81.682%
remaining. The worst cells remain(1,2,2) and(2,1,2); z behaves similarly.
The density maximum difference falls from1.88627e-13 to4.21885e-15kg m^-3.
The total-energy maximum difference instead increases by1.54%, consistent
with retaining the actual energy exchange rather than replacing it. The
combined energy-budget residual is-3.524e-19J. These observations establish
sensitivity to predictor energy forcing, not an additive decomposition of all
momentum error, reference parity, or physical validation. Pressure/face-flux
tracing remains necessary for the residual. C2 and C3 remain open.


## Directional hydro-flux accounting — 2026-10-29

`build/evaporation-face-flux-20261029/` records unmodified candidate and native
hydro momentum rates and six surrounding directional face fluxes for every
cell in each of four hydro evaluations over two100ns steps. Isolated module
copies emit the traces. Both instrumented executables reproduce their ordinary
final gas/particle outputs exactly; candidate budget CSVs are also byte-identical.
Compile/link/run/export commands and source identities are in `record.json`.

`analyze.py` checks768 cell/component records per solver, unique matched keys,
and directional divergence closure. Candidate fluxes are time-integrated
momentum flux density(Pa s); native fluxes are extensive force(dyn). For this
1mm cell, candidate face differences divide by dt*0.001; native differences
multiply by1e4 to obtain kg m^-2 s^-2. Maximum reconstruction residuals are
1.489e-8 and7.840e-9kg m^-2 s^-2, within the explicit face-roundoff bound
32ULP(maximum face magnitude), propagated through the divergence scaling.
This arithmetic bound is not a trajectory acceptance tolerance.

At the first hydro evaluation's worst x/y cells, the rate difference is
-1.230655e-4kg m^-2 s^-2 and comes entirely from the matching normal-direction
flux divergence; the other directional differences are zero. Subsequent
maxima remain normal-direction dominated. The captured final hydro-rate
difference times dt at the first-step endpoint x/y maximum is-1.118606e-11,
versus total endpoint difference-6.721348e-11kg m^-2 s^-1(16.64%). The
arithmetic remainder is-5.602742e-11. This16.64% captured-rate accounting is
a different quantity from the18.32% predictor-intervention response above.
The remainder is not yet assigned to transport: downstream limiting,
non-hydro source integration and stage combination must be traced explicitly.
Endpoint accounting verifies zero initial gas momentum from both recorded
stage inputs. No production numerical policy or qualification gate changes.


## Accepted update decomposition — 2026-10-29

`build/evaporation-update-ledger-20261029/` records each nonreacting candidate
update immediately after old+dt*rate, including accepted hydro, old/new
transport, particle and external rates. Native records after computeTemp
include the corresponding old/new state and hydro/transport/particle terms.
Both isolated builds complete, reproduce ordinary final outputs byte for byte,
and retain source/binary identities and compiler commands in `record.json`.

The768 records per solver close their update expressions: maximum state
reconstruction residual is0 for candidate and2.711e-20kg m^-2 s^-1 for native.
Candidate accepted hydro rates equal all768 pre-correction face-trace rates
exactly. Thus flux correction does not explain this fixture's discrepancy.
At the first-step x/y maximum, hydro contributes-1.118606e-11 and transport
-5.602742e-11kg m^-2 s^-1, respectively16.64% and83.36% of the total
-6.721348e-11. Particle and external contributions are zero there; remaining
arithmetic closure is-2.045e-21. This assigns the previously unattributed
first-step residual to the actual transport update, not merely an endpoint
correlation. It does not yet distinguish transport inputs from its operator.

These are pre-publication candidate updates, not universally the final CSV
state. At the second-step x/y/z worst cells, candidate publication changes
momentum by+4.163336e-14kg m^-2 s^-1; native recorded and exported values
agree there. `analysis.json` separately measures candidate post-update
changes rather than counting them as transport. Next attribution must compare
transport input stencils and viscous fluxes. No production numerical changes,
new acceptance gate, C2/C3 closure, or physical-validation claim is made.


## Viscous face capture — 2026-10-29

`build/evaporation-viscous-face-20261029/` records1440 face evaluations per
solver over six transport evaluations. Captures include left/right primitives,
temperatures, shear/bulk viscosity, all nine velocity-gradient components and
three momentum fluxes in SI. Both isolated executables reproduce ordinary
final gas/particle outputs; candidate budget CSVs also remain identical.
Each recorded flux satisfies its recorded Newtonian stress formula.

Initial uniform-state shear/bulk viscosities differ only by1.016e-20Pa s.
At the first nonzero transport evaluation's largest face-flux discrepancy,
left/right velocities are identical and both viscosity differences are
1.016e-20Pa s, whereas a tangential velocity-gradient difference is
0.0365646s^-1 and shear-flux difference is-8.48148e-7Pa. The face mapping in
`analyze.py` follows source loop order (interior axes then boundary axis/side);
normal ordering and local stress closure are checked. Full stencil-location
validation remains separate from this mapping.

Source inspection shows that absent physical-face labels leave periodic flags
false: interior cells at the patch edge use one-sided tangential derivatives,
and donor boundary faces use the interior tangential derivative. Native regular
periodic diffusion uses centered, periodically extended neighbors. This is a
candidate mechanism for the measured transport discrepancy. An isolated4^3
periodic control is being checked before any general MPI/AMR change. Global
periodic wrapping is not an acceptable fix for arbitrary patches or coarse/fine
interfaces. C2/C3 qualification remains open.

The isolated control completed at
`build/evaporation-periodic-transport-control-20261029/`. Zero mode reproduces
ordinary one-step gas/budget/particle CSVs byte for byte. Enabling periodic
transport gradients only in the regular4^3 single-patch diagnostic reduces
x/y maximum momentum difference from6.721348e-11 to1.230786e-11kg m^-2 s^-1,
an81.688% reduction; z decreases similarly. Combined energy-budget residual
is-2.439e-19J, with all eight published budget expressions independently
reconstructed. The total-energy maximum difference changes by only0.0352%.
This intervention supports patch-edge tangential stencils as the dominant
momentum-error pathway in this fixture. It neither removes the separate
conservative/native energy-source difference nor qualifies general MPI/AMR
transport. The next implementation must obtain both tangential stencils from
actual neighboring/ghost cells instead of globally wrapping a local patch.

## Full-neighborhood CTU implementation candidate — 2026-10-29

The first full-halo pilot refused its first step because the copied MOL call
requested coarse-linear halos without the required held-energy field; its
failed record remains at `build/evaporation-full-halo-transport-pilot-20261029/`.
The corrected pilot retains CTU's existing interpolation/energy convention
while using the already implemented complete transport neighborhood.
`build/evaporation-full-halo-transport-pilot2-20261029/` completes one/two steps
and MPI2. MPI1/MPI2 gas and budget CSVs match byte for byte. First-step transverse
maximum difference is1.230786e-11kg m^-2 s^-1; at200ns it is4.1427e-11,
about85.43% below the previous ordinary result. Component budgets reconstruct;
two-step combined energy residual remains1.897e-19J.

Matched serial/MPI source changes now route ordinary CTU diffusion through
complete neighborhoods while retaining prior coarse interpolation and energy
conventions. Existing MOL routing remains unchanged. The source change and
pre-change copies are at `build/ctu-transport-full-stencil-20261029/`.
Focused Debug regressions and current-source Release reproduction are pending;
this paragraph does not count them as passed or declare general adoption.

A separate one-cell x-translation check uses frozen before/after executables,
two100ns steps and the same periodic4^3 fixture. It rolls output indices back
before comparing momentum. The preset1e-12kg m^-2 s^-1 bound fails before:
maximum x/y/z differences2.424e-10/2.500e-10/2.500e-10. The corrected pilot
has exactly zero difference in all three components. This independently
exposes and removes the boundary-location dependence in the fixture. Script,
inputs, logs and binary hashes are in `ctu-transport-full-stencil-20261029/translation/`
and its parent `check_translation.py`. This is a bounded metamorphic check,
not MPI partition, AMR convergence, fixed-reference parity or physical validation.

The first current-source Debug run passes auxiliary-source coupling,
rank-neutral SDC checkpoint and spatial-halo MPI2/MPI3 tests (four passes,
no skips). The diffusion-callback comparison fails: its callback still invokes
the old narrow-stencil transport helper while the serial built-in path now
uses the full neighborhood. The failed JUnit/log remain in
`build/ctu-transport-full-stencil-20261029/`. The recheck adds independent exact
built-in serial/MPI state, temperature and reaction-history assertions before
running the callback, and supplies the full-stencil operator through that
callback; its deliberate density corruption and transactional refusal remain.
No exact comparison or tolerance was relaxed. Recheck/Release results are
pending at `build/ctu-transport-full-stencil-recheck-20261029/`.

The reusable bounded regression is
`tools/check_evaporation_periodic_translation.py`. It validates the prepared
fixture, runs both translations, rejects nonfinite momentum, records binary,
input/output hashes and commands, and uses the preset absolute momentum bound.
For example:

```sh
python3 tools/check_evaporation_periodic_translation.py \
  --executable build/evaporation-full-halo-transport-pilot2-20261029/candidate-probe \
  --case-dir build/evaporation-current-reference-20261029/candidate \
  --output build/translation-check-new-directory --ranks 2
```

The final tool independently fails the frozen original executable (exit1)
and passes the corrected frozen pilot with MPI2 (zero differences). Exact
artifacts are `translation-tool-before/` and `translation-tool-final/` under
the recheck directory. This tool result does not replace the pending
current-source Release, callback or wider AMR evidence.

## Current-source verification completed — 2026-10-29

The recheck completes the previously pending evidence. The matched serial/MPI
CTU full-neighborhood implementation passes all five focused tests, with zero
skips: diffusion callback, auxiliary coupling, SDC checkpoint, and spatial
halo MPI2/MPI3. The diffusion test independently asserts exact built-in
serial/MPI state, temperature and reaction history before its callback check.
The first failed callback attempt remains in the original JUnit record and
is not counted as a pass. `verified.json` independently reconciles both runs.

The GNU/OpenMPI selected-heptane Release rebuild completes in710.07s and
reproduces the isolated corrected pilot's two-step gas and budget CSVs byte
for byte. MPI1/MPI2 particle rows match each other and the pilot exactly by ID.
The reusable periodic-translation check also passes this newly built Release
binary with MPI2, with exactly zero momentum difference after shifting indices
back. Its executable SHA256 is
`cf290ac3c79844d2b8e71909f396c7ce075820831736076c87c49ceaeb1de3bf`.
Debug recheck build/test times are70.01/189.25s; commands, source hashes,
logs and outputs are in `build/ctu-transport-full-stencil-recheck-20261029/`.
Provenance, Python syntax, project contract and diff whitespace checks pass.

This closes implementation and bounded regression verification for the CTU
patch-edge transport-stencil correction. It does not close C3 fixed-reference
spray parity: the conservative/native particle energy-source difference and
stored-temperature discrepancy remain. The regular single-patch MPI2 run has
an empty rank and is not a particle-migration qualification. C2 EB maximum-norm
convergence and C3-C5 remain open; no physical validation is claimed.

## Combined transport/source diagnostic — 2026-10-29

`build/evaporation-combined-control-20261029/` combines the implemented
full-neighborhood transport with the earlier frozen old particle-energy-source
correction applied only inside the CTU predictor. Physical particle exchange
remains conservative. Zero correction reproduces corrected one-step gas,
budget and particle outputs exactly. The combined run completes; all eight
budget expressions reconstruct. Production sources are unchanged by this
additional diagnostic.

A four-way comparison retains original, energy-only, transport-only and
combined fields. First-step x/y maximum differences are respectively
6.721348e-11,5.490127e-11,1.230786e-11 and4.001729e-15kg m^-2 s^-1.
The combined error is5.954e-5 of the original; z is2.547e-15. Maximum algebraic
interaction(combined-energy_only-transport_only+original) is1.4552e-15.
These measured controls account for nearly all first-step momentum discrepancy
in this fixture. They do not establish exact additivity for other states or
later times. `factorial.json` pins all four field files and definitions.

Only the transport correction is implemented in production. The frozen source
substitution is a diagnostic and cannot serve as a general coupling model.
The actual accepted liquid-energy change remains the gas-feedback authority;
changing it to erase the native explicit-source difference would violate the
existing discrete conservation contract. C3 trajectory closure and longer-time
source consistency still require separate evidence. Python syntax, project
contract and whitespace checks pass; no new physical qualification is inferred.

## Matched-endpoint time study — 2026-10-29

`build/evaporation-time-convergence-20261029/` independently repeats native
and current corrected candidate at100/50/25ns steps through the same200ns
endpoint(2/4/8 steps). All six runs and native exports complete; candidate
budget expressions reconstruct at each setting. The two-step binary is the
same verified corrected Release. No source model, grid or reference is changed.
The particle comparison report's hardcoded "Two-step" scope label was removed
from `tools/compare_evaporation_particles.py`; its numeric logic and gates are
unchanged. Reports were regenerated with recorded commands.

The candidate-minus-native x-momentum maximum difference grows from
4.143e-11 to1.286e-8 and1.687e-8kg m^-2 s^-1 as dt decreases. Total-energy
maxima are2.446e-6,7.156e-6 and8.637e-6J m^-3. Thus the short two-step result
is not trajectory-convergence closure. Particle maximum temperature differences
remain3.98e-13/6.25e-13/9.66e-13K; axial-velocity differences are at most
3.33e-15m s^-1. The growing discrepancy is predominantly in the gas fields.

Separately, each solver's temporal self-differences decrease in L1 and Linf for
all nine recorded gas/species observables. Observed orders are approximately
1.00-1.09 in both solvers. This self-convergence does not imply convergence of
their mutual difference or excuse the unresolved fixed-reference trajectory.
`analysis.json` reports these two evidence classes separately. Prefix exports
are being used to locate the first onset; no causal attribution is made yet.

Twelve independent candidate prefixes and matching native plotfile exports
complete at the requested clocks. In both shorter-step runs the first large
increase occurs at step3, not at a shared physical time. With50ns steps,
x-momentum maximum differences at steps1/2/3 are1.540e-12,4.772e-12 and
2.888e-9kg m^-2 s^-1; with25ns steps they are1.909e-13,6.301e-13 and
1.851e-10. The worst cell remains(1,2,2). This identifies the third update
and its carried state/history as the next diagnostic target, without proving
a history bug. Prefix inputs/logs/exports and `prefix/analysis.json` retain
the evidence. Source/pressure/transport decomposition must be extended through
this update before changing another numerical policy.

## Temperature-policy controls after step2 — 2026-10-29

The historical native `CONVERGENCE` diagnostic binary was first anchored by
repeating its recorded two-step particle and gas outputs byte for byte, then
run for four50ns steps. Against the unchanged strict-EOS candidate, its x-momentum
maximum difference at200ns is2.072e-11, compared with1.286e-8 for the pinned
default native. The temperature maximum difference is1.878e-9K. Records:
`build/evaporation-third-step-tight-20261029/`. This is evidence that the native
inversion stopping policy strongly affects the later trajectory, not permission
to replace the pinned default authority. The earlier two-step observation that
tightening native inversion barely changed momentum does not extend to step3+.

A distinct candidate-only control enables the existing
`pinned_eos_extrapolation=.true.` option, which selects the pinned default
inversion behavior as well as its polynomial policy. No solver code or stopping
tolerance is changed. At100/50/25ns dt and200ns endpoint, x-momentum maximum
differences against the unchanged default native become1.000e-9,3.610e-10 and
9.922e-11kg m^-2 s^-1; temperature differences are8.806e-8,2.098e-8 and
4.761e-9K. This improves the refined-step comparisons but worsens the two-step
momentum comparison, so neither configuration is silently substituted for
the other. All three candidate runs complete and budget expressions reconstruct.
Inputs and `verified.json` at `build/evaporation-pinned-eos-time-20261029/`
keep the controls separate. Four-step stage-ledger capture remains in progress
at `build/evaporation-third-step-ledger-20261029/`.

## Third-step ledger verification — 2026-10-29

The four50ns-step instrumentation completes and reproduces ordinary final
candidate gas/budget/particle and native gas/particle outputs byte for byte.
All1536 cell/component samples per solver close the recorded update:
maximum residual0 for candidate and1.3553e-20kg m^-2 s^-1 for native.
Candidate final published momenta equal the last captured update exactly.
At step3, cell(1,2,2), x-momentum difference-2.887959e-9 consists of prior
-4.772254e-12, new hydro-2.883208e-9, transport+2.135455e-14 and zero local
particle/external differences; remaining arithmetic residual is-1.72e-21.
The y component behaves similarly. Thus the short-step increase is in hydro,
not a reappearance of the corrected transport-stencil error. Combined with
the independent inversion controls, this narrows the next pressure/primitive
comparison; it does not establish a source-history implementation defect.
Artifacts: `build/evaporation-third-step-ledger-20261029/analysis.json`, source
copies, exact commands and hashes in `record.json`.

`tools/prepare_evaporation_reference.py` now exposes
`--eos-policy strict|pinned`, recording `candidate_eos_policy` in its metadata.
The default strict input is byte-identical to the prior fixture. Preparing
both modes verifies that only the candidate EOS line differs in input files;
native gas/particle inputs, spray and geometry are identical. This makes
configuration selection explicit rather than silently changing the reference.
Preparation artifacts are `build/evaporation-eos-preparation-{strict,pinned}-20261029/`.
Python compilation, provenance checks embedded in preparation, project contract
and whitespace checks pass. Neither policy is declared a complete C3 gate.

## Third-step face divergence — 2026-10-29

`build/evaporation-third-step-faces-20261029/` extends the output-preserving
face capture to four50ns steps. Both1536-row momentum captures reconstruct
hydro divergence within propagated face roundoff; maximum residuals are
1.656e-8(candidate) and8.178e-9(native)kg m^-2 s^-2. At the third accepted
hydro evaluation, cell(1,2,2), x-rate difference is-0.0576641651kg m^-2 s^-2.
The x-face divergence contributes-0.0576641648, versus approximately
-3.91e-11 each from y and z. This agrees with the accepted-update ledger and
localizes the large increment to the normal momentum-flux pair. A separate
pre-predictor primitive capture records both temperature seeds, recovered
values, pressure and conserved inputs; its result must be checked before
assigning the entire face difference to the EOS stopping rule.

## Actual third-step primitive recovery — 2026-10-29

`build/evaporation-third-step-primitives-20261029/` completes output-preserving
pre-predictor primitive capture:512 physical-cell records per solver across
eight hydro evaluations. Mapping uses physical cell centers and validates all
64 cells per evaluation. At the third-step pressure maximum, cell(4,2,2),
candidate seed/recovered temperature is399.9999998463949K; native seed and
recovered temperature both remain400K. Pressure differs by-3.8910133e-5Pa.
The ideal-gas algebraic temperature contribution is-3.8910085e-5Pa, while
the rho*R contribution is-5.6843e-11Pa. At step4 the corresponding pressure
difference grows to-1.3612329e-4Pa. This directly locates a sub-microkelvin
stored/recovered temperature discrepancy in the hydro input adjacent to the
periodic momentum-error maximum, consistent with the separate inversion
interventions. This split compares each solver's own state; it is not a
same-input EOS replay or a proof that every face error has one cause.
The fixed default reference, production EOS options and acceptance gates
remain unchanged. Python syntax, project contract and whitespace checks pass.

## Same-input third-step EOS replay — 2026-10-29

`build/evaporation-third-step-eos-replay-20261029/replay.py` rebuilds the
standalone pinned heptane inversion harness in default and `CONVERGENCE`
modes. For each of512 captured primitive inputs, it crosses native/candidate
conserved state with native/candidate stored temperature seed (2048 inputs
per mode). Native inputs retain CGS arithmetic; candidate inputs are converted
from SI. Native kinetic-energy reconstruction follows `pc_ctoprim`. All4096
inversions return success. The512 native-state/native-seed default replays
reproduce the captured native temperatures exactly.

At evaluation5, cell(4,2,2), native state and400K seed require a first Newton
correction of-1.53628765e-7K. The default1e-6K stopping rule exits before
applying that correction, returning400K. With the same conserved state and
candidate seed399.9999998463949K, it instead returns that seed. Tight replay
returns399.9999998463712/399.9999998463713K for the two seeds. Changing only
the conserved state to candidate changes the tight result by about1.2e-12K
at this cell. Thus this local primitive temperature discrepancy is dominated
by seed retention, rather than the conserved-input discrepancy.

Across the complete crossed input set, default results have zero state effect
at matched seed, while maximum seed effect is5.41613e-7K. Tight replay reduces
maximum seed effect to1.02318e-12K; its maximum state effect is3.67936e-9K.
These are replay effects on captured inputs, not predictions of an intervened
full trajectory. Neither the pinned baseline nor production EOS defaults,
conservative exchange or acceptance tolerances were changed. Full C3 reference
agreement and C2 convergence remain open. `analysis.json` records commands,
compiler, HEAD and hashes of harness, mechanism header, inputs and binaries.
Mechanism provenance, Python syntax, project contract and whitespace checks
pass; no application test matrix was rerun for this diagnostic-only addition.

## Old-temperature seed ownership intervention — 2026-10-29

Pinned `Advance.cpp::construct_Snew` copies all of `S_old`, including UTEMP,
before each nonreacting source update. Candidate total-energy CTU recovery
instead retained the temperature of the preceding SDC iterate. With the
pinned inverse's pre-update stopping rule this changes the returned primitive
state even when conserved inputs are nearly identical.

An isolated MPI module resets the seed to old temperature before predictor
and nonreacting corrector recovery. Existing pinned EOS is explicitly enabled;
all other input bytes and the default native reference are retained. At200ns,
x-momentum maximum differences for100/50/25ns steps fall from
1.00034e-9/3.61034e-10/9.92157e-11 to
4.71372e-11/2.32803e-11/1.15889e-11kg m^-2 s^-1. Temperature maxima become
3.76144e-9/1.88214e-9/9.41213e-10K. Energy maxima are
2.43681e-6/1.21922e-6/6.09709e-7J m^-3. These decreasing errors do not
constitute a complete trajectory or mesh-convergence qualification.

All four application runs complete their prescribed endpoints. At50ns steps,
MPI1/MPI2 gas and budget CSVs are byte-identical and all four particles match
exactly by ID. Eight budget rows per run independently reconstruct their
recorded residuals within eight ULPs of term scale. Evidence, hashes and the
verification script are in `build/evaporation-old-seed-pilot-20261029/`.

The scoped serial/MPI implementation now resets only active-cell seeds for
pinned-EOS, nonreacting CTU without auxiliary internal energy, at predictor
and corrector recovery. Strict EOS, MOL and reacting/held paths retain their
existing behavior. A new `seed_aux` test exercises nonreacting pinned serial/
MPI source coupling, distributed and idle-rank layouts, conservation and
late-failure rollback. It does not independently establish native parity.
Debug tests and Release reproduction complete under
`build/ctu-old-seed-integration-20261029/`: three tests pass (seed auxiliary,
ordinary auxiliary and diffusion), zero failures or skips. The rebuilt Release
reproduces all four isolated cases, including gas, budgets and particle CSVs.
Source hashes remain stable; `verified.json` records independent XML and
output-identity checks. Release SHA256 is
`6dadb4f8203d81ad13119813fe5cb3b8cf0fd655467cae73537e62aa6b301489`. No tolerance, mechanism or conservative exchange change
is included.

The evaporation diagnostic checker now accepts an explicit positive finite
`--expected-time`, retaining200ns by default and the same endpoint check.
Default report bytes are unchanged; NaN, infinity, zero and negative values
are rejected. The three-timestep extension to1us completes in
`build/evaporation-long-trajectory-20261029/`: candidate and unchanged native
finish10/20/40 steps at100/50/25ns. Native exported clocks and candidate
completion markers agree at1us; all24 candidate budget rows reconstruct.
At this endpoint, maximum x-momentum differences are
8.88852e-10/4.40327e-10/2.19237e-10kg m^-2 s^-1; energy differences are
1.11626e-5/5.58434e-6/2.79258e-6J m^-3; temperature differences are
1.70458e-8/8.52810e-9/4.26485e-9K. Particle temperature differences stay
below6.60e-12K. Energy-budget residual magnitudes stay below5.15e-19J.
Both solvers have decreasing L1 and Linf temporal self-differences for all
nine recorded gas/species quantities on this fixed4-cubed mesh. Native
self-differences retain CGS units in the JSON; intersolver reported differences
are SI. This extension supplies bounded temporal evidence, not spatial
convergence, full C3 acceptance, or physical validation. Source-energy
feedback differences remain a separate unresolved contribution. `verified.json`
records endpoint, budget and self-difference checks; no reference tolerance
was introduced or relaxed.

## Independent1us combined-energy accounting — 2026-10-29

`build/evaporation-long-energy-audit-20261029/audit.py` freshly exports the
three native initial plotfiles and reconstructs gas plus represented parcel
energy from initial/final CSV and particle data. Gas volume is1e-9m3 per cell.
Liquid specific energy is `h_NC7H16(300K)-363000+2248.3*(T-300)`J/kg, plus
kinetic energy; Cantera3.2.0 evaluates the pinned YAML reference enthalpy.
Initial native gas energy and candidate budget initial energy agree within
32ULPs; injected parcel energy agrees exactly. All input hashes are recorded.

For100/50/25ns steps, native combined energy changes are
-1.2533653e-14/-6.2621620e-15/-3.1244372e-15J. Candidate-minus-native final
combined-energy differences are
1.2533486e-14/6.2620536e-15/3.1242912e-15J. Their sums leave only
-1.67e-19/-1.08e-19/-1.46e-19J. Independently reconstructed candidate
energy changes have magnitude below2.44e-19J; its separately accumulated
budget residuals remain below5.15e-19J. Different summation paths explain
roundoff-scale differences between these two candidate diagnostics.

The pinned `Drag.H` deposits old-state drag work and convective heating,
then mass-related kinetic energy and vapor enthalpy. Candidate
`advance_mixture_evaporation` returns the negative actual accepted parcel
energy change, including its finite Euler temperature update. The measured
endpoint accounting confirms a native combined-energy drift decreasing
approximately linearly with timestep. It does not prove every local momentum
or energy error originates there. No production source-term change is made:
matching this finite-step global drift by breaking candidate conservation
would violate the completion contract. Fixed-reference error, conservation,
and physical validation remain separate evidence classes.

## Resolution-configurable reference fixture — 2026-10-29

Preparation and comparison now accept explicit4/8/16 cells per side in the
same4mm cube. Default five solver-input files and comparison JSON remain
byte-identical; an8-cubed request checks729 level-set nodes and unchanged
parcel/spray input bytes. An8-cubed checker rejects4-cubed data. Preparation
still validates the pinned YAML and regenerated Simple transport coefficients.

`build/evaporation-mesh8-20261029/` completes candidate/native8-cubed runs
through200ns with25/12.5ns steps. All four endpoints and sixteen candidate
budget rows are checked. Maximum x-momentum differences are
1.72114e-10/8.57262e-11kg m^-2 s^-1; energy differences are
4.75385e-6/2.37874e-6J m^-3. Particle temperature differences are below
3.98e-13K. Both fuel-density maxima map to cell(4,4,4), consistent with
matching source placement; this endpoint observation is not stage-local
ownership capture. The two timesteps provide one self-difference, not an
independently measured convergence order. At matched25ns, maximum errors
are larger than on4-cubed; fixed-position point deposition changes its local
volume normalization, so this is not a spatial-convergence pass.

The16-cubed12.5/6.25ns attempt under
`build/evaporation-mesh16-20261029/` is partially completed, not passed as a
study. Candidate and native complete16 steps at12.5ns through200ns. The
candidate then rejects step7 at37.5ns with6.25ns steps and zero retries;
no native6.25ns run was reached by that sequential driver. The original
failed record is retained. `verify_partial.py` independently checks only the
completed12.5ns case:4096 native cells at200ns, candidate completion marker,
and all eight budget rows. Maximum differences are1.19125e-9kg m^-2 s^-1
in x momentum,1.82121e-5J m^-3 in energy, and2.78488e-8K in temperature.
At matched12.5ns these maxima grow13.90/7.66/7.61 times from8-cubed,
respectively; no spatial-convergence pass is established.

An isolated gas-SDC and particle-wrapper refusal trace is running under
`build/evaporation-mesh16-refusal-20261029/`, together with the previously
unreached native32-step case. This retains the same inputs, zero-retry policy,
and production guards. Floating underflow/denormal flags in the failed log
are not by themselves a diagnosis of the refusal. No numerical
solver or tolerance change accompanies this fixture/checker extension.

###16-cubed refusal localization

The first isolated trace reproduces step7 rejection at37.5ns. The gas SDC
caller receives failure from `compute_sdc_mpi_eb_tree_hydro_3d` (production
line2901), and the particle wrapper propagates `gas_accepted=false` without
publishing the candidate. This localizes the refusal to hydro computation,
not the later combined inventory test. It does not yet identify the failing
face or primitive guard. The previously unrun default-native6.25ns case now
completes32 steps; all4096 exported cells have endpoint200ns, and final
particle output is present. `build/evaporation-mesh16-refusal-20261029/verified.json`
records source stability, endpoint checks and the reproduced trace.
A second isolated trace now targets hydro reconstruction and face correction
in `build/evaporation-mesh16-hydro-refusal-20261029/`. No production guard
has been changed.

###PLM gradient-product underflow identified and corrected

The deeper trace reaches normal PLM reconstruction before transverse work:
`reconstruct_pelec_plm_cell` rejects a face produced from positive conserved
trace species. Captured sound speed400.6283m/s and dt/dx2.5e-5s/m give
CFL0.0100157. Density and pressure are uniform and positive. For fuel mass
fraction, the five-point stencil is approximately
`[1.96699e-139,2.13816e-183,2.63558e-210,4.24605e-166,2.90721e-132]`.
The central opposite-sign gradients multiply to negative zero. The former
`dl*dr>=0` condition then admits slope4.27632e-183, generating a negative
open-face weight-2.13816e-183. Exact input and same-input passive replay are
in `build/evaporation-mesh16-plm-input-20261029/`.

`compute_pelec_plm_slopes` now compares signs directly at its central and
neighbor slope checks, avoiding underflow/overflow in the sign test. It does
not clip species, change a tolerance, or permit a negative face. The captured
minimum then has exactly zero slope and retains its positive center weight.
A standalone extracted-function negative control fails before and passes
after the change. The production regression adds this captured minimum,
nonnegative face reconstruction, and a1e-200-scale linear same-sign stencil.
This corrects a finite-precision limiter defect; native `PLM.H` also contains
the product predicate, so identical behavior for every extreme trace value
is not claimed. Conservation and nonnegative admission remain authoritative.
The two focused Debug CTests pass (PLM/PPM halo and PCM/PLM viscosity),
with zero failures or skips. The rebuilt Release completes the formerly
failing32-step16-cubed case at6.25ns on MPI1 and MPI2, with byte-identical
gas and budget CSVs and identical four-parcel rows by ID. The16-step12.5ns
case also completes. All24 budget rows independently reconstruct; source
hashes remain stable. Evidence is in
`build/plm-sign-underflow-fix-20261029/verified.json`. Release SHA256 is
`5c8a86a533c4da7d3f7baa8937a601ce850e6022a3296bc18946b90400b5f364`.

At6.25ns steps, fixed-native maximum differences are5.93122e-10kg m^-2 s^-1
for x momentum and9.09456e-6J m^-3 for energy, approximately half the
12.5ns results. Particle-temperature differences remain below1.14e-13K;
combined candidate energy residual is2.90024e-18J. This proves the bounded
failure recovery and MPI equivalence, not full reference qualification.
The first step6 MPI2 checkpoint attempt in
`build/plm-sign-underflow-restart-20261029/` fails before stepping with
`Spray context capacity exceeded`; it is not a restart pass. At16-cubed,
serialized17^3 nodal values exceed the old130000-character assembly bound.
The same application already admits16000000 characters for long profiles,
uses dynamically allocated context strings, and the MPI particle reader
accepts16008192 characters including generation overhead. The assembly
bound is now aligned to16000000 without changing serialized bytes or numeric
guards. This distinct application-capacity repair and a fresh MPI2-to-MPI3
attempt are recorded in `build/eb-spray-context-capacity-20261029/`.
The repaired application writes a150347-character particle context line and
completes MPI2 step6 checkpoint to MPI3 step32 continuation. Final gas and
budget CSVs are byte-identical to the uninterrupted MPI1 run; all four
particles match by ID. Changing one positive level-set node from1 to2
(still all-fluid geometry) rejects restart with `EB spray context mismatch`
before CSV publication, confirming that the enlarged context remains fully
compared. Original failed records are retained. Source/checkpoint hashes,
negative-refusal evidence and exact identities are in `verified.json`.
The rebuilt application SHA256 is
`62c98ffe29d174bb208a4326eff8dbeebc2babc3e4ad3f99c386e73905c908ad`.
Python syntax, project contract and whitespace checks pass. This is bounded
restart and refusal evidence, not a complete release matrix.
