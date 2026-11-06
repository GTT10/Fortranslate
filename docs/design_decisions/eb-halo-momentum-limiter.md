# Separate momentum limiting for EB coarse halos

Status: opt-in application policies; historical validation stages are recorded below.

The default `interpolate_eb_linear_halo_3d` behavior remains the pinned
shared-factor interpolation. The optional `separate_momentum=.true.` path
addresses a measured discontinuity near exact rest: a roundoff-sized momentum
extremum can set the common factor to zero and flatten resolved density,
species and energy slopes. The same behavior exists in the pinned AMReX
algorithm; this is a policy candidate, not a claim of a translation error.
See [the stationary-control evidence](../validation/eb-stationary-control-20261021.md).

The candidate retains the existing centered/monotonized slope construction.
Density, species, energies and temperature share one directional factor.
Each of the three momentum components instead uses its own limited slope.
Unlike the causal ablation, it does not remove the momentum slope bound.
It introduces no absolute threshold for deciding whether momentum is noise.

Common scalar scaling preserves linear identities such as species-density
sum equals total density when the stencil satisfies those identities.
Symmetric child offsets retain the parent average of every component.
This does not prove positivity of internal energy after interpolation:
independently changing momentum and total energy still requires the existing
EOS admission and transactional failure handling in the caller.

The candidate requires the conserved ordering `rho, momentum(3), energy, ...`;
vectors shorter than five entries are refused without changing the caller.
The default remains usable for generic packed vectors. Nonfinite stencils
and invalid offsets continue to be rejected before publishing a candidate.

Focused tests cover unchanged explicit/default legacy results, resolved scalar
gradients in the presence of tiny transverse momentum extrema, limited
nonaffine momentum, nonnegative species and their sum at every child,
conservative child averages, and refusal/rollback. The implementation and
test identities are in `build/eb-separated-momentum-20261021/`.
The Debug transfer CTest passed 1/1 with zero failures or skips; build time
was 91.02 s and the test took 0.07 s. Mechanism provenance, project contract
and whitespace checks also passed. This checks the kernel, not application
dispatch or MPI. The isolated default-trajectory comparison now passes:
all six field/budget files at 0, 50 and 100 ns on the saved smooth 8-grid
control are byte-identical. The probe run took 40.74 s; exact build/link
commands are in `build/eb-separated-momentum-20261021/default-parity/record.json`.
The matched bounded-candidate 32-grid trajectory through 100 ns passes all
unchanged field and budget limits. At 100 ns the maximum density difference
is `1.2640e-11 kg/m3`, energy difference `4.5171e-5 J/m3`, and temperature
difference `1.1258e-5 K`; the relative energy-budget residual is `5.97e-17`.
Exact copied-source identities, commands and reports are in
`build/eb-bounded-momentum-trajectory-20261021/record.json`. Both solvers use
the bounded candidate in this diagnostic; it is not original-PeleC parity.
Their first
diagnostic builds incorrectly reused a preprocessed-source compiler flag for
raw source; corrected raw-source compilation succeeds and the failed logs
are retained.

The serial and MPI prepared-halo APIs now accept the same optional selector.
Selecting it requires the prepared coarse geometry used by held MOL;
other halo paths are unchanged. MPI normalizes omitted/false and collectively
compares the resolved selector before publishing output. Focused serial/MPI
tests passed 4/4 in `build/eb-momentum-halo-api-20261021/`. They cover
two owner layouts, two/three ranks, and collective disagreement rollback.
A stale serial source-halo wall expectation was corrected to the existing
mirrored-donor implementation; no source-halo numerical code changed.

Only four spatial entry points currently request this prepared linear
interpolation: serial/MPI MOL hydro and serial/MPI prepared SDC transport.
These four entry points now forward the optional selector and reject its
use outside their prepared-halo path. MPI agrees the selector before dispatch.
Operator regression checks passed 12/12 (GNU 13.3/OpenMPI 4.1.6),
including both policies at held coarse-time hydro halos and mismatched-policy
refusal. Build took 147.75 s and tests 129.71 s; no failures or skips.
Exact identities and commands are in `build/eb-momentum-spatial-api-20261021/record.json`.
Candidate transport success subsequently passed 3/3 MPI2/MPI3/constant-transport
CTests (186.33 s), covering three owner layouts, all/selected levels, both
prepared policies and mismatched/unsupported selection refusal. Exact commands
and identities: `build/eb-momentum-transport-parity-20261021/record.json`. Provenance, project contract and whitespace checks passed. Application integration must thread it through
common/recursive SDC calls and their recovery contexts. A process-global
switch would obscure these dependencies and is not used. Restart identity
must distinguish the policy while retaining the old default identity.

Before any application adoption, compare the bounded candidate against a
separately identified native diagnostic using the same policy, retain the
original-reference results, and check successful/refused EOS transitions,
serial/MPI parity, regrid and changed-rank restart. If exposed in application
inputs, the choice must be collectively agreed and bound to checkpoint
identity. Focused kernel tests alone do not satisfy these gates.

The native source-prefix record now carries a normalized halo-policy identity.
Its validator rejects a different policy before reusing packets; omission means
the original policy. Integrator producers still use the original default.
Serial and one-owner MPI mismatch regressions passed as part of 3/3
source-prefix capture, consumer and recursive retry CTests. No callbacks or
physical state publication occur on policy mismatch. Build/test times were
155.65/146.74 s; provenance, project contract and whitespace checks passed.
Exact identities and commands: `build/eb-momentum-prefix-identity-20261021/record.json`.
The selector is not yet threaded into time integration.

The common serial/MPI SDC entry points now accept the selector for native MOL
with independent internal energy and redistribution `none`. They forward it
into held hydro and linear transport, and capture/compare it in native source
prefixes and full interval histories. MPI agrees the flag before numerical work.
Prepared redistribution is explicitly rejected with the candidate policy.
Common-integrator checks passed 4/4 (189.59 s), including candidate capture
and source replay, policy mismatch refusal and the original recursive retry
regression. Provenance, project contract and whitespace checks passed. Exact
identities and commands: `build/eb-momentum-common-sdc-20261021/record.json`.
Recursive dispatch now forwards the selector through common calls and every
internal retry/recursive call; native replay windows require matching policy.
Candidate recursive retry/parity, history identity and rank-disagreement
checks passed as part of 5/5 CTests (303.29 s). Provenance, contract and
whitespace checks passed; exact identities and commands are in
`build/eb-momentum-recursive-sdc-20261021/record.json`. Full-window replay checks with late policy mismatch passed 2/2
(221.25 s): seven-clock exact gas/temperature/held-energy/source/face/budget
replay and serial/idle-rank late policy refusal. Provenance, contract and
whitespace checks passed; commands and identities are in
`build/eb-momentum-window-replay-20261021/record.json`.

The gas EB-tree application now exposes `separate_momentum=.true.` (default
false), requires held MOL/SDC with redistribution `none`, and rejects spray
or enabled LES. MPI agrees the input before initialization. The shared
checkpoint envelope emits `PELEF_EB_TREE_SEPARATE_MOMENTUM_V1` only when
selected and refuses either direction of policy mismatch before reading state.
This retains the original default envelope. Selected reacting application regrid/rebalance and changed-rank restart
checks pass 3/3, including default regression and candidate common/recursive
stepping. See [application evidence](../validation/eb-momentum-application-20261021.md).
Fixed-mechanism default and common/recursive candidate lifecycle checks also
pass, including unsupported input admission and late restart rollback; see the
same evidence record for the corrected legacy-fixture failure and retained logs.

## Density-coupled velocity option

`consistent_velocity=true` is an additional explicit policy requiring
`separate_momentum=true`. Independent momentum limiting can destroy a constant
velocity when a scalar limiter clips density. The added option limits
momentum relative to the donor velocity and adds the limited density slope
multiplied by that velocity. It retains symmetric-child momentum averages.
The previous two policies remain unchanged when this option is false.

The option is now forwarded through ordinary serial/MPI held-MOL prepared
hydro/transport, common/recursive stages and saved source/history identity.
Application checkpoints carry a separate opt-in marker and refuse changing
it across restart. It inherits the held-MOL, no-redistribution/no-LES/no-spray
boundary. The diagnostic pressure representation is a separate experiment;
this option does not enable it. See the
[application evidence](../validation/eb-velocity-application-20261022.md).

## Resolved-flow qualification boundary

A same-state 8-cubed Taylor–Green stencil audit confirms that momentum
separation changes scalar factors in 48/384 directions; 24 lift a zero
shared factor to a nonzero scalar factor. One affected packed temperature
changes by 8.17 K before evolution. The independent formula replay agrees
with the production kernel to 1.12e-16. Thus the policy applies to resolved
velocity extrema as well as the rest-state roundoff motivation. Initial
coarse/fine flux mismatch decreases, but reacting pinned-reference field
agreement worsens. Keep the default unchanged and require independent
variable-velocity accuracy evidence for the opt-in policy. See the
[recorded-stencil audit](../validation/eb-momentum-attribution-20261022.md).

The completed analytic steady-shear comparison now supplies independent
8/16/32 evidence for all four policies. Momentum separation improves most
integrated errors, while transverse maximum-norm errors increase with
refinement and z-momentum L2 is worse than default. Consistent velocity and
paired pressure add little in this constant-density shear. Keep the policy
opt-in: the analytic evidence supports a bounded benefit, not uniform
all-field convergence. See [the complete results](../validation/eb-shear-analytic-20261022.md).

## Connected-direction halo kernel candidate

The planar analytic-shear diagnosis traced the first pressure drop to the
continuity update at a coarse/fine face. Near the embedded wall, the pinned
all-27-connectivity requirement disables even tangential interpolation.
An isolated per-axis candidate improved all nonzero cut-field norms at 50 ns
with momentum separation, but worsened default y-momentum L2 by 5.47%; see
[the shear evidence](../validation/eb-shear-analytic-20261022.md).

The planner now has an optional `connected_directions` kernel selector,
default false. Omission and explicit false retain the original all-connected
guard. True permits each axis only if both immediate neighbors are connected;
a blocked axis requests the donor twice, contributing zero slope without
reading solid state. The retained subspace must also have connected diagonal
entries. A covered center, no usable axis, or disconnected retained diagonal
falls back to PCM before donor reads. Existing unavailable-neighbor and inlet fallback rules
remain. Tangential interpolation uses the existing limiter and symmetric
child offsets. This changes the numerical policy and is not PeleC parity.

Application integration now exposes this selector only with separate momentum,
while preserving false defaults. Serial/MPI common and recursive stepping,
collective agreement, saved-history policy identity and checkpoint admission
carry the choice. All nine focused integration checks pass, including fixed/selected
common/recursive lifecycle and policy refusal; see the
[application record](../validation/eb-connected-application-20261022.md).

The completed frozen-candidate 8/16/32 comparison improves integrated errors
with momentum separation but leaves local transverse convergence open.
Sphere/corner donor-enumeration tests cover topology admission only; general
curved/corner flow accuracy, positivity and rollback qualification remain open.
Symmetric-child arithmetic-average preservation alone does not establish
volume-weighted conservation in arbitrary EB geometry.
