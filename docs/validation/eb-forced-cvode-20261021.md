# Forced CVODE substep admissibility, 2026-10-21

The tighter 32-cubed EB diagnostic stalled in a coarse corrector because
CVODE constrained every internal species-density step to be nonnegative.
A frozen spatial correction can remove an initially absent trace species
before chemistry produces it later in the interval. Its auxiliary forced
trajectory can therefore have a small negative intermediate value and a
valid nonnegative endpoint. The constraint made CVODE repeatedly reduce its
step rather than complete that interval.

`sundials_forced_reactor_mod.F90` now applies CVODE's per-substep constraints
only when all external species sources are nonnegative. Negative frozen
sources use the existing signed RHS/Jacobian path. The existing signed
composition bound (`Y >= -5e-13`), thermodynamic checks, endpoint
nonnegativity admission, machine-scale endpoint repair, analytic total mass
and transactional rollback are unchanged. No acceptance tolerance is widened.
Pinned PelePhysics does not install CVODE species constraints.

## Evidence from saved states

The previous moderate-tolerance attempt retained two Fortran spatial
evaluations before its 900-second timeout. An operator-only replay evaluated
those saved states without chemistry and additionally exported all species
rates. Its density, momentum and energy rates are exactly identical to both
original hydro and transport captures. Mean corrector forcing is negative
for initially absent water in 10,669 coarse fluid cells; the largest is
`-5.02590e-5 kg/m3/s` at `(1,9,17)`. Other initially absent radicals also
have negative frozen forcing. The original timeout remains a failed run.

Twelve saved coarse-cell inputs were reconstructed, including the mean of
both operators and the predictor energy correction recovered from stored
energies. These are isolated forced-reactor diagnostics, not a rerun of the
full spatial trajectory. Four tolerance settings were applied to the same
inputs with a 5000-step limit:

| Relative / absolute density tolerance (kg/m3) | Old accepted | Fixed accepted | Old CPU s | Fixed CPU s |
| --- | ---: | ---: | ---: | ---: |
| 1e-10 / 1e-15 | 12/12 | 12/12 | 0.0599 | 0.0469 |
| 1e-12 / 1e-16 | 12/12 | 12/12 | 0.2200 | 0.0635 |
| 1e-12 / 1e-17 | 12/12 | 12/12 | 2.5831 | 0.0834 |
| 1e-12 / 1e-18 | 4/12 | 12/12 | 9.6222 | 0.2867 |

Across all 48 completed fixed solves, the maximum absolute mass-balance
residual is `1.94e-17 kg/m3`; the largest relative H/O/N balance residual
is `1.89e-13`. Refused old solves are not included as completed balances.

## Analytic regression and coupled check

The existing forced-CVODE unit executable now includes an independent
three-species analytic problem: `A -> B`, rate 10/s, inert `C=1`, and
frozen sources `1e-10*(1,-0.1,0)`. Identical species thermodynamics keep
temperature constant. The exact solution is
`A=1e-10*(1-exp(-10*t))/10`, `B=0.9e-10*t-A`, `C=1`.
The intermediate B minimum is about `-5.18e-14`, inside the unchanged
signed-trial bound; B is positive at 0.1 s. The old backend fails the
5000-step limit. The fixed backend passes the analytic endpoint check and
the existing unit checks. Ending the same problem at 0.001 s, where B is
still inadmissibly negative, is refused with exact input/temperature rollback.

The 8-cubed root/fine coupled diagnostic also passes all original field,
species, inventory and budget gates at 0, 50 and 100 ns with the native
activation mechanism and held-EOS native comparator. Native takes 2.121 s
and PeleF 31.367 s. The largest numerical energy-budget residual is
`7.23e-17`. This is `diagnostic_passed`, not unchanged-PeleC qualification.

## Reproduction and limits

- `build/eb-forcing-replay-20261021/`: operator replay, input reconstruction,
  before/after cell outputs, elapsed CPU times and element audit.
- `build/eb-forced-signed-20261021/`: isolated corrected reactor object,
  generated fixture, exact build commands, old/new unit logs and source hashes.
- `build/eb-forced-signed-20261021/mesh8/matched_eb_amr_ssg9l2uk/checks.json`:
  completed coupled comparison and unchanged acceptance limits.
- `build/eb-forced-signed-prefix-20261021/`: completed five-evaluation 32-cubed
  PeleF capture (734.275 s), reusing the already completed native prefix at
  matched density tolerance 1e-18 kg/m3 and relative tolerance 1e-12.
  Wrapper exit 0, deliberate solver exit 86, completion sentinel and all five
  hydro/transport captures are verified. The completion notification and logs
  are in `job/`. This is an intentionally truncated diagnostic, not a full step.

## Completed prefix and remaining amplification

The water-controlled limiter factors now closely agree: at coarse parent
`(1,0,14)`, native/PeleF are 0.31372898/0.31373980, and at `(2,11,17)` they
are 0.02983908/0.02983882. The former default-tolerance differences were
0.31429522/0.08775575 and 0.02985883/0.07789892. Those comparisons change
both the chemical tolerances and the forced backend; they do not isolate
the backend's numerical effect from the tolerance effect.

The first fine corrector still differs. Its largest spatial energy-rate gap
is `1.805e6 J/m3/s` at `(21,32,24)`, while valid input conserved states agree
to roundoff. Other halo maxima are larger than in the former capture:
total energy differs by `8.922 J/m3` at `(26,27,35)`. The shared limiter is
now controlled by tiny z-momentum variations: its x factors at parent
`(15,0,15)` are 0.42500654/0.60627244. Do not call the prefix a field pass
or launch a full trajectory merely because the forced solve completed.

Saved-state caloric evaluation with the unchanged native mechanism identifies
an additional precision effect. At first coarse predictor cell `(2,5,16)`,
the captured temperature difference is `-6.75255e-7 K`; independently
converging the caloric equation for each saved state reduces it to
`2.27e-13 K`. Across all 22,528 common fluid cells the largest converged
temperature difference is `8.37e-11 K`. Both production inverses deliberately
retain the pinned `1e-6 K` Newton stopping rule, but different temperature
seeds can leave different sub-tolerance residuals. A direct caloric check
does not itself prove the later halo-error cause. The residual audit and
operator-only intervention are in `build/eb-eos-residual-20261021/`.
Its first analysis assumed identical exported cell sets and was refused;
the corrected analysis explicitly records unmatched cells and compares the
common fluid cells, as the stage comparator does.

The completed operator replay holds that entire saved coarse predictor fixed
and changes only PeleF's pinned caloric Newton stop from `1e-6` to `1e-12 K`
in a private object. The default replay reproduces both hydro and transport
rates exactly. No chemistry is integrated in either replay. Against the
unchanged native capture, maximum absolute rate differences decrease as follows:

| Rate | Default inverse | Tighter inverse |
| --- | ---: | ---: |
| Density, kg/m3/s | 1.8182e-4 | 3.2019e-6 |
| z-momentum, kg/m2/s2 | 0.10692 | 0.0022470 |
| Total energy, J/m3/s | 669.872 | 11.2858 |

Each replay takes about 25.3 s. This isolates an EOS stopping-error contribution
to the first coarse spatial rates; the later halo effect still needs a coupled
intervention. `build/eb-tight-eos-prefix-20261021/` completed that bounded
five-evaluation intervention with private physical and signed-trial EOS objects.
Its native prefix remains unchanged at the original `1e-6 K` stopping rule,
so this is explicitly a one-sided causal diagnostic. Production sources,
acceptance limits, old runtime archives and native results remain unchanged.
It completes in 425.590 s within the 900-second solver limit; the verified
completion notification is in `job/`. The resulting halo improvement and
explicit precision-mode implementation are in the
[caloric convergence record](eb-eos-convergence-20261021.md).

GNU Fortran 13.3.0, Debug bounds/FPE checks and SUNDIALS 7.2 are used. The
reactor object is compiled and linked explicitly in the new artifact
directory; old runtime archives and executables are preserved. The analytic
test uses a fixture generated with the existing mechanism generator; no
generated chemistry source is hand-edited. An initial unit build lacked the
fixture module and was retried after generation. Early analytic fixture
drafts were corrected to stay inside the existing signed-trial bound and
to use machine-scale error for the unit-density carrier; only the final
`verified-unit-*` run is the completed unit evidence.

Full 32-cubed field comparison, mesh order, full configuration matrix,
MPI/restart, installed-binary audit and physical validation are not completed
by this change. The original spatial-EOS reference gate has not been rerun.
