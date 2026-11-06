# Recorded inflection-screen admission audit

Scope: isolated serial 16-cell/50-ns paired-policy candidate, not production
integration or fixed-PeleC qualification. HEAD is
`0a0e00ebe9138145fbb0e0e745665d0d768df009` with existing uncommitted work.
Frozen source, GNU Fortran Debug compile/link commands and binary identities
are in `build/eb-inflection-admission-20261023/build-record.json`.

The completed diagnostic receipt reports 50.570 seconds, exit zero, normalized
budget residual 5.83314e-14, and six CSV tables byte-identical to the preceding
inflection evolution candidate. This turn independently inspected that saved
receipt and replayed its log; it did not rerun the solver.

`tools/analyze_eb_quartic_admission.py` replays all 64 recorded 5-cubed stencils
at parents x=16, y=1..8, z=6 or 8 over four fills. All admission decisions
match runtime exactly. The first three fills each accept 8/16 sampled parents;
the fourth accepts 4/16. These are sampled counts, not total halo adoptions.

Initially, all eight rejections first occur in component 5 (total energy),
y direction. For parent y=2 its five values are
145435.926368411, 145435.926368411, 145440.993023961,
145448.158356957, 145453.225012507. The first difference is zero while the
second differences change sign. The strict-monotonicity screen rejects this
resolved smooth energy extremum. Thus the remaining linear band at fine y=3/4
is not explained by an EOS family rejection or cut-cell geometry alone.

In the fourth fill, previously accepted sampled parents y=1 and 8 switch to
rejection. Three first fail y momentum along x: second differences reach
3.23353e-15 against noise 2.61226e-15. The fourth (16,8,8) first fails x momentum
along x: maximum second difference 3.10862e-15 against noise 2.61226e-15.
These are near the floating-point allowance, unlike the initial energy
rejections. They establish sampled admission switching, not the complete
cause of the observed y-momentum L2 increase. The first failure can also hide
later failures; simply increasing the allowance is not justified by this audit.

Reproduce without a solver run:

```sh
python3 tools/analyze_eb_quartic_admission.py \
  --input build/eb-inflection-admission-20261023/grid16.log \
  --output build/eb-inflection-admission-20261023/admission-audit.json
```

The output binds the input SHA256 and preserves every decision, first failure,
per-fill count and transition. Independent checks reproduced the earlier
64-record replay exactly; truncated input and a flipped runtime decision were
both rejected. Python syntax, mechanism provenance, project contract and
`git diff --check` passed. No numerical thresholds or production code changed.

Next: assess smooth-extremum preservation separately from numerical-noise
sensitivity using saved stencils, including discontinuity controls, before
another evolution experiment. The candidate remains unpromoted; C2-C5 remain
open. No diagnostic job remains active.

## All-rejection extension and rejected extremum rule

The analyzer now accepts `--all-failures` and records every rejecting line,
component and axis. Replaying the same 64 blocks preserves all original
first failures and decisions exactly. In the initial fill, each rejected
sampled family has 25 energy/y lines and no other failing component. In fill
four, energy failures persist at y=2,3,6,7 while x/y momentum failures appear
along x. At (16,8,8), the first x-momentum failure hides a further y-momentum
failure: fixing only the first reported component would not restore admission.
Full output: `build/eb-extremum-screen-20261023/all-failures.json`.

An offline alternative removed strict monotonicity and required quartic
children inside adjacent-parent bounds plus max|d3| <= max|d2| + noise,
in addition to the existing curvature rule. It accepts all sine/sine-squared
parents at 16/32/64 (none at 8), including all eight initially rejected
recorded energy lines. However, it also accepts the centered cell-average
step [0,0,0.5,1,1], and one of nine offsets for each tanh width 0.125, 0.25,
0.5 and 1 cell. This fails the intended discontinuity-screening criterion;
it is not evidence that those bounded child values alone cause instability.
The rule is rejected without a Fortran change or evolution run.

`build/eb-extremum-screen-20261023/evaluate.py` and `evaluation.json` retain
the complete scalar study. The local pinned quartic header still supplies
arithmetic, not an EB smooth-extremum admission rule. Next work must address
both the energy-extremum criterion and the momentum-noise sensitivity;
a fitted threshold based only on this sine case is not a justified correction.
