# Vector-scaled momentum screen evolution — 2026-10-23

The optional momentum roundoff scale retains all 960 initially adopted halos
through every recorded fine-stage fill, removing the previously observed
960-to-464 switch. The 16-cell/50-ns cut-shear run completes in 144.321 seconds.
All reported nonzero error norms decrease relative to the saved connected
production baseline; default outputs remain unchanged. **This is a bounded
serial candidate result, not production adoption, convergence or C2 closure.**

## What changed and what was verified

Only the isolated smooth-screen cancellation scale for the momentum components
changes from the component magnitude to the momentum-vector magnitude. The
16-epsilon factor, curvature-sign rule, non-momentum predicates, EOS admission,
physical child-volume weighting and common-family contraction are unchanged.
The flag defaults false in the helper; the frozen candidate enables it.
The [admission audit](eb-quartic-admission-20261023.md) explains the causal
evidence and the 84 focused assertions, including preserved jump rejection.

Artifacts: `build/eb-quartic-roundoff-20261023/`. `build-record.json` and
`record.json` retain source/executable/input identities, commands, compiler
options and timings. `compare_roundoff.py` compares both the saved connected
baseline and the earlier scalar-scale quartic candidate at matched time and
component. `roundoff-comparison.json`, `comparison.json` and
`stage-comparison.json` contain the complete comparisons.

The four default ordinary tables and default spatial trace are unchanged.
All policies pass directional face-rate closure, 60 pre/post synchronization
probe pairs and two-stage fine momentum trajectory reconstruction. The maximum
normalized conservation residual across the policies is 5.93e-14. These tests
do not establish general stability, MPI correctness or restart compatibility.

## Paired-policy results at 50 ns

| Quantity | Candidate/baseline L1 | L2 | Linf |
| --- | ---: | ---: | ---: |
| y momentum | 0.13011 | 0.08247 | 0.06738 |
| z momentum | 0.25676 | 0.19960 | 0.14002 |
| Pressure | 0.25742 | 0.19787 | 0.13888 |

The z-momentum L1 ratio was 1.38403 with the scalar-scale screen; the increase
is absent in this run. The candidate's pressure maximum is 0.214495 Pa.
Maximum y/z errors now occur at fine full-volume cells (1,1,15)/(1,1,16),
still beside coarse/fine interfaces. Maximum pressure occurs at cut cell
(1,1,11). The error has not disappeared from the interface.

At the full-volume tracked probe, the second interval begins with a pressure
offset of -0.000704742 Pa, rather than -0.395790 Pa with the scalar-scale
candidate. The following predictor reaches -0.00141911 Pa, rather than
-1.170281 Pa. Together with unchanged admission counts, this supports the
identified fallback mechanism without claiming it explains every error.

## Remaining decision

The identical-binary regular-geometry control in
`build/eb-quartic-regular-20261023/` exceeded its 240-second bound. Its preserved
record remains failed. Default, momentum and velocity emitted complete 50-ns
tables; `partial-audit.json` checks these tables separately. Default is exact;
momentum/velocity y/z Linf ratios are 0.06738/0.14002, with pressure Linf
0.13887. Paired output is incomplete and is not counted as passed.

A recovery probe in `build/eb-quartic-paired-only-20261023/` changes only the
diagnostic mode loop to execute paired. It links the same frozen solver/helper
objects. The bounded recovery first requires its four cut-case tables to be
byte-identical to the completed candidate, then runs regular paired alone.
The recovery completes: cut takes 50.679 seconds and reproduces all four
paired tables byte-for-byte; regular paired takes 91.849 seconds. Its y/z
Linf ratios are 0.06738/0.14002 and pressure Linf is 0.13887. No reported
nonzero norm increases. The largest regular conservation residual across
the four settings is 6.06e-14. `recovered-comparison.json` in the regular
directory identifies all 16 accepted table hashes and their original or
recovery source; the original failed receipt is preserved. No full regular
sweep was relaunched.

After this control, a matched mesh study is needed to determine
whether the reduced interface error actually converges. General geometry,
serial/MPI lifecycle, restart, fixed-PeleC trajectory comparison and physical
qualification remain open; no production selector or checkpoint was changed.
