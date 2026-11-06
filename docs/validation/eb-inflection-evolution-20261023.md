# General inflection-screen candidate evolution

Scope: isolated serial paired-policy 16-cell/50-ns shear. The optional general
inflection hierarchy is enabled for all available candidate donor families;
there is no coordinate override. Production is unchanged and C2 remains open.
Exact source/binary identity, commands, logs and tables are in
`build/eb-inflection-evolution-20261023/`.

The run completed in 49.871 seconds, return code zero, with normalized budget
residual 5.83314e-14. Stage face-rate closure, 27 sync pairs and all six fine
momentum trajectories pass. Initial analytic rows match the baseline. The
baseline here is the saved vector-roundoff quartic candidate, not fixed PeleC.

Relative to that baseline, y momentum L1/L2/Linf ratios are
0.841111/1.014755/0.823037; z ratios are 0.601108/0.715284/0.854323; pressure
ratios are 0.672744/0.735447/0.853966. Of all 17 components and three norms,
only y L2 increases. No tolerance was relaxed to hide that increase.

The y maximum moves to fine (1,3,12), z to (1,3,16), pressure to (1,3,11).
The original target predictor pressure still drops from 0.111649534 to
0.00354680694 Pa. Thus local causal improvement survives without the coordinate
override, but does not establish uniformly improved error. Measured candidate
halo adoption counts are 3264,3264,3264,2752 over the four fills: the last stage
still changes admission.

`initial_profile.py` predicts the initial halo error transition using analytic
sine and sine-squared parent averages: fine y=1/2 errors fall to +/-0.000873185
m/s, while y=3/4 retain linear errors +/-0.0233016 m/s. This is consistent with
the shifted maximum near the admission boundary; it is not a runtime
attribution of every donor or of the final-stage switch. Full comparison and
maximum locations are retained in `comparison.json`.

Decision: do not promote the candidate. Further work must explain the residual
admission transition and y L2 degradation, and establish matched-grid/longer-time
behavior before general integration. No new large calculation is justified
solely by the local maximum improvement. Python syntax, project contract and
`git diff --check` pass. No job remains active.

The subsequent [recorded admission audit](eb-inflection-admission-20261023.md)
separates initial total-energy extremum rejection from final-fill
near-roundoff momentum switching in the sampled stencils.
