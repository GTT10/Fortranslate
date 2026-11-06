# Falsification of a monotone inflection-screen extension

The bounded donor intervention established that rejecting a smooth inflection
can contribute to the remaining coarse/fine shear error. It did not establish
a replacement smoothness classifier. This record evaluates one proposed
extension before implementation into the Fortran candidate.

Proposed rule: keep the existing curvature-sign test, and additionally accept
a five-point line if all four first differences have the same resolved sign
and its ratio-two quartic child averages remain between the adjacent coarse
averages on their respective sides. Roundoff allowance is 16 machine epsilons
times the line magnitude; no physical tuning parameter is introduced.

`build/eb-inflection-screen-20261023/evaluate.py` uses analytic finite-volume
averages of sine, cubic and tanh functions, plus four offsets of a step and a
single-cell peak. The rule accepts the targeted 16-cell sine inflection and
rejects the four sharp step placements. However it also accepts tanh profiles
with widths 0.125, 0.25 and 0.5 coarse cells. These are subcell transition
witnesses: monotonicity and bounded child averages do not establish resolved
smoothness. Exact inputs, children and differences are saved in
`screen-evaluation.json`.

Decision: reject this extension as a smoothness classifier. It is not added
to the production solver or candidate Fortran. No trajectory calculation was
launched. This does not prove those tanh interpolations would cause instability;
it falsifies the intended ability to distinguish resolved inflections from
underresolved sharp transitions. A further proposal must address this evidence
and retain volume conservation, EOS and discontinuity checks. C2 is still open.

## Difference-hierarchy candidate

A second experimental extension retains strict monotonicity and adjacent-child
bounds, and additionally requires max(abs(third differences)) <=
min(abs(first differences)) plus the existing roundoff allowance. This is a
conservative finite-difference hierarchy heuristic, not a proven shock sensor.
No coordinate, mesh size or target-cell selector appears in the rule.

`evaluate_hierarchy.py` evaluates all parent positions at 8/16/32/64 cells using
both sine momentum and sine-squared kinetic-energy averages. The combined
screen admits 0/8/32/64 parents, respectively; at 16 it adds parents 1,8,9,16
to the original four. All nine shifted tanh transitions at each width
0.125/0.25/0.5/1 cells are rejected; widths 2 and 4 are admitted at all nine
positions. These finite samples do not establish universal shock robustness.

The isolated helper `tools/reference_probes/eb_guarded_quartic_mod.F90` now
accepts optional `inflection_hierarchy`, default false. Existing callers retain
the original curvature screen. The change is not connected to production halo
filling. Existing geometry, weighted conservation and both EOS admission checks
remain unchanged. The standalone Fortran probe passes 164 assertions, including
all previous 94, two smooth-shear inflection families with each conserved mean
checked, unit rescaling and 36 shifted underresolved tanh states. Exact commands,
compiler and source/library hashes are recorded in
`build/eb-inflection-screen-20261023/fortran/guarded-record.json`.

Command: `python3 tools/check_eb_guarded_quartic.py --build-dir
build/eb-eos-application-20261021/isolated-openmpi/registration --work-dir
build/eb-inflection-screen-20261023/fortran`.
Mechanism provenance, Python syntax, project contract and diff checks pass.
The next gate is a bounded same-condition candidate trajectory comparison
using this general rule, with all norms, conservation and maximum locations
checked. No trajectory or MPI/restart qualification is claimed here.
