# Screened/admissible quartic family candidate — 2026-10-23

## Decision and boundary

An isolated gas-family candidate now screens for abrupt curvature changes,
falls back to the existing separate-momentum limited interpolation, and
contracts all eight children toward one admissible baseline with a single
coefficient when necessary. The smooth shear improvement survives the focused
checks. **This is not a production halo policy or a C2 completion result.**

The family uses equal child volumes. General cut-cell volume weighting,
125-cell connected donor availability, boundary/ancestor-time sampling,
serial/MPI integration and evolved shear errors are still untested. In
particular, an equal-volume mean check must not be presented as conservation
for children with different fluid volumes. The curvature screen is an
explicit candidate heuristic, not an established general shock detector.

## Implementation

`tools/reference_probes/eb_quartic_kernel_mod.F90` contains the isolated raw
kernel shared with the previous reference probe. Refactoring it out of that
probe retained the **119-case exact C++ agreement**, parent-average checks and
polynomial reproduction, independently rerun in this artifact directory.

`eb_guarded_quartic_mod.F90` adds two layers:

1. On every five-point coordinate line of every packed gas/held-energy
   component, the three second differences must be all consistent in sign,
   or all within a floating-point cancellation allowance. A mixture of zero
   and resolved curvature or a sign change declines the quartic candidate.
   The allowance is 16 machine epsilons times the largest line magnitude;
   it is not a tuned physical scale. Declined stencils use the current
   separate-momentum limited eight-child family. If that family is not
   admissible, it falls back to the parent state in every child.
2. The candidate is blended with its admissible baseline by a **single theta
   for every component and every child**. Trial theta is halved until the
   complete family passes, with an explicit theta=0 fallback. This is an
   admissible backtracking choice, not a maximal-theta optimization. Families
   with different component means are rejected before blending. Nonnegative
   species and positive density are required. Each child must pass both the
   existing held-energy EOS and the existing total-energy EOS, so excessive
   kinetic energy cannot be hidden behind a valid held-energy temperature.
   No positive-total-energy assumption is imposed on formation energies.

The packed layout is conserved gas followed by held internal-energy density;
this is not yet the production halo's gas/temperature/held packing. Results
and theta publish only on successful admission. Nonfinite data or mismatched
family means leave caller outputs unchanged. No production numerical source,
application selector or checkpoint format was changed.

## Completed evidence

Artifacts: `build/eb-guarded-quartic-20261023/`.
The fixed-core GNU Fortran 13.3.0 Debug library supplies the existing EOS and
limited baseline routine. Its SHA256, probe/kernel sources, compiler, exact
commands and executable identity are recorded in `guarded-record.json`.
The direct diagnostic compilation enabled bounds and floating-point checks.

The guarded probe passes **37 focused assertions**:

- At the 16-cell analytic shear extremum, theta=1 and the velocity-average
  error remains **1.736873665e-4 m/s**. All component parent means are retained.
- A velocity jump declines the smooth screen and returns the limited baseline.
- An excessive-kinetic-energy trial contracts the entire family to
  theta=**0.00048828125** and passes both existing EOS paths.
- An excessive thermal-energy contrast contracts to theta=**0.0625**.
- Negative perturbations to an absent species force exact theta=0 baseline
  restoration, even when density/species sum and sibling mean are preserved.
- NaN input and unequal baseline/candidate means are refused transactionally.

The runtime reports an array temporary when extracting one component's 3D
stencil; this is a diagnostic performance limitation, not an accuracy failure.
No CFD run was launched. These assertions do not establish shock accuracy,
general positivity under evolution, EB conservation or MPI correctness.

```
python3 tools/check_eb_quartic_stencil.py \
  --header build/eb-quartic-stencil-20261023/AMReX_Interp_3D_C.H \
  --work-dir build/eb-guarded-quartic-20261023
python3 tools/check_eb_guarded_quartic.py \
  --build-dir build/eb-eos-application-20261021/isolated-openmpi/registration \
  --work-dir build/eb-guarded-quartic-20261023
```

Mechanism provenance, project contract, Python syntax and whitespace checks
pass. Next integration work must explicitly provide physical child weights and
available connected stencils before a bounded shear evolution can count as
an EB candidate test. Default pinned-reference behavior remains unchanged.
