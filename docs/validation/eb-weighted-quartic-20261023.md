# Volume-weighted connected quartic candidate — 2026-10-23

The isolated guarded candidate now accepts physical child-volume weights and
checks connected donor geometry. **This is a local candidate test, not a
production halo policy, evolved convergence result or C2 qualification.**
It extends the equal-volume evidence in
`eb-guarded-quartic-20261023.md`; that earlier record remains historical.

## Scope and implementation

`tools/reference_probes/eb_guarded_quartic_mod.F90` recenters both the limited
baseline and quartic family to the parent's volume-weighted component mean.
One contraction coefficient applies to all eight children and all components.
Both existing EOS paths must admit every positive-volume child. Zero-volume
children receive the parent placeholder. Nonfinite, negative and all-zero
weights are rejected transactionally. Recentring establishes conservation;
it does not establish interpolation accuracy at arbitrary cut-cell centroids.

`eb_quartic_stencil_geometry_mod.F90` retains an axis only when all five axial
donors and four connecting faces are available and fluid-connected. It checks
every donor and face in the resulting rectangular subspace. Missing diagonal
coverage collapses the candidate to the parent; blocked axes sample the center
plane. Solid-cell NaNs are not read. Nonfinite active samples fail without
publishing outputs. Periodic mapping, distributed ownership and ancestor-time
sampling remain caller responsibilities.

## Independently completed checks

The direct probe passed **78 assertions** using GNU Fortran 13.3.0, bounds
checking and floating-point traps, linked against the existing fixed-core
Debug archive. Exact source, archive and executable SHA256 values, HEAD,
compiler and commands are recorded in
`build/eb-weighted-quartic-20261023/guarded-record.json`; output is in
`guarded.log`. The worktree was dirty at HEAD
`0a0e00ebe9138145fbb0e0e745665d0d768df009`; HEAD alone is not source identity.

- A plane wall leaves parent fluid volume 0.72 and child volume fractions
  0.44 or 1.0. Only the two wall-parallel axes remain active. Weighted child
  inventories reproduce the parent inventories while ignoring solid NaNs.
- The smooth 16-cell shear velocity error remains 1.736873665e-4 m/s with
  theta=1. This is a local analytic-average error, not evolved CFD error.
- Arbitrary unequal weights, including a zero-volume child, preserve weighted
  means. Invalid weights, active NaNs and missing diagonal donors exercise
  rejection or fallback as appropriate.
- Previous jump fallback, thermal/kinetic contraction, absent-species and
  transactional checks remain included in the 78 assertions.

Reproduction:

```sh
python3 tools/check_eb_guarded_quartic.py \
  --build-dir build/eb-eos-application-20261021/isolated-openmpi/registration \
  --work-dir build/eb-weighted-quartic-20261023
```

The runtime reports an array temporary for a component stencil. General
curved-wall accuracy, evolution, MPI, restart and fixed-PeleC trajectory
agreement are untested by these assertions. No production selector or default
was changed. A separate frozen-copy halo pilot evaluates integration cost and
initial rates before any evolved candidate run.

Mechanism byte provenance (`--bundle mechanisms/h2o2_full.json`), project
contract (`--source .`), Python compilation of both checkers and the pilot
drivers, and `git diff --check` pass. Initial invocations of the first two
checkers omitted required arguments and exited with usage status 2; corrected
commands passed without changing source or tolerances.
