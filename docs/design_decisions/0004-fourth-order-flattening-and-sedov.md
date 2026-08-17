# Decision 0004: expose fourth-order slopes and flattening as independent controls

## Status

Accepted for the one-dimensional constant-`gamma` milestone.

## Context

PeleC applies a five-point limited PLM slope and an optional shock-flattening coefficient before characteristic tracing. Implementing both unconditionally would erase the already verified order-2/no-flattening baseline and make changes in smooth accuracy difficult to distinguish from changes in shock stabilization.

## Decision

Keep `reconstruction = "pelec_plm"` and add two orthogonal controls:

```fortran
plm_order = 2          ! or 4
use_flattening = .false.  ! or .true.
```

The defaults preserve all existing characteristic-PLM results. The order-4 option ports the regular-cell `plm_slope` formula. The flattening option ports the one-dimensional regular-cell pressure-jump and compression tests from `Godunov.H`.

Add a symmetric planar Sedov-type blast as the first strong-shock gate. It is deliberately checked through positivity, symmetry, conservation, shock location, and deterministic signatures rather than mislabelled as comparison with the spherical analytic Sedov-Taylor solution.

## Consequences

Benefits:

- order, flattening, and Riemann effects remain independently selectable;
- smooth convergence is measured without flattening;
- strong-shock robustness is measured with flattening;
- all older baselines remain in CI;
- future EB or multidimensional flattening paths can be compared against a stable 1D reference.

Limitations:

- outflow boundaries use constant primitive extension instead of AMReX patch data;
- no embedded-boundary neighborhood test is present;
- general-EOS thermodynamic characteristic terms are absent;
- the blast regression is planar and numerical, not a spherical analytic parity claim.

The next milestone should establish two-dimensional uniform-grid infrastructure and transverse Godunov verification before AMR is attempted.
