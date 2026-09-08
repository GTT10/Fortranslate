# Historical PR #11: pressure-dependent chemistry and dense implicit verification

This preserves the unique decision note from commit
`24cf4ec35017b27ed9903fb1b45b7e4f3ca8543d`, originally named
`docs/design_decisions/0009-pressure-dependent-chemistry-and-implicit-reactor.md`.
The current decision 0009 covers a different topic; no file is overwritten.
This is historical rationale, not the current capability or qualification page.

PR #11 was closed as superseded after verifying that all 17 chemistry,
mechanism, generator, application/case and numerical verification files were
byte-identical in `788be768a53b7fb6119b67215345e18f47975e0b`, which is an
ancestor of main `bc584fc65e7e945b98475cd1c0394b483c5cc685`. The old PR commit
itself is not an ancestor and is not represented as merged. Its branch is kept.

## Original context

The four-reaction H2/O2 subset verifies reversible elementary kinetics but
excludes the third-body and falloff reactions that make the complete small
mechanism stiff. Extending the flow solver before a complete zero-dimensional
parity gate would make chemistry, thermodynamics and transport errors difficult
to separate.

## Original decision

1. Preserve the seven-species/four-reaction executable and its high-precision Cantera gate.
2. Generate a separate ten-species/29-reaction mechanism from normalized JSON.
3. Represent elementary, third-body and falloff reactions through one runtime reaction record.
4. Evaluate Troe broadening and collider efficiencies inside the shared kinetic kernel.
5. Assemble the fixed-temperature concentration Jacobian analytically from reaction products and rate derivatives.
6. Form a reduced constant-energy Jacobian by eliminating the closure species and including the temperature response implied by `u(Y,T)=u0`.
7. Use dense backward Euler/Newton with line search as an in-tree verification integrator.
8. Estimate time-discretization error with one full step versus two half steps; accept a Richardson-extrapolated state when physical, otherwise retain the two-half-step state.

## Original consequences

The complete Cantera H2/O2 mechanism can be verified without introducing a
third-party ODE dependency. This dense solver is intentionally limited to
small mechanisms and is not a performance architecture for hydrocarbon
chemistry. CVODE/ARKODE, sparse Jacobians, direct mechanism parsing and
chemistry-flow coupling were separate milestones at the time. The full
mechanism must not silently replace the elementary regression; both remain
in CI. See [current status](../current_status.md) for subsequent capabilities.
