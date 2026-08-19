# Decision 0009: pressure-dependent chemistry and dense implicit verification

## Context

The four-reaction H2/O2 subset verifies reversible elementary kinetics but excludes the third-body and falloff reactions that make the complete small mechanism stiff. Extending the flow solver before a complete zero-dimensional parity gate would make chemistry, thermodynamics, and transport errors difficult to separate.

## Decision

1. Preserve the seven-species/four-reaction executable and its high-precision Cantera gate.
2. Generate a separate ten-species/29-reaction mechanism from normalized JSON.
3. Represent elementary, third-body, and falloff reactions through one runtime reaction record.
4. Evaluate Troe broadening and collider efficiencies inside the shared kinetic kernel.
5. Assemble the fixed-temperature concentration Jacobian analytically from reaction products and rate derivatives.
6. Form a reduced constant-energy Jacobian by eliminating the closure species and including the temperature response implied by `u(Y,T)=u0`.
7. Use dense backward Euler/Newton with line search as an in-tree verification integrator.
8. Estimate time-discretization error with one full step versus two half steps; accept a Richardson-extrapolated state when it remains physical, otherwise retain the two-half-step state.

## Consequences

- The complete Cantera H2/O2 mechanism can be verified without introducing a third-party ODE dependency.
- The dense solver is intentionally limited to small mechanisms and is not a performance architecture for hydrocarbon chemistry.
- CVODE/ARKODE, sparse Jacobians, direct mechanism parsing, and chemistry-flow coupling remain separate milestones.
- The full mechanism cannot silently replace the elementary regression; both remain in CI.
