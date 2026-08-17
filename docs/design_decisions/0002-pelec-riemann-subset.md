# Decision 0002: port a selectable ideal-gas PeleC Riemann subset

## Status

Accepted as an intermediate, explicitly qualified port.

## Context

PeleC's `Source/Riemann.H` is coupled to PelePhysics thermodynamics, species arrays, rotating-frame terms, and the surrounding Godunov method. Waiting for all of those dependencies before testing any of the Riemann algorithm would leave reconstruction and flux changes entangled.

The existing Rusanov solver is robust and verified but overly diffusive for contact and shock-density-wave problems.

## Decision

Implement a separate `riemann_pelec_mod` containing the constant-`gamma`, single-species reduction of PeleC's acoustic star-pressure/star-velocity estimate, star-density correction, shock treatment, and inward/outward wave interpolation.

Add an explicit dispatch layer with two selections:

- `rusanov`;
- `pelec`.

Do not silently fall back to Rusanov when the selected PeleC-style solver fails. Return failure so the calling timestep can reject the state.

## Verification

The new path is gated by:

- equal-state consistency with the physical Euler flux;
- stationary-contact behavior;
- pinned canonical Sod interface state and flux values;
- invalid-selection rejection;
- second-order periodic entropy-wave convergence;
- exact-solution Sod errors and conservation;
- positive, conservative Shu-Osher evolution with pinned field signatures.

## Consequences

Positive consequences:

- Riemann behavior can be tested independently of characteristic reconstruction;
- Rusanov remains a differential-diagnosis baseline;
- the Shu-Osher problem exposes excessive dissipation and oscillation loss;
- future general-EOS and multispecies work has a stable dispatch interface.

Limitations:

- no species mass-fraction state;
- no PelePhysics general EOS;
- no rotating-frame energy terms;
- no boundary-condition velocity multiplier;
- no characteristic tracing or multidimensional transverse correction.

The next step is characteristic-variable projection and tracing behind a third reconstruction path, not an expansion of this decision's parity claim.
