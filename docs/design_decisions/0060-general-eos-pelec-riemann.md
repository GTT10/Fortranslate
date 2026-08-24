# Decision 0060: extend the PeleC acoustic solver through the NASA7 mixture EOS

## Context

The constant-`gamma` `riemann_pelec_mod` qualified the acoustic star-state and
wave-interpolation algorithm in isolation. Reactive flow later added NASA7
thermodynamics and a separate HLLC solver, but `pelec` was still unavailable to
the conserved `(rho, rho*u, rho*v, rho*w, rho*E, rho*Y_k)` path.

PeleC `Source/Riemann.H` evaluates sound speed and internal energy through its
EOS, carries species densities through the origin, star, and interpolated
interface states, and closes the final flux from that interface state.

## Decision

Add `reactive_pelec_flux_x` beside the independent reactive Rusanov and HLLC
solvers. Recover both primitive states and frozen sound speeds through the
NASA7 mixture EOS, form the acoustic impedance estimate of star pressure and
velocity, and choose the origin species-density vector from the upwind side.
Average the two species-density vectors when the interface is stationary.

Apply the pressure correction to every origin species density, rebuild and
validate the star state through the EOS, interpolate between origin and star
states using the inward/outward wave speeds, and rebuild total energy from the
final interface density, composition, pressure, and velocity. Force the last
species flux to close exactly to total mass flux.

The public dispatcher and reactive 1D/2D configuration readers accept
`pelec`. The existing y-direction momentum rotation reuses the x-normal
kernel. Because this Fortran API recovers primitive variables from total
energy rather than receiving already reconstructed primitives, its stationary
test includes the EOS inversion roundoff scale.

## Consequences

The reactive path now exercises the same acoustic star-state, species-density,
shock/rarefaction, and wave-interpolation structure as PeleC while retaining
Rusanov and HLLC as independent diagnostic baselines. The claim is limited to
the implemented NASA7 ideal-gas mixture. Rotating-frame energy, auxiliary and
linear advected variables, embedded-boundary coupling, and other PelePhysics
EOS models are not included.

Verification covers equal states, stationary and moving composition contacts,
a nonuniform shock pair, species closure, invalid input, x/y rotation, and a
full 40/80/160-cell composition-wave convergence sequence.
