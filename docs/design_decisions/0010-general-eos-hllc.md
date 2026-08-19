# 0010 — Use HLLC as the contact-resolving general-EOS intermediate

## Decision

Add a selectable HLLC flux to the reactive one-dimensional solver while
retaining Rusanov as an explicit robustness baseline. Do not label the HLLC
path as PeleC Riemann parity.

## Rationale

The NASA7 reactive state already supplies a thermodynamic pressure and frozen
sound speed for each reconstructed face. HLLC can use those side-local values
without requiring a constant-gamma assumption across the interface, and it
restores an explicit contact wave that Rusanov smears.

The implementation uses conservative Davis outer-wave bounds, the standard
HLLC contact-speed relation, and Rankine--Hugoniot star states. Species mass
fractions and tangential velocities remain frozen across each outer wave. The
last species flux is closed against the total mass flux.

Every star state must pass the existing NASA7 conserved-to-primitive recovery.
If the selected HLLC path cannot construct a physical state, the step fails;
there is no silent fallback to Rusanov.

## Verification boundary

The accepted evidence is:

- equal-state physical-flux identity;
- stationary and moving contacts with different H2/N2 composition;
- species-flux closure;
- second-order smooth composition-wave convergence;
- lower discontinuous material-contact error than Rusanov;
- a positive, conservative reactive-hotspot application run.

The smooth composition-wave test does not show HLLC outperforming Rusanov at
every resolution. The claim is therefore contact resolution, not universal
error reduction.

This remains a frozen-composition ideal-gas-mixture HLLC solver. It is not the
PeleC/PelePhysics general-EOS Riemann algorithm, and it does not include
real-gas thermodynamics, characteristic thermochemical coupling, or molecular
transport.
