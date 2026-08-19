# 0010 — Keep three reactive general-EOS Riemann baselines

## Decision

Add selectable `rusanov`, `hllc`, and `pelec` interface solvers to the
one-dimensional NASA7 reactive-flow path.

Rusanov remains the robustness baseline. HLLC uses frozen-composition outer
wave speeds and a contact wave, with every star state validated by the shared
NASA7 conserved-to-primitive conversion. The `pelec` option extends the
already verified constant-gamma acoustic-impedance/star-state interpolation to
an ideal-gas mixture by evaluating impedances and density corrections from the
local frozen-composition sound speed.

No solver silently falls back to another. Invalid wave ordering, non-positive
star states, composition failure, or thermodynamic inversion failure is
reported to the caller.

## Rationale

Rusanov is reliable but diffuses material interfaces and weak reaction-driven
waves. HLLC supplies a standard contact-resolving comparison independent of
PeleC. The qualified PeleC-style path provides a responsibility-level port of
the acoustic algorithm without pretending that the current seven-species
ideal-gas-mixture closure reproduces every PelePhysics general-EOS branch.

Keeping all three selectable allows errors from reconstruction, interface
physics, and chemistry splitting to be separated rather than hidden inside a
single preferred configuration.

## Verification consequences

- Equal states must return the physical flux for HLLC and PeleC-style paths.
- Stationary material contacts must have zero mass, energy, and species flux
  and momentum flux equal to pressure.
- Moving material contacts must return the upwind physical flux.
- A finite pressure jump must produce finite, correctly directed acoustic
  fluxes with exact species-mass closure.
- Smooth entropy waves must retain approximately second-order convergence.
- Contact advection and reactive-hotspot errors are compared directly against
  Rusanov at equal grid and reconstruction settings.

## Scope boundary

The HLLC and PeleC-style solvers use frozen composition across each acoustic
solve. They do not include a fully coupled thermodynamic eigenstructure,
real-fluid EOS, PeleC boundary scaling, rotating-frame terms, embedded
boundaries, or multidimensional transverse coupling.
