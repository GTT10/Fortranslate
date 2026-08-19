# 0009 — Isolate reactive general-EOS characteristic PLM

## Decision

Introduce composition-dependent reactive Euler flow through a new `pelef_reactive_1d` driver and `reactive_1d_mod`, leaving the established constant-`gamma` and passive-multispecies solvers unchanged.

The first high-order reactive reconstruction uses a frozen-composition ideal-gas-mixture characteristic basis. Hydrodynamic differences are limited in the `u-c`, contact/shear, and `u+c` waves using the local NASA7 mixture sound speed. Species mass fractions are limited componentwise and traced with the contact velocity. Rusanov remains the only reactive Riemann solver in this milestone.

Chemistry is coupled by reaction-hydro-reaction Strang splitting. The reaction substep keeps density, momentum, and total-energy density fixed, and obtains temperature from the unchanged specific internal energy after composition changes.

## Rationale

A direct replacement of the existing constant-`gamma` kernels would make it difficult to distinguish general-EOS errors from regressions in already verified hydro algorithms. A separate path provides differential tests and preserves lower-order baselines.

For an ideal-gas mixture at frozen composition, the local acoustic structure retains `u-c`, `u`, and `u+c`. This gives a useful second-order intermediate implementation without claiming the complete general-EOS characteristic formulation of PelePhysics/PeleC.

Keeping `rhoE` fixed during the adiabatic chemistry substep is required because the NASA7 species internal energies already contain formation-energy offsets. Adding a separate heat-release source would double-count chemical energy.

## Consequences

- The reactive path has independent conserved/primitive conversion, CFL, reconstruction, and flux tests.
- PCM and characteristic PLM remain selectable for direct error comparison.
- Homogeneous flow must reduce to the zero-dimensional reactor.
- The nonuniform hotspot must conserve global Euler invariants and produce a finite pressure/velocity response.
- General-EOS PeleC-style Riemann parity, pressure-dependent full chemistry, molecular transport, and multidimensional reacting flow remain separate milestones.
