# Parity and verification strategy

## Verification levels

PeleF uses four gates:

1. unit verification of algebraic kernels;
2. analytical or manufactured-solution verification;
3. parity against a pinned PeleC reference case when available;
4. conservation and deterministic field-signature verification.

Visual agreement is supplementary and never the sole acceptance criterion.

## One-dimensional gates

The existing suite retains independent checks for:

- EOS and primitive/conserved conversion;
- Rusanov and qualified PeleC-style Riemann fluxes;
- componentwise and characteristic reconstruction;
- order-2 and order-4 slopes;
- shock flattening;
- smooth entropy-wave convergence;
- Sod exact-solution error;
- Shu-Osher oscillation retention and field signatures;
- planar Sedov-type positivity, symmetry, conservation, and shock location.

Higher-order options never replace the lower-order baselines in CI.

## Two-dimensional algebraic gates

Directional flux verification checks that:

- rotating y-normal states to the x solver and rotating fluxes back is an exact involution;
- analytical y-direction Euler fluxes are recovered;
- equal-state Rusanov and PeleC-style y fluxes equal the physical flux.

The transverse-correction unit test checks:

- zero transverse divergence leaves a face state unchanged;
- a physical correction is applied without scaling;
- an excessive correction is reduced to a physical state with `0 < theta < 1`;
- invalid negative correction scales are rejected.

## Dimensional-reduction gate

A periodic entropy wave that varies only in x is repeated across y. One 2D CTU step is compared with the existing 1D characteristic Godunov step using the same timestep. Agreement to roundoff demonstrates that:

- y-direction flux divergence vanishes;
- transverse corrections do not contaminate dimensionally reduced flow;
- directional rotation preserves transverse momentum consistently.

## Isentropic-vortex convergence gate

The periodic isentropic vortex is run at `24 x 24`, `48 x 48`, and `96 x 96`. Both density refinement pairs must show at least order `1.8`. Current observed orders are approximately `2.278` and `2.276`.

The test also requires all periodic conserved integrals to remain within `2e-10`; current errors are approximately `6e-14`.

A differential run at `48 x 48` disables the transverse correction. The corrected solution must have density L1 error below 65% of the uncorrected error. Current values are:

```text
with transverse correction       5.3334803639e-4
without transverse correction    1.1493796865e-3
```

This prevents the transverse path from existing only nominally while having no measurable multidimensional effect.

## Application-level 2D gate

`pelef2d` runs the 64² vortex case from a namelist and writes CSV. An independent Python checker recomputes the translated analytical solution and verifies:

- density, pressure, and both velocity L1 errors;
- positive density and pressure;
- periodic mass, x-momentum, y-momentum, and energy conservation;
- the expected number of grid rows and finite output values.

## Reference-data policy

Pinned numerical signatures may be updated only with an explained numerical-method change. Analytical thresholds and conservation limits must not be relaxed merely to accept a regression. Future direct PeleC comparisons must record upstream commit SHA, input, build options, variable definitions, and comparison time.
