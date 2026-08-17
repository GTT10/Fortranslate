# Parity and verification strategy

## Verification levels

PeleF uses four gates:

1. unit verification of algebraic kernels;
2. analytical or manufactured-solution verification;
3. parity against a pinned PeleC reference case when available;
4. conservation and deterministic field-signature verification.

Visual agreement is supplementary and never the sole acceptance criterion.

## Fourth-order slope and flattening gates

The slope unit test checks:

- exact reconstruction of linear data with the five-point order-4 stencil;
- the expected PeleC order-2 expression;
- zero slope at a discrete extremum;
- direct multiplication by the flattening coefficient;
- rejection of unsupported orders.

The flattening test checks:

- coefficient one for smooth data;
- coefficient zero for a strong compressive pressure step;
- coefficient one for the same pressure step in expansion.

## Smooth-convergence gate

The periodic entropy wave is run at 40, 80, and 160 cells with order-4 characteristic PLM and no flattening. Both refinement pairs must show at least order 1.9. The observed orders are approximately 2.17 and 2.18; time integration limits the complete method to second order.

## Strong-shock gate

The planar Sedov-type regression verifies:

- all CSV values finite;
- positive density and pressure;
- exact discrete reflection symmetry within tolerance;
- mass, momentum, and energy errors below `2e-10`;
- sufficiently strong compression and pressure peaks;
- shock radius within a pinned interval;
- density-squared and pressure integrals within a `5e-6` relative tolerance.

Pinned signatures for the current 800-cell case are:

```text
density-squared integral = 1.4513811919127926
pressure integral        = 2.0606157390942492
shock radius             = 0.131875
```

The checker accepts explicit replacement references on its command line, but reference changes must be reviewed and explained rather than silently regenerated.

## Non-regression policy

New higher-order or PeleC-style options do not replace old paths. Existing PCM, componentwise PLM, order-2 characteristic PLM, and Rusanov cases continue to run in CI so regressions can be localized to reconstruction, flux evaluation, or time integration.
