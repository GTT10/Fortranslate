# Parity and verification strategy

## Verification levels

PeleF uses four distinct gates:

1. **Unit verification** checks algebraic kernels such as EOS conversion, limiter behavior, linear reconstruction, and flux consistency.
2. **Analytical verification** compares against exact or manufactured solutions where available.
3. **Reference parity** compares against a pinned PeleC case and configuration.
4. **Conservation verification** checks integral balances including known boundary fluxes.

The current milestone implements levels 1, 2, and 4. Direct PeleC numerical parity is deferred until a reproducible PeleC reference dataset is pinned in `reference/`.

## Sod acceptance criteria

For the canonical 400-cell case at `t = 0.2`, both reconstruction paths require:

- all values finite;
- positive density and pressure;
- total mass error at most `2e-12`;
- total energy error at most `2e-12`;
- total momentum agreement with the integrated boundary-pressure impulse to `2e-12`.

The piecewise-constant baseline permits density and pressure L1 errors up to `3e-2`. The PLM/MC path has tightened limits of `4e-3` and `3e-3`, respectively. These tighter limits prevent a silent fallback of the entire second-order path to first order.

The expected momentum is not zero. Before waves reach the domain edges, the fixed left and right exterior states exert a net pressure impulse

\[
\Delta P_x=(p_L-p_R)t.
\]

## Smooth convergence gate

A periodic entropy wave with constant pressure and velocity is advanced to `t = 0.1` on 40, 80, and 160 cells. The density L1 convergence order must be at least `1.8` on both refinement pairs. The verified orders are approximately `2.05` and `2.09`.

Periodic verification wraps both ghost-cell states and PLM slopes. Wrapping only the states introduces a first-order defect at the periodic interface and is therefore explicitly covered by the convergence test.

## Reference data policy

Future PeleC reference artifacts must record:

- upstream repository and commit SHA;
- exact input file;
- build options and compiler;
- dimensionality and mechanism;
- output time and variable definitions;
- comparison script version.

Reference results must never be silently regenerated to make a failing test pass.
