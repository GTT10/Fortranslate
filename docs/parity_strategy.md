# Parity and verification strategy

## Verification levels

PeleF uses four distinct gates:

1. **Unit verification** checks algebraic kernels such as EOS conversion and flux consistency.
2. **Analytical verification** compares against exact or manufactured solutions where available.
3. **Reference parity** compares against a pinned PeleC case and configuration.
4. **Conservation verification** checks integral balances including known boundary fluxes.

The current milestone implements levels 1, 2, and 4 for the Sod problem. Direct PeleC numerical parity is deferred until a reproducible PeleC reference dataset is pinned in `reference/`.

## Current Sod acceptance criteria

For the canonical 400-cell case at `t = 0.2`:

- all values are finite;
- density and pressure remain positive;
- density L1 error against the exact Riemann solution is at most `3e-2`;
- pressure L1 error is at most `3e-2`;
- total mass error is at most `2e-12`;
- total energy error is at most `2e-12`;
- total momentum agrees with the integrated boundary-pressure impulse to `2e-12`.

The expected momentum is not zero. Before waves reach the domain edges, the fixed left and right exterior states exert a net pressure impulse

\[
\Delta P_x=(p_L-p_R)t.
\]

## Reference data policy

Future PeleC reference artifacts must record:

- upstream repository and commit SHA;
- exact input file;
- build options and compiler;
- dimensionality and mechanism;
- output time and variable definitions;
- comparison script version.

Reference results must never be silently regenerated to make a failing test pass.
