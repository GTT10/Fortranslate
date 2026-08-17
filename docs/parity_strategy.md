# Parity and verification strategy

## Verification levels

PeleF uses four distinct gates:

1. **Unit verification** checks algebraic kernels such as EOS conversion, reconstruction, flux consistency, and solver dispatch.
2. **Analytical verification** compares against exact or manufactured solutions where available.
3. **Reference parity** compares against a pinned PeleC case and configuration.
4. **Conservation verification** checks integral balances including known boundary fluxes.

The current milestone implements levels 1, 2, and 4 for Sod and the entropy wave. Shu-Osher uses unit, deterministic-signature, positivity, oscillation-retention, and boundary-balance gates. Direct PeleC field parity remains deferred until a reproducible reference artifact is pinned in `reference/`.

## Sod acceptance

For the canonical 400-cell case at `t = 0.2`, every solver path must produce finite positive states and conservative integral balances. The PLM/MC/PeleC-style path additionally requires:

- density L1 error against the exact Riemann solution at most `2e-3`;
- pressure L1 error at most `1.5e-3`;
- mass and energy errors at most `2e-12`;
- momentum agreement with the integrated boundary-pressure impulse to `2e-12`.

Before waves reach the domain edges, expected momentum change is

\[
\Delta P_x=(p_L-p_R)t.
\]

## Smooth convergence acceptance

The periodic entropy-wave test runs 40, 80, and 160 cells. Both the Rusanov and PeleC-style PLM paths must show observed density-L1 order of at least `1.8` across each refinement pair.

## Shu-Osher acceptance

The 800-cell `t = 1.8` regression requires:

- finite, positive density and pressure;
- mass, momentum, and energy balances including fixed boundary fluxes;
- sufficient maximum density and interaction-window density range;
- at least 15 local extrema in the interaction window;
- pinned density-squared and weighted-density integral signatures.

The two signatures detect broad field changes without storing a large generated reference file. They are regression evidence, not proof of exactness or direct PeleC parity.

## Reference data policy

Future PeleC reference artifacts must record:

- upstream repository and commit SHA;
- exact input file;
- build options and compiler;
- dimensionality and mechanism;
- output time and variable definitions;
- comparison script version.

Reference results must never be silently regenerated to make a failing test pass.
