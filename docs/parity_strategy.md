# Parity and verification strategy

## Verification levels

1. **Unit verification**: algebraic identities, projection inverses, flux consistency, and invalid-selection rejection.
2. **Analytical verification**: exact Riemann solutions and smooth convergence.
3. **Deterministic regression**: pinned field signatures for problems without a compact exact solution.
4. **Reference parity**: comparison to a pinned PeleC build and output dataset.
5. **Conservation verification**: integral balances including known boundary fluxes.

The current milestone implements levels 1, 2, 3, and 5. Direct PeleC dataset parity remains pending a pinned upstream executable configuration and reference artifact.

## Characteristic-PLM gates

- primitive-slope to characteristic-wave projection and exact inverse;
- isolated acoustic-wave synthesis;
- zero-Courant face-state identity;
- contact/shear preservation of normal velocity and pressure;
- rejection of invalid negative `dt/dx`;
- periodic entropy-wave order greater than `1.7` on two refinement pairs;
- 400-cell Sod density L1 at most `1.6e-3` and pressure L1 at most `1.0e-3`;
- positive density and pressure;
- roundoff-scale Sod integral balances;
- Shu-Osher positivity, oscillation retention, extrema count, boundary-flux balances, and pinned density signatures.

## Reference data policy

Future PeleC reference artifacts must record upstream commit SHA, exact inputs, compiler and build options, EOS/mechanism, output variable definitions, output time, and comparison-tool version. Reference data must not be regenerated merely to make a failing test pass.

The current deterministic signatures are implementation-regression gates, not evidence that the entire Fortran result is numerically identical to PeleC.
