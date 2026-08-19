# Decision 0011: add monotone primitive PPM as an independent reactive path

## Context

The verified reactive one-dimensional path already had a frozen-composition
characteristic PLM predictor and general-EOS HLLC flux. That path has a small
error coefficient on smooth waves, but a separate parabolic reconstruction is
useful for sharper material interfaces and as groundwork for later PeleC-style
PPM characteristic tracing.

## Decision

Add `reconstruction = "ppm"` without replacing the existing PCM or
characteristic-PLM paths.

The PPM path:

- constructs fourth-order interface candidates from four neighboring
  cell-centered primitive states;
- clips each interface candidate to the adjacent-cell range;
- applies the standard cellwise parabolic monotonicity constraints;
- re-normalizes bounded species mass fractions at both cell edges;
- converts the limited general-EOS primitive edge states back to conserved
  variables;
- evolves the semidiscrete flux divergence with SSPRK3;
- can be paired independently with Rusanov or HLLC.

The wider stencil is filled explicitly by periodic wrapping or outflow constant
extension. PPM never silently changes the selected Riemann solver.

## Consequences

The implementation resolves the moving H2/N2 material contact more sharply
than characteristic PLM with the same HLLC solver. It remains monotone and
conservative on smooth entropy/composition waves and reacting hotspots.

The current PPM path is componentwise in primitive variables and does not yet
include PeleC's full characteristic integration, contact steepening, shock
flattening, or general-EOS internal-energy wave. Consequently, characteristic
PLM remains more accurate for the present smooth Gaussian hotspot and entropy
wave. The regression suite records that distinction instead of assuming PPM is
universally superior.
