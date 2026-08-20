# Decision 0013: establish a separate periodic reactive 2D CTU path

## Context

The constant-gamma solver already had a verified 2D CTU-style update, while the composition-dependent reactive solver was limited to one dimension. Directly modifying the old 2D path would risk changing its pinned constant-gamma regressions and would still leave species and energy transverse corrections split across incompatible closures.

## Decision

Add `reactive_2d_mod` and `pelef_reactive_2d` as a separate uniform-grid path. Reuse the qualified 1D NASA7 conserved/primitive conversions, frozen-composition characteristic PLM relations, and Rusanov/HLLC kernels. Rotate x/y momentum and velocity components for the y-normal solve.

Compute provisional x/y fluxes and apply the transverse half-step correction to the complete conserved face vector. Use an EOS-admissibility bisection to scale that vector only when necessary. Do not correct species separately from mass or total energy.

Keep the first milestone periodic and restrict reconstruction to PCM or characteristic PLM. Multidimensional PPM tracing, contact steepening, shock flattening, and physical boundaries remain separate work. Chemistry is the existing seven-species, four-reaction cell reactor and is Strang split around the CTU hydro step.

## Consequences

A y-uniform field reduces to the 1D characteristic-PLM/HLLC update at roundoff. The exact oblique entropy wave converges and records a measurable transverse-correction improvement. The periodic vortex and reacting hotspot remain positive, conservative, and composition-closed.

This establishes a qualified reactive regular-grid CTU subset, not full PeleC multidimensional Godunov/PPM parity.
