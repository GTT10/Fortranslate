# Decision 0014: reuse the qualified characteristic PPM normal predictor in 2D CTU

## Context

PeleF 0.14.0 established a periodic reactive two-dimensional CTU path using
PCM or frozen-composition characteristic PLM. PeleF 0.13.0 had already
qualified a one-dimensional characteristic PPM normal predictor with optional
PeleC-style shock flattening and bounded Colella--Woodward contact steepening.
The next step was to verify that predictor in both coordinate directions
without claiming the complete multidimensional PPM algorithm used by PeleC.

## Decision

Add `reconstruction = "characteristic_ppm"` to the reactive 2D driver.
For each cell and coordinate direction:

1. gather the periodic seven-cell primitive stencil;
2. rotate the y-normal stencil into the existing x-normal ordering when needed;
3. apply the same five-point PPM edge reconstruction used by the qualified 1D
   path;
4. optionally apply the one-dimensional PeleC shock-flattening coefficient;
5. optionally apply the bounded contact steepener to density and species;
6. integrate each parabola over the `u-c`, `u`, and `u+c` domains of
   dependence;
7. project density, normal velocity, and pressure over the frozen-composition
   acoustic/contact basis while carrying species and transverse velocities on
   the middle wave;
8. rotate y-normal states back and pass all face states through the NASA7 EOS;
9. apply the existing full-conserved-state CTU transverse half-step correction;
10. recompute the final directional HLLC or Rusanov flux.

The PPM-only controls remain invalid for PCM and characteristic PLM. A local
characteristic Courant number above one is rejected. The bounded steepening
strength remains capped at `0.5`, matching the qualified 1D path.

## Verification contract

The 2D characteristic-PPM path is accepted only when:

- x-normal and y-normal reductions match the corresponding 1D update below
  `3e-12` relative difference;
- exact oblique density and H2/N2 composition waves converge under refinement;
- the transverse correction has a measurable signature and does not degrade
  the 32-square density result;
- synthetic cell stencils activate shock flattening and contact steepening
  without losing positivity or composition closure;
- bounded steepening sharpens a periodic 2D material contact;
- an oblique pressure-ratio-three shock remains positive and conservative,
  and the flattening switch has a deterministic solution signature;
- a reacting 2D hotspot remains positive, conservative, and composition-closed.

## Consequences

The new path verifies a PeleC-style **normal** characteristic PPM predictor in
both coordinate directions and its coupling to PeleF's conservative full-state
CTU correction. It does not reproduce PeleC's complete multidimensional PPM
transverse/corner characteristic tracing. The distinction is recorded in the
README, numerical-method documentation, mapping table, and parity gates.

The next numerical-physics milestone is molecular transport. Full
multidimensional PPM corner coupling can be added later against a dedicated
PeleC parity case without blocking viscosity, heat conduction, or species
diffusion work.
