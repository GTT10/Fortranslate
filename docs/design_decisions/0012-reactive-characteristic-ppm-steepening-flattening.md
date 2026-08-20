# Decision 0012: add time-traced reactive characteristic PPM as a separate path

## Context

The reactive solver already had two qualified higher-order alternatives:

- frozen-composition characteristic PLM with a time-centered normal predictor;
- componentwise primitive PPM used as a semidiscrete operator with SSPRK3.

The latter verified monotone parabolic interpolation but did not exercise the
normal characteristic profile integration in PeleC `Source/PPM.H` and
`Source/PPM.cpp`.  A separate path was needed so the PeleC-style predictor
could be tested without redefining the existing PPM regression baseline.

## Decision

Add

```text
reconstruction = "characteristic_ppm"
```

while retaining `pcm`, `characteristic_plm`, and `ppm` unchanged.

For each primitive component, the new path:

1. constructs five-point PPM edge values with the van-Leer and
   Colella--Sekora limiters used by PeleC;
2. optionally blends the parabola toward the cell center using the
   one-dimensional PeleC shock-flattening coefficient;
3. integrates the parabola over the regions swept by `u-c`, `u`, and `u+c`;
4. carries species mass fractions and transverse velocities on the middle
   wave;
5. projects density, normal velocity, and pressure over the existing
   frozen-composition acoustic/contact basis;
6. normalizes species and reconstructs face energy through the NASA7 EOS;
7. performs one conservative time-centered Godunov update.

A profile is rejected when a local characteristic Courant number exceeds one.
No selected HLLC solve is silently replaced by Rusanov.

## Contact steepening

Add an opt-in Colella--Woodward-style detector based on:

- density change relative to pressure change;
- opposite-sign neighboring density curvatures;
- a minimum one-percent relative density jump.

The detector blends density and mass-fraction edges toward neighboring
MC-reconstructed values and clips them to the adjacent-cell range.  The
applied strength is capped at `0.5`.  Full-strength simultaneous density and
composition steepening over-compressed the moving material-contact regression
with the current frozen-composition HLLC star states; the cap is therefore an
explicit qualification rather than an undocumented stability adjustment.

## Shock flattening

Add the one-dimensional regular-cell formula from PeleC `Godunov.H`, including
its shifted neighboring detector.  It is opt-in and applies only to
`characteristic_ppm`.  Smooth regions and expansions retain an unflattened
parabola, while strong compression can reduce the coefficient to zero.

## Consequences

The new path is second order on both the smooth entropy and composition waves.
It sharpens the 200-cell H2/N2 material contact relative to componentwise PPM,
and bounded steepening further reduces the H2 L1 error from
`7.35878653e-5` to `2.50998077e-5`.

The periodic pressure-ratio-three shock remains positive, conservative, and
inside the initial pressure extrema with flattening enabled.  The smooth
Gaussian reacting hotspot converges and remains far better than PCM, but its
current error is larger than both characteristic PLM and componentwise PPM.
That result is retained as a regression: this decision establishes the
qualified normal characteristic predictor, not universal superiority.

Still excluded are the complete PeleC/PelePhysics general-EOS eigensystem,
independent internal-energy characteristic tracing, multidimensional
transverse corrections, and a pressure-dependent stiff chemistry mechanism in
CFD.
