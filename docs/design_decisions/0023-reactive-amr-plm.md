# Decision 0023: flux-consistent PLM for reactive AMR

## Status

Accepted for PeleF `0.28.0`.

## Context

PCM established a robust conservation baseline for the runnable reactive AMR
path but diffuses moving contacts. Higher spatial order must preserve mixture
state validity and must supply a time-integrated flux consistent with reflux.

## Decision

Reconstruct primitive density, three velocities, pressure, and species mass
fractions with the existing MC or minmod limiter. Fall back to the cell center
for invalid density or pressure and clip/renormalize species before converting
face states to conserved variables.

Advance each level with SSPRK2. Return the mean of the two stage fluxes as the
effective conservative flux accumulated by the AMR flux register. At periodic
physical boundaries use one shared numerical flux. For each fine substep use
coarse ghost data interpolated to that substep's midpoint.

## Consequences

- smooth contacts are less diffusive than the PCM AMR baseline;
- face thermodynamics and species closure remain explicit gates;
- reflux matches the actual SSPRK2 conservative update;
- midpoint fine ghosts improve coarse-time boundary accuracy;
- characteristic tracing, PPM/WENO, and higher-order coarse-time variation
  within both SSPRK stages remain follow-on work.
