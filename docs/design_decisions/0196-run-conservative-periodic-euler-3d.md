# Decision 0196: Run conservative periodic Euler in three dimensions

## Context

Milestone `0.202.0` established independent x/y/z coordinates and directional
Riemann semantics, but no field update exercised the three flux directions
together. Introducing reacting physics, AMR, or EB before a regular-grid
conservation and time-integration gate would combine too many failure modes.

## Decision

Represent a uniform periodic 3D field without ghost cells. Store the upper
face flux for every cell in each direction and obtain the lower flux by
periodic predecessor indexing. Use PCM face states and the existing Rusanov or
PeleC-style Riemann selector. Choose the timestep from the sum of directional
spectral radii and compose two spatial evaluations with SSPRK2.

Keep both stage states private and publish only a completely physical second
stage. Qualify the implementation by exact x/y/z reduction to the established
1D PCM/SSPRK2 solver, a smooth diagonal entropy-wave convergence study,
roundoff conservation, and a separately checked public application CSV.

## Consequences

PeleF now has a runnable conservative 3D Euler application with a precise
serial single-level boundary. PCM is intentionally first order; no claim is
made for a multidimensional high-order Godunov predictor. General-EOS reacting
evolution, chemistry, transport, AMR, restart, MPI decomposition, EB, LES,
particles, and production I/O remain later workstreams.
