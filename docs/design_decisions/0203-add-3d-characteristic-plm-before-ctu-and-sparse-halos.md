# Decision 0203: add 3D characteristic PLM before CTU and sparse halos

## Context

Milestone `0.209.0` qualified periodic PCM/SSPRK2 hydro on a regular 3D grid
and a static two-level hierarchy, including deterministic replicated-state MPI
arithmetic. Moving directly to a PeleC-style multidimensional PPM/CTU update
would combine reconstruction, time tracing, transverse corrections,
reconstruction-width MPI halos, coarse/fine ghost fill, and flux-register
changes. A failure would not identify which contract was wrong.

## Decision

Add frozen-composition characteristic PLM as a selectable method-of-lines
spatial operator. Recover the full primitive field, rotate y- and z-normal
differences into the established x-normal characteristic basis, apply MC or
minmod limiting, scale the complete slope to preserve density, pressure, and
species, and convert both face states through the NASA7 EOS. Retain SSPRK2 and
publish the arithmetic mean of its two stage fluxes.

Use the same operator on the serial static hierarchy. Construct two fine ghost
layers at each SSPRK2 stage from coarse endpoint states interpolated in time.
Recover coarse primitives, form limited component slopes, evaluate them at the
fine-cell-center offset, apply one physical-state scale, normalize species,
and convert back to conserved state. Keep the existing substep-time and
face-area average, six-face reflux, and average-down sequence unchanged.

Fingerprint reconstruction and limiter in checkpoint schema 2. Keep the
replicated MPI AMR kernel PCM-only until a reconstruction-width halo and owner
work contract can be qualified independently.

## Consequences

The regular solver gains a second-order smooth-wave path, and the serial
static hierarchy can demonstrate materially lower error while retaining
roundoff conservation, exact synchronization, rollback, and restart identity.
PCM remains the default, so prior field arithmetic is unchanged. This is not a
claim of PeleC's unsplit PPM/CTU algorithm, distributed high-order AMR,
physical-boundary reconstruction, source/diffusion coupling across levels,
dynamic topology, or 3D EB.
