# 0236: Add static 3D AMR molecular transport

## Status

Accepted for `0.244.0`.

## Context

The `0.243.0` applications qualify fixed and configure-time selected chemistry
for a static, strictly-interior, periodic, two-level Cartesian 3D AMR
hierarchy in serial and rank-local sparse MPI execution.  Regular-grid
molecular transport already has the mixture-transport flux contract, but the
AMR operators did not advance that transport across a coarse/fine interface.
The missing path must preserve the existing transactional state ownership,
temperature recovery, sparse ownership, and chemistry/hydrodynamics ordering.

The qualification case is one rectangular fine patch with an integer refinement
ratio.  It does not change the static-topology or periodic-boundary assumptions
of the existing AMR applications.

## Decision

Add mixture-averaged molecular transport to the serial and sparse-MPI static
two-level AMR operators.

- Use the established mixture transport records and fluxes: viscosity,
  Fourier heat conduction, mixture-averaged species diffusion with the
  correction velocity, species enthalpy transport, and the optional
  barodiffusion contribution (which requires species diffusion).  The
  transport timestep is the minimum of the coarse stable timestep and
  `r**2` times the fine stable timestep, where `r` is the integer refinement
  ratio.
- Advance transport with transactional SSPRK2.  Each Euler stage advances the
  coarse level and then takes `r**2` fine subcycles of `dt/r**2`; the final
  transport candidate is the SSPRK2 average.  EOS temperature recovery,
  physical-state checks, and candidate publication remain part of the same
  acceptance transaction.
- At every fine subcycle, construct coarse/fine ghost states by linear time
  interpolation between the provisional coarse start and end states followed
  by limited componentwise primitive PLM spatial reconstruction. Molecular
  transport uses this PLM ghost fill at the interface unconditionally; the
  hydro operator's existing PCM/characteristic-PLM selection is separate.
- Accumulate each of the lower and upper `x`, `y`, and `z` interface fluxes
  over all fine subcycles and child faces.  Divide by `r**2` for the time
  average, apply six-face time-averaged diffusive reflux to the coarse level,
  and conservatively average down the covered coarse cells from the fine
  level before the final temperature recovery.
- Apply a common coarse/fine full-face species-positivity limiter on each of
  the six AMR faces.  The common face theta is constrained to `[0,1]` from
  the signed replacement-flux loss and the coarse species mass left after the
  other coarse faces, the averaged fine interface flux is scaled by that
  theta, and the removed interface flux receives a compensating,
  equal-and-opposite correction on the adjacent fine boundary cells.  This
  limiter and its compensating fine-side flux correction are applied
  independently to every lower/upper `x`, `y`, and `z` face before reflux and
  average-down.
- In sparse MPI, face ownership remains rank-local for state and flux work.
  The six full-face theta arrays are combined with collective minima, and
  boundary corrections are applied only by the ranks that own those fine
  cells.  A rank with zero fine planes still participates in every consensus,
  limiter reduction, validity/acceptance reduction, and transport collective;
  it performs no empty fine-state dereference.
- Compose the full reactive operator as
  `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`, with `T` the transactional AMR
  transport operator, `H` the hydro operator, and `R` the chemistry half-step.
  Any failed stage, EOS recovery, limiter bound, conservation check, or
  collective agreement discards private candidates and publishes none of the
  state, temperature, clock, or diagnostics.
- Route fixed and generated selected-mechanism frontends through the same
  mechanism-independent serial and sparse-MPI application drivers.  Selected
  bundle, ordered species/thermodynamic/reaction/transport records, and the
  generated integrator are validated and agreed collectively before
  data-dependent sparse work.
- Reject any transport-enabled checkpoint interval or restart request before
  opening an output file.  The fixed schema-2 and selected schema-3
  chemistry-context contracts do not bind the gas-transport database and
  complete parabolic operator as a qualified restart contract, so transport
  plus checkpoint/restart is an explicit startup rejection rather than an
  incomplete restart mode.

## Consequences

Serial and sparse-MPI static AMR now have one mixture-transport contract,
including six-face interface positivity, reflux, average-down, temperature
recovery, and transactional rollback.  The sparse implementation can retain
its rank-local ownership and still give the same accepted result when some
ranks own no fine planes.  The selected dispatch path covers transport as
well as the fixed path, while transport-enabled checkpoint/restart remains
clearly unsupported.

This decision qualifies only a static, strictly-interior, periodic, two-level
Cartesian hierarchy with one rectangular fine patch and the existing
mixture-averaged transport model.  Dynamic topology or regridding, physical
coarse/fine boundary conditions, 3D embedded-boundary AMR, deeper or general
multi-patch hierarchies, scalable/distributed I/O, external PeleC or other
cross-code field parity, performance/scalability qualification, and physical
validation remain outside this decision.  Transport model extensions such as
Stefan--Maxwell, Soret, or Dufour transport are also separate work.
