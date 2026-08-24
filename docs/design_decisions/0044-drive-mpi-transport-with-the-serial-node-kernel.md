# Decision 0044: drive MPI transport with the serial node kernel

## Context

Recursive molecular transport has the same hierarchy coupling as hydro, but
its stable fine schedule uses `r^2` substeps and its face flux combines two
SSPRK2 stages. The `0.51.0` MPI hydro scheduler already establishes owner-only
patch execution without duplicating finite-volume formulas.

## Decision

Extract the serial one-patch SSPRK2 molecular-transport update as a shared AMR
kernel. It retains viscosity, Fourier conduction, mixture-averaged species
diffusion, barodiffusion, correction velocity, species-enthalpy flux, physical
boundary handling, thermodynamic recovery, and effective face-flux output.

Walk the patch tree collectively, but call that kernel only on each patch
owner. Broadcast the owner's interval-start state, complete effective face
flux, accepted patch, and transport counter. Replicas then execute the
existing `r^2` child schedule, parent and sibling ghost construction, shared
fine/fine diffusive flux correction, flux-register accumulation, reflux,
average-down, and temperature recovery. Restore the synchronized pre-call
solution on every rank after any collective rejection.

## Consequences

Serial and MPI transport use one numerical kernel, while patch work and
parabolic subcycle counts are exclusive to owners. Deep branched and adjacent
cross-owner cases retain serial field parity and composite conservation.

The bridge still broadcasts complete patch and flux arrays and repeats
hierarchy synchronization on replicas. This decision does not provide a
combined distributed full-physics transaction, sparse rank-local storage,
stage-synchronous point-to-point exchange, or regrid migration.
