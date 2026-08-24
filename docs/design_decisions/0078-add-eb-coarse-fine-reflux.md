# Decision 0078: add EB coarse/fine reflux

## Context

EB-volume-weighted average-down synchronizes state inside a refined rectangle,
but it cannot correct the different numerical fluxes used by independently
advanced coarse and fine levels. A regular flux register stores that mismatch
on cells outside the fine patch. At an embedded boundary, both face aperture and
cell fluid volume enter the correction, and applying the full volume-normalized
increment to a small cut cell recreates a stiffness problem.

AMReX's EB flux register accumulates aperture-weighted coarse and fine fluxes.
During reflux it keeps a cut cell's volume-fraction share of the correction and
redistributes the remaining share to connected neighbors. Contributions that
land geometrically inside the fine region are transferred back to the fine
level.

## Decision

Add a serial `amr_eb_flux_register_2d` for one aligned rectangular patch. Store
the accumulated state correction on the full coarse layout and require coarse
and restricted-fine open-face measures to agree along all four patch sides.
Coarse and fine calls carry their own `dt`, so any valid subcycle schedule can
accumulate before reflux.

Use the exterior coarse cell's fluid volume to normalize both its coarse face
term and the physical sum over fine subfaces. Opposite signs on low and high
sides reproduce the finite-volume orientation. Skip a side coincident with the
physical domain boundary because there is no coarse/fine interface there.

Apply a regular-cell correction directly. For a cut cell, multiply the raw
state correction by `kappa`, add that stable increment to the source cell, and
scatter its `1-kappa` share equally in state after normalizing by connected
neighbor fluid volume. Map a neighbor below the fine rectangle to all of its
fine children. Reset the register only after both candidate level arrays are
finite. The reactive wrapper copies the register, restores covered states,
recovers both temperature fields, and commits only after every active EOS call
succeeds.

## Consequences

Matching coarse/fine fluxes cancel across cut and regular interfaces even with
fine subcycling. A nonmatching correction remains globally conservative after
cut-cell stabilization and fine-recipient transfer. A failed reactive EOS check
retains both input levels and the still-populated register for diagnosis or
retry.

This does not yet implement EB prolongation, coarse/fine ghost interpolation,
level advancement, dynamic regridding, multiple patches, deeper levels, or MPI
ownership.

Primary references:

- [AMReX `EBFluxRegister`](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EBFluxRegister.cpp)
- [AMReX two-dimensional EB flux-register kernels](https://github.com/AMReX-Codes/amrex/blob/development/Src/EB/AMReX_EBFluxRegister_2D_C.H)
