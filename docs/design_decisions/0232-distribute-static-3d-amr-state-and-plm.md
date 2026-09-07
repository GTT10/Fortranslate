# 0232: Distribute static 3D AMR state and characteristic PLM

## Status

Accepted for `0.240.0`.

## Context

The `0.209.0` MPI 3D AMR path assigned arithmetic to x slabs but replicated
both numerical levels on every rank. The `0.210.0` serial hierarchy added
characteristic PLM, which requires two reconstruction planes and
time-interpolated coarse context around fine-patch boundaries. Removing the
replicated fields therefore requires an ownership and halo contract before
dynamic topology, AMR source terms, or 3D EB can be added safely.

## Decision

Introduce `mpi_amr_sparse_reactive_3d_mod` as the production MPI stepping
boundary for the qualified static, periodic, two-level hierarchy.

- Coarse cells use deterministic uneven contiguous x slabs.
- Fine x planes belong to the owner of their parent coarse plane. A rank may
  own no fine cells and still participates in every collective contract.
- Each hierarchy stores only local conserved state and temperature. Two
  periodic coarse ghost planes are propagated across ranks, including the
  one-cell-slab case. Active fine ranks exchange internal x ghosts; exterior
  fine ghosts are reconstructed from local coarse start/end halos.
- PCM and characteristic PLM use rank-local SSPRK2 candidates. PLM reuses the
  serial characteristic reconstruction and limited coarse prolongation on
  only the required slab planes.
- Fine interface fluxes are accumulated in deterministic owner order, then
  applied through local six-face reflux, average-down, and temperature
  recovery. Caller state is published only after collective acceptance.
- Public operations fingerprint patch, ownership arrays, NASA7 data, floating
  controls, and solver/reconstruction/limiter text before data-dependent
  communication. Rank-dependent valid metadata is a collective rejection.
- Existing schema-2 checkpoint and deterministic CSV contracts remain
  unchanged. Root gathers only at initialization, checkpoint, restart, or
  output; rank count and ownership do not enter persistent state.

The older replicated module remains an internal compatibility oracle. The
public MPI application uses the sparse hierarchy.

## Consequences

The public PCM and characteristic-PLM cases are byte-identical to serial at
one, two, four, and eight ranks. A two-rank checkpoint resumes independently
on four and eight ranks with the same exact final coarse and fine files. At
eight ranks four ranks own no fine cells, directly exercising empty fine
payloads and active-neighbor discovery.

Root gather/scatter is rank-neutral compatibility I/O, not scalable I/O.
This decision does not add dynamic regridding, AMR chemistry or transport,
physical coarse boundaries, CTU/PPM, 3D EB AMR, performance scaling, or
external PeleC field parity.
