# Decision 0119: rebuild sparse EB topology transactionally

## Context

The public sparse MPI EB clock can advance a fixed multipatch hierarchy, while
the established serial EB regrid already provides conservative average-down,
PCM initialization, and old/new overlap retention. Sparse callers have no safe
way to replace the child topology and its ownership metadata together.

## Decision

Add an explicit-plan collective regrid that treats the distribution, sparse
payload container, and replicated geometry template as one transaction.
Validate and agree on the refinement ratio and collection controls, materialize
the current owner fields, and apply the qualified serial multipatch EB regrid.
Build a fresh deterministic subcycle-weighted distribution from the resulting
ordered children, scatter only each new owner's numerical payload, validate the
one-copy sparse result, and then commit all three objects together.

Any invalid control, materialization failure, serial regrid failure, ownership
failure, or sparse scatter failure returns without changing the caller's old
distribution, fields, or template. Report `changed` only from committed patch
count, bounds, or refinement-ratio differences.

## Consequences

Sparse MPI EB applications can now move, resize, create, or remove child
patches through an explicit serial-compatible plan and resume direct owner
timesteps on the rebuilt hierarchy. Normal physics intervals remain free of
replicated child fields.

This first correctness bridge materializes the hierarchy on every rank during
the infrequent regrid event. Owner-local tagging, direct old/new overlap
transfer, scheduled regrid integration with the public clock, checkpointing,
and parallel output remain future lifecycle work.
