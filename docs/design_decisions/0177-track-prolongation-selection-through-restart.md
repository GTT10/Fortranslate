# Decision 0177: Track prolongation selection through restart

## Context

Fixed-depth fresh runs can select limited-linear prolongation, but the method
is absent from checkpoint identity. Arbitrary-depth serial and sparse-MPI
regrids therefore remain PCM-only. Lifting either restriction without tracking
the selection could restart a solution under a different mesh-initialization
contract.

## Decision

Advance all fixed-depth checkpoint schemas to version 3 and the shared
serial/sparse patch-tree fingerprint to schema 4. Store and compare
`prolongation_method` before accepting topology or fields. Pass the selected
dispatcher through arbitrary-depth iterative planning and final rebuilding.

In sparse MPI, include the method in collective control consensus. The parent
owner performs the selected prolongation once, after which the existing direct
child-owner transfer and transactional overlap-retention path continue.

## Consequences

Linear prolongation is restart-safe across every public reactive 2D EB AMR
lifecycle, including cross-rank sparse restart. Older fixed-depth schema-2 and
fingerprinted patch-tree schema-3 files are intentionally rejected. PCM stays
the default, and cut/topology-mismatched or inadmissible linear parents still
fall back to PCM.
