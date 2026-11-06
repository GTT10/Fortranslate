# Pressure-consistent spatial representations, 2026-10-21

A diagnostic pressure reconstruction now evaluates two alternatives with the
production held-MOL spatial operator. This follows the conservative-stencil
[pressure audit](eb-initial-stencil-20261021.md). Both alternatives retain
interpolated density, momentum and species, without the nearly-PCM slope
contraction rejected in that audit.

The fine-only alternative interpolates donor EOS pressure with the existing
scalar monotonic kernel and rebuilds ghost internal/total energy with the
mixture EOS. It leaves all fine valid cells byte-identical. This is a spatial
ghost representation, **not conservative energy prolongation**.

The paired alternative additionally volume-averages fine EOS pressures into
covered coarse pressure donors. Both coarse and fine spatial operators use
this same pressure representation. Coarse covered/ghost energy is rebuilt
only in operator-local arrays; persistent coarse restriction and fine state
remain unchanged. Uncovered coarse valid gas is asserted byte-identical.
No application flag, checkpoint identity, production integrator or MPI path
has been changed.

## Initial planar result

| Maximum interface quantity | Original | Fine-only | Paired |
| --- | ---: | ---: | ---: |
| Fine mass flux (kg/m2/s) | 1.18231e-4 | 3.83938e-14 | 3.83938e-14 |
| Area-weighted mass flux mismatch (kg/m2/s) | 1.95225e-5 | 6.12462e-5 | 2.86588e-14 |
| Normal momentum flux mismatch (Pa) | 0.127029 | 0.0700050 | 1.45519e-11 |
| Total-energy flux mismatch (W/m2) | 64.9799 | 203.569 | 9.51696e-8 |

Fine-only pressure repair removes the fine disturbance but worsens the
coarse/fine transfer mismatch. The paired representation removes the
stationary interface mismatch to roundoff while retaining density, momentum
and species reconstruction. Thus coarse and fine pressure reconstruction
must be treated together. This supports further investigation of a paired
spatial-pressure policy, not immediate production adoption.

The test is at time zero, before any chemistry or time integration. Actual
cell inventories are unchanged at this point because no step is taken.
Global evolved conservation, pressure-work/held-energy consistency, re-reflux,
reaction/transport coupling, retry/rollback, MPI parity, restart identity and
asymptotic convergence still require evidence. Altering auxiliary coarse
energy inside a spatial operator is a policy change relative to the pinned
reference. Do not claim that initial flux agreement establishes those gates.

## Reproduction and evidence

The extended reference probe accepts
`--audit-initial --separate-momentum --pressure-pair`.
It writes the original, `_pressure_candidate`, and
`_pressure_pair_candidate` pressure/face CSVs. The paired pressure CSV
contains the actual paired operator representation, including covered coarse
cells. Fluxes are intensive SI, divided by dt, positive-coordinate; the
comparator uses open-area weighting and sums four fine faces per coarse face.

The final verification uses an immutable copied probe and four controls:
planar imported/uniform and regular-restricted imported/uniform. The regular
control redoes covered coarse restriction for its equal-volume geometry.
`tools/compare_eb_pressure_ghosts.py --run-record ... --output ...` verifies
binary/input identity and stationary pressure/mass/energy interface bounds.
It refuses incomplete records; it does not count earlier exploratory runs
as final evidence.

Artifacts: `build/eb-pressure-ghost-20261021/`, with `final-run.json`,
`final-comparison.json`, compiler/build logs, frozen binary and source hashes.
The earlier 12-control exploratory run used the shared build target and is
excluded from the final qualification count. Final acceptance belongs only
to the immutable copied binary.

All four immutable-binary controls passed. Across them, the maximum paired
pressure jump is7.28e-11 Pa and maximum energy-flux mismatch is1.32e-7 W/m2.
Six Python checks passed (four interface/restriction checks plus two evidence
admission checks), with zero failures/skips. Provenance, syntax, project
contract and whitespace checks passed. GNU Fortran13.3.0 Debug was used;
no toolchain/configuration change or full suite is claimed.
