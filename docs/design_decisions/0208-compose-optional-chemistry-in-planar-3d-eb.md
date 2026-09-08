# Decision 0208: compose optional chemistry in the planar 3D EB application

## Context

The `0.214.0` application qualified a serial single-level exact-axis-plane EB
hydro workflow with an active-cell CFL selector and either FluxRedist or
StateRedist. The regular 3D reactive path already provided transactional
cell-local chemistry, but the EB application had no source split. Treating
species components as conserved in that application would also make its
invariant gate incorrect once reactions were enabled.

The selected 3D geometry has covered storage that must not be altered by a
source stage, and cut-cell hydro stabilization must remain the same as the
existing inert path. A coupled implementation therefore needs a private
candidate, an active-cell mask, and a gate that distinguishes fluid invariants
from chemical species conversion and the physical wall-normal impulse.

## Decision

Extend the 3D cell-local chemistry wrapper with an optional `[nx,ny,nz]`
active mask. It slices each z plane into the existing 2D chemistry API and
rejects a shape mismatch transactionally. Every chemistry plane operates on a
private state and temperature candidate; the caller is published only after
all planes succeed.

Add an EB split driver with the composition

```text
R(dt/2)-H(dt)-R(dt/2)
```

where `R` is the masked elementary chemistry update and `H` is the existing
PCM EB Euler update selected by `flux_redist` or `state_redist`. Construct the
active mask from `cell_type /= covered`, check covered state and temperature
identity after each stage, and publish the candidate only after final
physicality checks. Invalid solver, tolerance, reaction, geometry, or stage
requests return the original state and temperature unchanged. When chemistry
is disabled, call the same hydro routine with the same inputs so the old
StateRedist/FluxRedist arithmetic remains byte-identical.

Add validated namelist controls for `chemistry_enabled`, relative tolerance,
and absolute tolerance. Load the frozen seven-species elementary H2/O2/N2
mechanism used by the existing regular reactive path. In the public invariant
gate, retain mass, total energy, and the two tangential momenta; leave the
wall-normal momentum as a reported wall exchange. With chemistry enabled,
replace individual-species conservation checks by H/O/N elemental totals and
species-closure/physicality checks.

Provide a focused x/y/z-mask unit gate and separate public inert/reacting
StateRedist cases. The inert CSV must equal the existing `0.214.0`
StateRedist CSV byte-for-byte. The checker validates the common 480-row
geometry/EOS contract, exact covered rows, active species change, fluid
invariants, and H/O/N totals against the baseline.

## Consequences

The focused Debug unit reports maximum split mass/energy/tangential drift
`2.5597e-16`, elemental drift `1.1045e-15`, and active species change
`1.5629e-06`. The public reacting StateRedist run takes two steps to
`5.0e-5`, changes an active species by `5.3319005211e-02` in conserved
species density, and reports mass/energy/tangential drift `4.718e-16` and
H/O/N drift `8.366e-15`. The baseline and inert CSVs are byte-identical.

This qualifies only the serial, single-level, inviscid, first-order PCM
exact-axis-plane workflow using the frozen elementary mechanism and public
StateRedist case. Its timestep is selected from the state before the first
reaction half step, so no general reaction-aware CFL guarantee is made. It
does not qualify molecular transport, configurable or stiff production
mechanisms, higher-order EB reconstruction, general geometry,
AMR/reflux/regridding, checkpoint/restart, MPI ownership, or physical
validation against PeleC or experiment.
