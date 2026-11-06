# Recursive EB MOL stage capture — 2026-10-23

## Scope

The serial common and recursive SDC entry points now accept a read-only
`observe_mol` callback. On the native held-energy MOL spatial path it receives
stage, depth, interval start time, interval width, input tree, per-node local
diagnostics and the tree operator's returned integrated fluxes. Unselected
nodes have unallocated diagnostic arrays. No observer means no diagnostic
array allocation. The ordinary numerical policy and default options remain
unchanged.

These are **attempt** observations, not accepted-step receipts: later failure
or a retry does not retract callbacks. Stage zero evaluates the old state;
positive stages evaluate the new iterate. The interval-start time is relative
to the current recursive invocation; nested replay invocations restart that
origin. Native replay/cache-only paths need not emit fresh spatial samples.
This interface is currently serial and held-energy MOL only; it does not claim
complete MPI, CTU, transport, replay or retry trace coverage.

The tree diagnostic snapshots are taken before cross-patch limiting; returned
fluxes are after it. They are not coarse/fine synchronized time-integrated
registers. Explicit pre/post recursive reflux snapshots are still outstanding.

The reference probe's `--audit-shear-stage-flux` mode enables the existing
connected-directions shear policies and records six directional faces around
two fine-node cells: the x/cut-plane intersection at the shear extremum and
the x/z coarse/fine intersection. At root size 16 these are (1,8,11) and
(1,9,16), matching the previously localized paired-policy maxima. They are
fixed probe cells, not a new search for moving maxima. The CSV records EOS
cell pressure, both reconstructed face pressures, raw/centroid/tree momentum
fluxes, local spatial and wall rates, and the input momentum density.
Closed faces require aperture masking. Flux units are instantaneous in this
CSV; tree integrated fluxes are divided by the recorded interval width.

## Bounded verification

Artifacts: `build/eb-stage-capture-20261023/`.
Pre-edit source copies, build logs, bounded driver, independent CSV audit and
source/executable/input SHA256 identities are retained there.
The driver reuses the frozen 16-cell 50 ns cut-shear input, checks all sixteen
ordinary result tables against the completed production-path baseline, and
independently reconstructs directional momentum divergence plus wall rates
from each set of six recorded faces.

Build and bounded calculation completed successfully. The cut case took
55.05 seconds; the independent audit returned zero. All **16 ordinary CSV
tables match exactly** at parsed field level, and all four policies produced
8 probe-stage samples each (192 directional-face rows in total). The largest
relative momentum closure residual is 1.618e-16. Local-to-tree instantaneous
face-flux differences are at most 1.456e-11, consistent with multiplying and
dividing the integrated flux by dt; this does not measure recursive reflux.
The first build failed on an added Fortran line exceeding 132 columns; the
continuation was corrected and the failed log retained.

For both probe cells in all four policies, the first fine predictor lowers
cell pressure by 0.7908902 Pa. In its next spatial evaluation, the cut-plane
probe has y-momentum rate +2530.31251 from the y-face difference alone; the
full-volume intersection probe has y rate -2530.31251 and z rate -1265.42306.
The latter comes entirely from the z-face difference: the low-face normal
momentum flux is 101324.209110 Pa and the high-face flux 101324.604555 Pa.
The high face reconstructs 101324.209110 Pa on the fine side and 101325 Pa on
the donor side. Thus the observed initial transverse acceleration is already
present in the raw Riemann face flux, before centroid interpolation or the
tree boundary limiter. This identifies its location in the operator, not the
upstream interpolation defect that created the pressure mismatch.

At the cut probe, the large z-face pressure contribution is balanced by the
wall contribution; it is not a spurious net z acceleration there. The sampled
first predictor momentum increments agree with dt times the sampled local
rate (largest absolute residual 2.063e-16 in x; transverse residuals at most
1.356e-20). These are first-predictor checks, not a claim that corrected SDC
steps use forward Euler. Subsequent corrected stages are retained in
`audit.json`. Fine closure and root reflux remain to be measured separately.

Mechanism provenance, project contract, Python syntax and diff whitespace
checks passed. This is serial 16-grid, 50 ns diagnostic evidence, not MPI,
8/32 convergence or fixed-PeleC qualification. General source modifications
remain in dirty main at HEAD `0a0e00ebe9138145fbb0e0e745665d0d768df009`;
source/executable/input identities for this run are in `record.json`.
