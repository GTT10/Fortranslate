# Quartic candidate stage evolution — 2026-10-23

The frozen cached candidate completes the 16-grid, 50-ns cut-shear diagnostic
in **142.041 seconds**. It reduces both transverse momentum maximum errors,
but increases z-momentum L1 error. **Do not admit this candidate as a production
policy or claim C2 convergence.** This is one serial grid and interval, compared
with the saved production connected-direction baseline, not fixed-PeleC parity.

## Evidence and reproducibility

Artifacts: `build/eb-quartic-evolution-20261023/`. The driver `run.py` records
input/executable SHA256, the frozen candidate build record, command, process
identity, elapsed time and successful comparison status in `record.json`.
The binary is the same cached candidate used in the initial pilot. No
production numerical source, default, restart policy or tolerance was changed.

The default's four ordinary output tables and spatial traces are unchanged.
All four policies pass directional face-rate closure, 60 pre/post synchronization
probe pairs (separating raw register mismatch from applied redistribution), and
two-stage fine-cell momentum trajectory reconstruction. The executable also
enforces material, energy and momentum budgets. Exact receipts and residuals
are in `audit.json`, `sync-audit.json`, `trajectory-audit.json` and `audit.log`.
The largest normalized budget residual is 5.93e-14; the largest reconstructed
momentum-trajectory residual is 1.76e-16 in momentum-density units.

`compare.py` checks matched physical time and components, then compares norms
and maximum locations. `compare_stages.py` checks matched stage clocks before
comparing pressure and rates. Their outputs are `comparison.json` and
`stage-comparison.json`.

## Outcome at 50 ns

For paired pressure (the other two separated policies give the same ratios
to the shown precision), candidate/baseline errors are:

| Quantity | L1 ratio | L2 ratio | Linf ratio |
| --- | ---: | ---: | ---: |
| y momentum | 0.27725 | 0.30594 | 0.49745 |
| z momentum | **1.38403** | 0.88263 | 0.50597 |

The y-momentum maximum moves from cut cell (1,8,11) to full-volume cell
(1,8,15). The z-momentum maximum remains at full-volume cell (1,8,16).
All are fine cells beside coarse/fine interfaces. Reduced maxima therefore
do not establish uniformly reduced error, spatial convergence or EB accuracy.

## Stage attribution and next decision

At both tracked fine cells, the first predictor pressure offset falls from
-0.7908902 Pa to -0.000719674 Pa. At the cut-cell probe, the following y-rate
falls from 2530.31 to 2.30199 SI units. This supports the link from initial
halo mismatch to pressure generation and subsequent transverse forcing.

However, at full-volume probe (1,9,16), the next streamwise rate is again
-71867.58, despite an initial candidate rate of -65.39. At the start of the
second fine interval its pressure offset has reached -0.395790 Pa, compared
with -0.000711 Pa at the cut-cell probe. Each separated policy reports 960
adoptions initially and 464 on each subsequent fill. These observations
motivate tracing the candidate's admission/fallback decision at the specific
donors; aggregate counts alone do not prove which rejection rule caused the
local change. Do not weaken screening or EOS admission on this evidence.

The next useful check is donor-local candidate admission through the existing
stage clocks, followed by a regular-geometry control. A larger mesh campaign
is premature while this switching and the increased z L1 error are unresolved.
