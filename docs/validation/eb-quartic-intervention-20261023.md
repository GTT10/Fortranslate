# Bounded donor-screen causal intervention

This is an isolated 16-cell/50-ns serial experiment, not a production policy
or fixed-PeleC qualification. Artifacts are in
`build/eb-quartic-intervention-20261023/`; exact compile/link commands, source
hashes and binary identity are in `build-record.json`.

Only root donor families at x=N,y=1 bypass the curvature screen. Geometry
availability, volume-weighted conservation, held-energy EOS and total-energy
EOS admission remain active. The override is in artifact copies, not the
production source. The existing candidate is the comparison baseline.

The run completed in 50.171 seconds, return code zero. Maximum normalized
budget residual is 5.68781e-14. Finite values, stage face-rate closure, all 27
sync pairs and six fine-probe momentum trajectories pass. Initial analytic
rows match the baseline. Four recorded target fills use the override with
theta=1. Initial halo velocity is 0.979470818621283 m/s, matching the saved
quartic replay within 1e-12.

At the targeted first predictor, pressure offset decreases from 0.111649534
to 0.00354680694 Pa (ratio 0.0317673). At the former y maximum, fine (1,1,15),
final transverse momentum decreases from 8.52888765e-6 to 1.23026655e-7
(ratio 0.0144247). At fine (1,1,16), z momentum decreases from 8.72277621e-6
to 2.77786960e-7 (ratio 0.0318462). This intervention supports a causal role
for the rejected high-order donor family in the local pressure/transverse
error; it does not establish that it is the only contribution.

Global y/z/pressure Linf remains unchanged. Their maxima move to the unmodified
opposite side: (1,16,15), (1,16,16), and (1,16,11). y-momentum L1 improves to
0.961786 of baseline, but L2 increases to 1.005137. z L1/L2 ratios are
0.899748/0.936415; pressure L1/L2 ratios are 0.898674/0.936631. All 17 component
norm comparisons and maximum locations are retained in `comparison.json`.

The hardcoded override is not an adoptable fix. The next candidate needs a
general criterion that recognizes smooth inflections while rejecting
underresolved discontinuities, with the existing conservation and EOS guards
retained. That criterion needs independent smooth/discontinuous-state evidence
before another trajectory test. C2 remains open; MPI/restart and physical
qualification are not inferred. Python syntax, project contract and diff
whitespace checks pass for the tooling and record updates.
