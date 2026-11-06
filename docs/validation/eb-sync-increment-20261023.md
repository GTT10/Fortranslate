# Closed stage/closure momentum accounting — 2026-10-23

## Completed scope

The serial 16-grid, 50 ns cut-shear diagnostic now accounts separately for
spatial-stage face fluxes, raw reflux mismatch, and the applied closure
increment after re-reflux redistribution. The solver numerical policy is
unchanged. All four policies retain exact parsed equality of the 16 ordinary
output tables and all four spatial-stage traces with their saved baselines.
This is bounded operator attribution, not C2, MPI or physical qualification.

`close_held_prepared_eb_tree_3d` now optionally returns the raw gas increment
and the gas increment actually applied, publishing both only on successful
closure. The recursive closure observer forwards them to the probe. These
are conserved-density increments, not rates or fluxes. The pre-closure CSV
uses zero placeholders; only phase 1 contains the measured increments.

Artifacts: `build/eb-sync-increment-20261023/`. The GNU Fortran 13.3.0 Debug
selected-runtime build and probe both completed. Probe elapsed time was
53.93 seconds; the independent audit took 0.087 seconds. The source, binary,
input and command identities are in `record.json`. HEAD was dirty main
`0a0e00ebe9138145fbb0e0e745665d0d768df009`; the recorded file hashes, not HEAD
alone, identify the calculation. Build command:

```
cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration \
  --target test_eb_tree_reacting_reference -j 4
python3 build/eb-sync-increment-20261023/run.py
```

## Independent accounting

All 60 pre/post probe pairs pass two independent links at the existing
1e-12 gate: changed face-register divergence versus raw increment, and
applied increment versus observed leaf-state change. Maximum absolute
residuals are respectively 1.655e-24 and 1.788e-16. The full two-stage fine
momentum trajectories, including all closure increments, agree to at most
2.507e-16. No gate was relaxed. The earlier face-only audit failure remains
preserved in `build/eb-sync-capture-20261023/`; its missing redistribution
term is now measured directly, not silently absorbed into tolerance.

For paired policy, the fine cut probe's x-momentum change is -0.00354916123:
-0.00354903946 comes from integrated spatial rates and -1.21771621e-7 from
closure redistribution. At both fine probes, all closure y/z momentum
increments are zero in this first root interval. The fine cut y-momentum
change is +1.26587816e-4; the full-volume intersection y/z changes are
-1.25482436e-4 and -6.22949970e-5. These transverse changes are fully reproduced
by the sampled spatial rates. The coarse cell above the fine corner does
receive -7.74306966e-8 in z momentum during closure. Thus closure is not absent
from the solution, and later-time feedback through coarse halos remains open.

## Initial velocity mismatch identified without another CFD run

`tools/analyze_eb_shear_stage_origin.py` reads the saved initial halos and
fluxes in `build/eb-shear-initial-20261022/` and the current stage traces.
All 48 sampled initial momentum faces agree within 1e-12 relative tolerance.
The script preserves SHA256 identities of all inputs in `halo-origin.json`.
This establishes a local link to those older halo samples, not whole-run
identity between the historical and current configurations.

At both fine probes, the periodic coarse donor supplies u=9.744953584 m/s,
while the adjacent fine average is 9.935868511 m/s. Pressure has no initial
jump. The coarse velocity stencil at the lower-y probe is
[8.261372739, 9.744953584, 9.744953584] m/s. Its centered slope is nonzero,
but the monotonized slope is zero because one one-sided difference is zero.
At the upper-y probe the corresponding stencil is mirrored. Applying the
actual limited formula reproduces all eight recorded donor velocities.
Connected directions cannot restore a slope that the limiter itself removes.
At the cut probe a PCM fallback could produce the same local value; the
regular-volume intersection and the retained-direction formula provide the
independent smooth-extremum diagnosis.

The exact analytic finite-volume sine averages predict constant-donor velocity
mismatches of 0.7417904, 0.1909149 and 0.0480754 m/s at root sizes 8, 16 and 32.
These are local interpolation estimates, **not newly computed CFD errors or
convergence orders**. The observed sequence is now: initial donor velocity
mismatch, nonzero normal flux divergence, predictor pressure reduction,
then transverse normal-flux imbalance. The 50 ns trace excludes a direct
fine transverse closure increment as the explanation of that initial growth.
It does not yet demonstrate that a proposed alternative interpolation fixes
the later 200 ns Linf behavior.

Next: isolate smooth-extremum halo reconstruction on this recorded stencil,
checking conservative child means and admissible gas/held energy before any
candidate evolution. Preserve the default pinned policy and original-reference
failures; do not classify a new policy as a translation correction.

Reproduction of the saved-data analysis:

```
python3 tools/analyze_eb_shear_stage_origin.py \
  --initial build/eb-shear-initial-20261022 \
  --stages build/eb-sync-increment-20261023 \
  --output build/eb-sync-increment-20261023/halo-origin.json
```

Mechanism provenance, project contract, Python syntax and `git diff --check`
pass. No full configuration matrix, installed binary audit, new fixed-PeleC
comparison or physical validation was run in this increment.

## Current-source repeat, 2026-10-29

A fresh build of the current dirty main (HEAD still
`0a0e00ebe9138145fbb0e0e745665d0d768df009`) and the same 16-cell/50-ns input
completed in the clean GNU/OpenMPI Release configuration. Build took
420.84 seconds and the serial diagnostic took 14.99 seconds. Exact source,
input, cache and executable hashes are in
`build/eb-shear-current-stage-20261029/record.json`.

The original strict historical audit **failed** on the default-policy table;
that failure and log remain unchanged. The default is no longer the historical
numerical implementation, and this run does not isolate individual source
changes. Its maximum-pressure-offset diagnostic differs by 5.0415e-7 Pa.
All twelve ordinary tables for the momentum, velocity and paired experimental
policies remain exactly equal to the historical parsed tables. Their original
two fine-cell stage traces also match exactly. The new probe records additional
cells, so whole-file row-count identity is not applicable.

`audit_current_scope.py` explicitly selects the original two stage probes and
five synchronization probes, reports default historical inequality separately,
and retains the original local closure gates. All 60 pre/post closure pairs
and the two-stage momentum trajectories pass. Maximum directional relative
residual is 1.618e-16; maximum raw-register residual is 1.059e-22 and maximum
trajectory residual is 2.507e-16 kg m^-2 s^-1. Both fine probes again have zero
transverse closure increments. Paired first-predictor pressure drops by
0.79089020098 Pa; subsequent transverse rates and their integrated increments
reproduce the values above. This is a fresh bounded serial diagnostic, not a
new native comparison or 200-ns convergence result.

The halo-origin analyzer now selects the two original probe IDs explicitly
and checks six unique axis/side pairs, allowing expanded trace files without
silently accepting duplicated faces. It passes on both the fresh and historical
traces (8 halo probes, 48 momentum faces each). Reproduction and hashes are in
`build/eb-shear-current-stage-20261029/verify.py` and `verified.json`.
Mechanism provenance, project contract, Python compilation and diff whitespace
checks also pass. No production numerical policy, tolerance or golden output
was changed. The isolated quartic evolution and its increased z-momentum L1
error remain separate, unadopted evidence; see
[the candidate evolution record](eb-quartic-evolution-20261023.md).
