# Donor admission and momentum roundoff scale — 2026-10-23

## Completed attribution

The instrumented frozen candidate completes the 16-cell/50-ns run in
142.546 seconds. **All 24 CSV tables are byte-identical** to the prior cached
candidate, including stage and synchronization traces. Thus the extra donor
diagnostics have not changed these results. The diagnostic stores the prepared
five-cubed stencil, retained axes, smooth-screen result, EOS admission result,
contraction coefficient and ancestor-time fraction for four coarse donors.

Artifacts: `build/eb-quartic-admission-20261023/`. Build/source/executable
identity is recorded in `build-record.json`; run/input identity and the 24-table
comparison are in `record.json`. `replay_screen.py` exactly reproduces all
**48 smooth-screen decisions** from `cut.log`; the final receipt is
`admission-audit.json`. Earlier partial receipts were explicitly marked partial
and are not used as completed-run evidence.

The 18 rejections all concern donors (16,4,8) or (16,5,8), after the first
stage. The first failing component is z momentum, on an x-directed stencil.
At ancestor fraction 0.5, representative second differences are
`[0, 5.81948e-16, -5.81820e-16]`. The original component-local cancellation
allowance is only `2.06795e-30`, because the z-momentum component itself is
almost zero. At fraction 1 the values approximately double. The cut-cell-side
donors at z index 6 remain admitted. The relevant family EOS admission succeeds;
the smoothness predicate, not an EOS failure, rejects these donors.

This links the observed full-volume probe's fallback to a component-scale
problem in the experimental screen. It does not establish a general shock
detector or account for every spatial error in the candidate.

## Isolated proposed change and checks

`tools/reference_probes/eb_guarded_quartic_mod.F90` now has an optional
`momentum_roundoff_scale` argument, default false. When enabled, momentum
components use the maximum Euclidean momentum-vector magnitude along the
sampled line as their cancellation scale. The factor remains 16 machine
epsilons. Density, energy, species, curvature signs, physical-volume means,
EOS admission and common-family contraction are unchanged. No dimensional
absolute floor or fitted velocity threshold is introduced. This remains an
experimental diagnostic helper, outside production source/CMake.

Saved-stencil replay with this scale admits all 48 recorded cases. The direct
Fortran probe passes **84 assertions**, including the previous 78 and checks
for tiny transverse perturbations, conserved-unit rescaling, a resolved
transverse jump, a density jump and EOS-admissible family construction.
The resolved jumps still fail the smooth screen. Commands, GNU Fortran 13.3.0
toolchain and source/archive/binary hashes are in
`build/eb-quartic-roundoff-20261023/guarded-record.json`.

These checks justify a same-condition evolution comparison, not adoption.
The new frozen tree copy is being built in the roundoff artifact directory.
Its 50-ns run must retain default equality, conservation and stage balances,
and must evaluate all error norms, especially the previously increased
z-momentum L1 error. Regular geometry, mesh convergence, MPI/restart,
fixed-PeleC trajectory agreement and physical validation remain open.

## Independent rotation check while the mesh run proceeds

The probe now also rotates the stencil and momentum vector by 90 degrees
about z, permutes the eight children accordingly, and retains unequal z-layer
fluid-volume weights. Both original and rotated families are admitted with
theta=1; all eight transformed child states agree within 128 machine epsilons
using a componentwise scale. The case includes the near-zero transverse
perturbation that triggered the original screen problem.

The complete direct probe passes **94 assertions**. This is a discrete axis
rotation check for the isolated weighted family, not arbitrary-angle geometry
or full-solver rotational invariance. `build/eb-quartic-rotation-20261023/`
stores the previous probe source, exact new source/library/binary hashes,
compile/run commands and output. No numerical helper, solver or running mesh
binary was changed for this test. Mechanism provenance, project contract and
whitespace checks pass independently after the probe update.
