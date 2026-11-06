# Initial coarse/fine stencil reconstruction, 2026-10-21

The initial spatial probe now reconstructs every adjacent fine-halo donor
and its six neighbors with the public production planner/interpolator. It
compares gas and held energy with the actual operator halos, checks the
average of all eight reconstructed children, and independently restricts
fine states with fluid-volume weights. It records EOS pressures of the
center, neighbors and interpolated state. This extends the
[initial interface audit](eb-initial-interface-20261021.md).

## Results on the saved stationary 8-grid import

All 160 open adjacent fine-halo states in the planar case reproduce the
production conserved and held-energy values exactly. The maximum scaled
child-average residual is 2.49e-16; independently volume-restricted coarse
gas and held energy agree within 2.46e-16. Thus the measured disturbance
is not a donor indexing error or failed conservative restriction.

The worst ghost is fine `(3,0,7)`, supplied from coarse `(2,8,4)`. Its
center EOS pressure is 101325 Pa, and six neighbor pressures span
101325 to 101325.136899 Pa. The interpolated EOS pressure is
101324.724778 Pa: 0.275222 Pa below the input pressure envelope despite
conserved child averages. Replacing only the packed temperature guesses
with EOS-consistent temperatures changes none of the adjacent ghost pressures.
This separates the pressure disturbance from stored-temperature guessing and
identifies pressure non-monotonicity of the conservative interpolated state.
It does not imply that linear interpolation is mistranslated from PeleC.

The old geometry-only regular control retained planar volume-weighted coarse
states. The new `regular_restricted` control recomputes every covered coarse
state from eight equal-volume children after removing the plane. Its
restriction residual is 2.46e-16, and the same maximum adjacent ghost
pressure deficit remains. Its covered-coarse EOS pressure offset is larger,
3.57980 Pa, and its coarse/fine energy-flux mismatch reaches 5207.38 W/m2.
This is a controlled geometry/restriction ablation of imported states, not a
freshly initialized reference trajectory. It confirms that removing cut cells
does not by itself eliminate the pressure disturbance, while the precise
coarse/fine flux mismatch depends on the coarse restriction.

Both momentum policies agree at this initial rest state. Uniform controls
have zero mass flux and zero reported interface flux mismatch. The script
now verifies stencil replay, child-average conservation and restriction,
except for the intentionally inconsistent old geometry-only control, which
is explicitly marked `restriction_checked=false`.

## Diagnostic candidate

The probe also tests a pressure-envelope slope factor on all eight children
of a donor. One common factor contracts the conserved increments toward the
parent, preserving their average. It uses a pressure roundoff allowance of
128 machine epsilons and bounded bisection. This is an isolated stencil
experiment, not an application option or a change to the solver fluxes.
The candidate reduces the maximum adjacent pressure offset to 2.93e-9 Pa,
but limits 48 of 128 linear stencils in the planar control. All 48 retain
less than 1% of their slopes; the worst factor is 1.81e-8, effectively PCM.
The regular/restricted control limits 64 of 320 stencils with the same
minimum factor. It is not adopted: this experiment removes the pressure
undershoot by nearly flattening resolved gradients, not by establishing an
accurate thermodynamic reconstruction. No candidate time-evolution or flux
qualification is claimed. See `candidate-summary.json`.

## Verification

Artifacts and exact identities are in `build/eb-initial-stencil-20261021/`.
GNU Fortran 13.3.0 Debug, using the existing matching OpenMPI registration
build. No compiler/MPI configuration changed. The probe was rebuilt.
Twelve final application controls completed with the replay/conservation
assertions active. Four Python tests cover periodic/open-area pairing, regular conservative
restriction, nonfinite rejection and conservation-check refusal. All pass.
Provenance, Python syntax, project contract and whitespace checks pass.
No new trajectory, native PeleC, MPI, or convergence qualification is claimed.
