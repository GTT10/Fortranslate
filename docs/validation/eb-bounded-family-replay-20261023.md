# Saved gas-family bounded reconstruction

All 64 recorded prepared-state families pass the bounded eight-child
reconstruction with theta=1 (no EOS contraction). The maximum normalized
volume-weighted mean residual is 2.4027323721396817e-16. The 512 returned
children have positive density, nonnegative species densities and species
sum equal to density within the reported tolerance. No evolution was run.

The new `eb_bounded_family_replay.F90` loads the original selected LiDryer
thermodynamics and the case's converged inverse-EOS setting. It rebuilds the
native axis-plane geometry from the original V4 input (16 root cells,
0.000625 m spacing and z=0.0033 m plane). Parent z=6 has four child fractions
0.44 and four 1; parent z=8 has eight 1. The prepared state arrays are read
from the saved QBOUND_SAMPLE records without analytic replacement.

The original logger did not record child volume fractions or stage EOS
initial guesses. Therefore these weights are reconstructed from the native
geometry, and the temperature guess is recovered from the parent held-energy
EOS using 900 K as its initial guess. This is an admissibility/reconstruction
experiment on recorded states, not a byte-identical replay of the original
stage call, a complete halo path, or evidence of trajectory agreement.

Artifacts are in `build/eb-bounded-family-replay-20261023/`: report.json
records source/archive/binary/input hashes, command, each family's weights,
recovered temperature, theta and conservation residual. input.dat and
output.dat preserve every prepared input and resulting child state.

```sh
python3 tools/check_eb_bounded_family_replay.py \
  --input build/eb-inflection-admission-20261023/grid16.log \
  --case build/eb-inflection-admission-20261023/grid16.dat \
  --build-dir build/eb-eos-application-20261021/isolated-openmpi/registration \
  --work-dir build/eb-bounded-family-replay-20261023
```

The compile used GNU Fortran with -O0 and bounds/FPE checks. Python syntax,
mechanism provenance, project contract and whitespace checks pass. Next is
a bounded 16-cell/50-ns candidate evolution through the real halo call, with
stage/flux/budget and all-field error comparison against saved candidates.
The isolated candidate is not a production option. C2-C5 remain open.
