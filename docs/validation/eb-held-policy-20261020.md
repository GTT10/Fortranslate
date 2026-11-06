# Spatial held-energy policy isolation, 2026-10-20

The remaining twelve nonuniform reacting EB-AMR field failures are dominated
by the spatial EOS energy choice. A diagnostic native binary using stored
UEINT for the spatial EOS agrees with production PeleF within every original
field tolerance at all four compared times. **This is an intervention in the
native algorithm, not a passing comparison against the pinned algorithm.**
The original reference capture, its acceptance limits and production PeleF
remain unchanged; the original twelve field failures remain open.

## Controlled change

Pinned `PeleCUtilities.H::pc_ctoprim` computes specific internal energy as
`(UEDEN - kinetic)/rho`, then computes temperature, pressure, sound speed and
EOS derivatives from it. `Diffusion.cpp::getMOLSrcTerm` uses those primitives
for both MOL hydro and molecular transport. PeleF instead recovers these
quantities from actual stored UEINT, preserving its cancellation contract.
The two energy representations can differ on fine cells after parent reflux.

The diagnostic copies the pinned `pc_ctoprim` function under a different
name and changes only its specific-energy expression to `UEINT/rho`.
A copied `Diffusion.cpp` calls that helper. No reference header or source
file is edited; reaction integration, energy reset, reflux and transport
coefficients are unmodified. This changes both spatial pressure and
transport temperature consistently. The previous hydro-only total-energy
substitution in Fortran did not test this complete policy alignment.

The copied native binary runs the same three 100 ns root steps, recursive
MOL2, LiDryer CVODE, Simple transport, planar EB and ratio-two fine patch,
with the same `eta2=2` setting. No PeleF run is repeated: comparison uses the
completed [covered-parent production result](eb-covered-energy-20261020.md).
All 1024 initial gas/T/UEINT cells match exactly; composite comparisons use
the same 352 coarse and 192 fine fluid cells and SI conversion as the frozen
checker. The native diagnostic is explicitly excluded from the pinned-reference
qualification path.

## Results

| Maximum absolute difference at 300 ns | Original native algorithm | Native held-EOS diagnostic |
| --- | ---: | ---: |
| Density (kg/m3) | 5.967861497e-9 | 5.019581972e-12 |
| x-momentum (kg/m2/s) | 4.764959025e-6 | 3.536982263e-9 |
| Total energy (J/m3) | 0.02020956951 | 1.781203900e-5 |
| Held energy (J/m3) | 0.02018086310 | 1.781579340e-5 |
| Temperature (K) | 1.098005948e-4 | 7.451792794e-5 |

All original field limits are satisfied by the held-policy diagnostic at
0, 100, 200 and 300 ns, without changing any limit. The density/momentum/
energy discrepancy is reduced by roughly three orders of magnitude. Stored
temperature output still has a different refresh schedule on fine cells;
its error remains below the existing `1e-3 K` limit.

Both sides retain nonnegative species and satisfy species-sum and composite
mass/element checks. Maximum mass/element drift is `3.11e-15` for native held
EOS and `5.98e-14` for PeleF. Hydrogen-consumption relative error is `2.24e-5`,
below the original `1e-3` limit. The diagnostic does not introduce a new native
energy-budget audit; production PeleF's numerical-correction budget was
checked independently in the covered-parent record.

This intervention identifies the dominant remaining discrepancy as a
numerical-policy difference, rather than evidence for another unrestricted
reset or flux change. Production continues to use actual UEINT. In particular,
matching this low-speed case is not grounds for discarding the high-cancellation
or rotating-energy behavior. Full fixed-PeleC parity, longer trajectories,
mesh/time convergence and physical validation remain separate work.

## Reproduction

Artifacts: `build/eb-held-policy-20261020/`. `run_native.py` creates and builds
the isolated native translation unit, runs the three-step diagnostic and
exports the four snapshots. `compare.py` performs the field, composition,
inventory and consumption comparisons. `comparison.json` sets
`qualification: false` and records zero `diagnostic_limit_exceedances`.
Exact compile/link/run/export commands and source/binary hashes are retained.
The production capture is
`build/eb-covered-energy-20261020/reference-runs/matched_eb_amr_z41j_1e2/`.
No tests, full matrix or application runs were repeated for this diagnostic.
